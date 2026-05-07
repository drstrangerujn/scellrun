# scellrun

[![CI](https://github.com/drstrangerujn/scellrun/actions/workflows/ci.yml/badge.svg)](https://github.com/drstrangerujn/scellrun/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![PyPI](https://img.shields.io/pypi/v/scellrun.svg)](https://pypi.org/project/scellrun/)

> 一个带决策日志的scRNA-seq分析CLI——两个分析师、两个agent跑同一份数据，得到同一个答案。

[English README](README.md)

---

## 30秒了解

同一份scRNA数据交给两个分析师重新分析，cluster数量大概会差20%（[Lähnemann等，2020 / PMC9122178](https://pmc.ncbi.nlm.nih.gov/articles/PMC9122178/)）。差的不是科学，是几十个没人写下来的细节决定：mt%上限多少、HVG选几个、用什么整合方法、cluster分辨率怎么定、注释panel选哪个。半年后回头看，谁也不记得当时改了什么。

scellrun是套在[scanpy](https://scanpy.readthedocs.io)外面的一层CLI，每个这种细节决定都会被记到一个可grep的文件里，旁边附带一句话理由。同一份数据，同一套默认值，每次跑出来一致。审稿人问"为什么mt%是20"，你直接打开`00_decisions.jsonl`第14行念给他听。

scellrun不是工作流管理器（集群级编排请用[nf-core](https://nf-co.re/scrnaseq)），也不替代scanpy——底层调用的还是scanpy，只是把参数定死、把决定记下来。

## 快速上手

```bash
conda create -n scellrun python=3.11 -y
conda activate scellrun
pip install scellrun

scellrun analyze data.h5ad --tissue "OA cartilage"
# 输出：./scellrun_out/run-<时间戳>/05_report/index.html
```

只有cellranger输出没有h5ad？

```bash
scellrun scrna convert path/to/cellranger_outs -o data.h5ad
scellrun analyze data.h5ad --tissue "OA cartilage"
```

加`--lang zh`生成中文报告。加`--profile joint-disease`处理软骨/滑膜/软骨下骨样本。完整流程见[`docs/quickstart.md`](docs/quickstart.md)。

## 你能直接拿到什么

- **一键五阶段**：QC → 整合（Harmony）→ markers → annotate → 报告。一条命令一份能直接发邮件的HTML
- **决策日志** `<run>/00_decisions.jsonl`——每个非平凡决定都附一句话理由。可grep；标auto/user/ai三种来源
- **5个组织profile**——`default`、`joint-disease`（Fan 2024软骨亚型panel）、`tumor`、`brain`、`kidney`。每个一个Python文件，欢迎贡献
- **自检机制**——每个阶段检查异常（panel不匹配、所有分辨率都碎、单样本却开了批次校正等），并给出最便宜的修复建议。`--auto-fix`自动应用
- **审阅回路**——`scellrun review <run>`起一个本地Flask小站，手动改cluster标签、调阈值、写备注；`analyze --apply-overrides <json>`带着你的修改重跑，记成`source="user"`决策行
- **PDF导出**——`scellrun export <run> --format pdf`，发文用

## 谁该用这个

- **跑分析的LLM agent**（Claude Code、Hermes、Codex、OpenClaw）。这是主要用户。把[`skills/scellrun/SKILL.md`](skills/scellrun/SKILL.md)放到agent的skills目录，agent就知道哪个命令对应哪种用户意图、怎么读决策日志、自检触发时该说什么
- **临床生信团队**。希望每个项目报告长得一样——同样的QC布局、同样的决策表、同样的可追溯链路，跨样本、跨学生、跨轮转都不变样
- **审稿人**：问"为什么mt%是20"。打开`00_decisions.jsonl`第14行就有答案

[`docs/agent-demo.md`](docs/agent-demo.md)是Claude Code agent在真实OA软骨数据上端到端跑完scellrun的逐字记录，包括用户问"为啥选res=0.3"时agent直接引用决策日志回答。

## 决策日志长什么样

```jsonl
{"stage":"qc","key":"max_pct_mt","value":20.0,"default":20.0,"source":"auto",
 "rationale":"mt%上限20% — 关节组织本身应激高，教科书的10%会把真正的软骨细胞过滤掉（PI 2024-2026队列经验，AIO PM=20）"}
{"stage":"analyze","key":"method_downgrade","value":"none","default":"harmony","source":"auto",
 "rationale":"obs里没有sample/batch列——单样本输入；自动把--method从harmony降到none"}
{"stage":"analyze","key":"chosen_resolution_for_annotate","value":0.3,"source":"auto",
 "rationale":"最少singleton + 最均衡（所有分辨率都碎）——选res=0.3：n_clusters=13，最大cluster占31.5%，最小0.2%，singleton 2个"}
{"stage":"analyze","key":"annotate.auto_panel","value":"celltype_broad","source":"auto",
 "rationale":"切到celltype_broad：chondrocyte命中=2，broad命中=9；保留chondrocyte panel需要chondrocyte命中至少是broad的1.5倍"}
```

`source`三种：`auto`（内置启发式）、`user`（CLI或review手动覆盖）、`ai`（LLM调用）。`attempt_id`把同一次调用的所有行串起来；`fix_payload`只在自检的`*.suggest`行里有，记录orchestrator能机械应用的修复内容。完整schema见[`skills/scellrun/SKILL.md`](skills/scellrun/SKILL.md)。

v1.3.2版本起，`chosen_resolution_for_annotate`的理由、panel自动选择的依据、自检触发的具体含义都会直接显示在HTML报告的"At a glance"区——不用再去翻jsonl了。

## Profile

profile是把"某个组织/疾病的工作经验"沉淀成一个Python文件，里面包含默认阈值和marker panel。

| profile | mt% | hb% | panels | 备注 |
|---|---|---|---|---|
| `default` | 20% | — | — | 新鲜组织10x v3基线（OARSI上限） |
| `joint-disease` | 20% | 紧 | Fan 2024的11个软骨亚型 + 15组broad | 已冷验证；免疫细胞多时自动切broad |
| `tumor` | 20% | — | TISCH/Sun 2021泛癌TME（仅broad） | 暂未冷验证 |
| `brain` | 10% | — | Tasic/Hodge皮层-海马（仅broad） | 暂未冷验证 |
| `kidney` | 15% | — | KPMP/Stewart 2019肾单位+免疫（仅broad） | 暂未冷验证 |

```bash
scellrun profiles list
scellrun profiles show joint-disease   # 看阈值+panel
```

你的组织或疾病有跟默认值不一样的工作惯例？写一个profile贡献回来——`src/scellrun/profiles/`下加一个Python文件就够。详见[`docs/contributing.md`](docs/contributing.md)。

## 项目状态

**v1.3已冻结surface——只做scRNA。** v1.x进入维护模式：只修bug、只加scRNA profile，不增加新公开命令。bulk RNA-seq、代谢组学、蛋白组学这些扩展推到将来的v2.0，见[`ROADMAP.md`](ROADMAP.md)。

CLI命令（`qc`/`integrate`/`markers`/`annotate`/`analyze`/`review`/`export`/`profiles`）在整个v1.x期内不变。

## 安装渠道

- **PyPI**：`pip install scellrun`（[pypi.org/project/scellrun](https://pypi.org/project/scellrun/)）
- **ClawHub**（agent skill）：`clawhub install scellrun`（[clawhub.ai/skills/scellrun](https://clawhub.ai/skills/scellrun)）
- **Docker**：`docker pull ghcr.io/drstrangerujn/scellrun:latest`（v1.0起）

## 协议

MIT——见[`LICENSE`](LICENSE)。

## 致谢

默认值来自刘老师组内R AIO流程和临床生信团队在OARSI/骨肌研究上的工作惯例。`joint-disease` profile里的软骨亚型panel来自Fan 2024。
