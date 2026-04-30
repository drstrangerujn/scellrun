# scellrun skills

> SKILL.md v2 (2026-04-30): rewritten for scellrun v0.8+ — covers the
> `analyze` one-shot command, the `00_decisions.jsonl` log, self-check
> findings, profile selection by tissue, and an end-to-end agent
> dialogue example. Re-symlink if you installed an earlier copy.

`skills/scellrun/SKILL.md` is a single canonical instruction document that
teaches an LLM agent how to invoke scellrun. It uses Claude Code's frontmatter
format, but the body is agent-agnostic.

## Install for your agent

### Claude Code

```bash
mkdir -p ~/.claude/skills
ln -s "$(pwd)/skills/scellrun" ~/.claude/skills/scellrun
```

The skill auto-loads via the Skill tool. Test with: `/scellrun` or simply
ask the agent to "run QC on this .h5ad" — it should reach for scellrun.

### Hermes

```bash
mkdir -p ~/.hermes/skills
ln -s "$(pwd)/skills/scellrun/SKILL.md" ~/.hermes/skills/scellrun.md
```

Hermes loads markdown skills from its skills directory at startup; restart
the Hermes process to pick up new files.

### Codex CLI

```bash
mkdir -p ~/.codex/skills
ln -s "$(pwd)/skills/scellrun/SKILL.md" ~/.codex/skills/scellrun.md
```

Or, if you maintain a project-local `.codex/` directory, drop the
symlink there instead.

### Generic (any markdown-aware agent)

`skills/scellrun/SKILL.md` is just structured markdown with a YAML
frontmatter header. Most agent harnesses can ingest it as a system
prompt or skill file. If yours can't, paste the body into the system
message.

## Authoring a profile (not a skill change)

Want different defaults for tumor scRNA, kidney biopsy, or your own
tissue? Don't add an agent skill — add a scellrun **profile** instead:

```bash
# 1. Add src/scellrun/profiles/<your_profile>.py
# 2. Register the dash-name in src/scellrun/profiles/__init__.py:_REGISTRY
# 3. Open a PR
```

Profiles are the contributor-facing extensibility point; skills are
about teaching agents to *use* scellrun, not extending it.
