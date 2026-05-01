# scellrun skills

> SKILL.md v3 (2026-05-01): tracks scellrun **v1.1.2**. Covers the
> `analyze` one-shot command, the `00_decisions.jsonl` log, self-check
> findings, profile selection by tissue, when-NOT-to-use boundaries,
> the version compatibility ritual, long-running-job hygiene, and an
> end-to-end agent dialogue. Re-symlink if you installed an earlier
> copy. The skill's frontmatter declares `tested_against_version`; an
> agent loading the skill should read that field, run `scellrun --version`,
> and refuse / warn if the installed version is out of range.

`skills/scellrun/SKILL.md` is a single canonical instruction document that
teaches an LLM agent how to invoke scellrun. It uses Claude Code-style
YAML frontmatter, but the body is agent-agnostic — every agent harness
listed below loads it as plain markdown.

## Install for your agent

### Claude Code

```bash
mkdir -p ~/.claude/skills
ln -s "$(pwd)/skills/scellrun" ~/.claude/skills/scellrun
```

The skill auto-loads via the Skill tool. Test with: `/scellrun` or simply
ask the agent to "run QC on this .h5ad" — it should reach for scellrun.

### openclaw / Claude Code fork

If you run a Claude Code fork (openclaw, Chorus, or similar), the skill
loader path is the same as upstream Claude Code:

```bash
mkdir -p ~/.claude/skills      # or whatever your fork's skills dir is
ln -s "$(pwd)/skills/scellrun" ~/.claude/skills/scellrun
```

If your fork uses a different skills root (`~/.openclaw/skills`,
`~/.chorus/skills`, etc.), swap that path in. The frontmatter format and
markdown body are identical to upstream — no fork-specific changes
needed in the skill file itself.

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
message. The frontmatter keys `min_scellrun_version` and
`tested_against_version` are plain strings — any harness can read them
even without YAML support.

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
