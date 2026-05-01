# scellrun skills

> v3 (2026-05-01) — tracks scellrun **v1.1.2**. Agents check the
> frontmatter's `tested_against_version`, run `scellrun --version`,
> warn or refuse if out of range.

## Install

Don't install this yourself — **hand it to your LLM agent** and ask
the agent to install. Agents know where their own skills directory
lives.

If the agent asks where, the typical paths are:

```
Claude Code / openclaw / forks  →  ~/.claude/skills/scellrun/SKILL.md
Hermes                          →  ~/.hermes/skills/scellrun.md
Codex CLI                       →  ~/.codex/skills/scellrun.md
anything else                   →  paste the SKILL.md body into the
                                   agent's system prompt
```

The actual file is `skills/scellrun/SKILL.md`. It's plain markdown
with a YAML frontmatter header — every modern agent harness can
ingest it.

## Adding a profile (different from a skill)

Want different defaults for tumor / brain / kidney scRNA? Don't
write a skill — add a scellrun **profile**:

```
src/scellrun/profiles/<your_profile>.py    # one Python file
src/scellrun/profiles/__init__.py:_REGISTRY # add the dash-name
```

Open a PR. Profiles are how working practice gets shared.
