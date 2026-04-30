"""
Schema-level tests for skills/scellrun/SKILL.md (v1.0).

These check structure and presence of required topics, not phrasing —
the SKILL.md content is allowed to evolve as long as agents can still
find the operational hooks.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SKILL_PATH = REPO_ROOT / "skills" / "scellrun" / "SKILL.md"


@pytest.fixture(scope="module")
def skill_text() -> str:
    assert SKILL_PATH.exists(), f"SKILL.md missing at {SKILL_PATH}"
    return SKILL_PATH.read_text(encoding="utf-8")


def test_skill_has_yaml_frontmatter(skill_text: str) -> None:
    """SKILL.md must open with a YAML frontmatter block."""
    assert skill_text.startswith("---\n"), "SKILL.md must start with a YAML frontmatter '---' line"
    # Frontmatter ends at the next '---' on its own line.
    closing = skill_text.find("\n---\n", 4)
    assert closing > 0, "SKILL.md frontmatter must close with a '---' line"


def test_frontmatter_name_and_description(skill_text: str) -> None:
    """Frontmatter must declare name=scellrun and a description >100 chars."""
    closing = skill_text.find("\n---\n", 4)
    frontmatter = skill_text[4:closing]

    name_match = re.search(r"^name:\s*(\S+)\s*$", frontmatter, flags=re.MULTILINE)
    assert name_match is not None, "frontmatter missing 'name:' line"
    assert name_match.group(1) == "scellrun", f"name must be 'scellrun', got {name_match.group(1)!r}"

    # description: may span lines but the YAML scalar value lives on the first line for our doc.
    desc_match = re.search(r"^description:\s*(.+?)(?:\n[a-zA-Z_-]+:|\Z)", frontmatter, flags=re.MULTILINE | re.DOTALL)
    assert desc_match is not None, "frontmatter missing 'description:' line"
    description = desc_match.group(1).strip()
    assert len(description) > 100, f"description too short ({len(description)} chars); needs to actually instruct an agent"


@pytest.mark.parametrize(
    "needle",
    [
        "scellrun analyze",
        "00_decisions.jsonl",
        "self_check",
        "joint-disease",
        "--ai",
        "--auto-fix",
        # v0.9.1 decision-log schema fields that an agent must know about.
        "schema_version",
        "attempt_id",
        "fix_payload",
        # v0.9.1 raised the QC self-check trigger ceiling from 30% to 60%.
        "60%",
        # v0.9.1 preserves failed first-pass artifacts when --auto-fix retries.
        ".failed-1",
        # v0.9.1 single-sample auto-degrade for the integrate stage.
        "harmony→none",
        # v1.0 honesty: a "Verification status" section pinning to v1.0.0.
        "Verification status",
    ],
)
def test_body_mentions_required_topics(skill_text: str, needle: str) -> None:
    """The body must reference each topic an agent reading cold needs to find."""
    assert needle in skill_text, f"SKILL.md must reference {needle!r} (operational hook for agents)"


def test_body_has_fenced_code_block(skill_text: str) -> None:
    """At least one fenced code block — install one-liner, dialogue, or sample log."""
    # Match opening fences only ("```lang" or "```") at line start.
    fences = re.findall(r"^```", skill_text, flags=re.MULTILINE)
    # Each block has an opening + closing fence, so we need >= 2 lines.
    assert len(fences) >= 2, "SKILL.md must contain at least one fenced code block"
