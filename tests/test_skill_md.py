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
        # v1.0 honesty: a "Verification status" section pinning to a version.
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


# v1.1.1 additions: version-sync frontmatter, When-NOT-to-use section,
# expanded Environment hygiene subsections.


@pytest.mark.parametrize(
    "key",
    [
        "min_scellrun_version",
        "tested_against_version",
        "schema_version",
    ],
)
def test_frontmatter_has_version_sync_keys(skill_text: str, key: str) -> None:
    """v1.1.1 frontmatter must declare version-sync keys so the agent can compare against the installed CLI."""
    closing = skill_text.find("\n---\n", 4)
    assert closing > 0, "SKILL.md frontmatter must close with a '---' line"
    frontmatter = skill_text[4:closing]
    pattern = rf"^{re.escape(key)}\s*:\s*\S+"
    assert re.search(pattern, frontmatter, flags=re.MULTILINE) is not None, (
        f"frontmatter must declare {key!r} (v1.1.1 version-sync requirement)"
    )


def test_has_when_not_to_use_section(skill_text: str) -> None:
    """v1.1.1 must carry a 'When NOT to use scellrun' H2 section so the agent knows the tool's edges."""
    assert re.search(
        r"^##\s+When NOT to use scellrun\s*$", skill_text, flags=re.MULTILINE
    ) is not None, "SKILL.md must contain an H2 section titled 'When NOT to use scellrun'"


@pytest.mark.parametrize(
    "subsection",
    [
        "Remote-server execution",
        "Multi-sample analysis",
        "Auto-fix retry hygiene",
        "Reporting back to the user",
    ],
)
def test_environment_hygiene_subsections(skill_text: str, subsection: str) -> None:
    """v1.1.1 expands Environment hygiene with these H3 subsections; each must be present."""
    pattern = rf"^###\s+{re.escape(subsection)}"
    assert re.search(pattern, skill_text, flags=re.MULTILINE) is not None, (
        f"SKILL.md must contain an H3 subsection starting with {subsection!r} under Environment hygiene"
    )


# v1.1.2 additions: top-of-doc version-check ritual, expanded
# When-NOT-to-use boundary cases, four more Environment-hygiene
# subsections, profile-generic glossary structure.


def test_has_version_compatibility_check_section(skill_text: str) -> None:
    """v1.1.2 must carry a 'Version compatibility check' H2 placed near the top so the agent runs it BEFORE any scellrun command."""
    assert re.search(
        r"^##\s+Version compatibility check",
        skill_text,
        flags=re.MULTILINE,
    ) is not None, "SKILL.md must contain an H2 section starting with 'Version compatibility check'"


@pytest.mark.parametrize(
    "needle",
    [
        "multimodal",
        "already integrated",
        "missing raw counts",
    ],
)
def test_when_not_to_use_boundary_cases(skill_text: str, needle: str) -> None:
    """v1.1.2 expands When-NOT-to-use with explicit multimodal / already-integrated / raw-count boundary cases."""
    # Locate the When-NOT-to-use section and the next H2 (its end).
    start_match = re.search(
        r"^##\s+When NOT to use scellrun\s*$", skill_text, flags=re.MULTILINE
    )
    assert start_match is not None, "When-NOT-to-use section missing (precondition)"
    section_start = start_match.end()
    next_h2 = re.search(r"^##\s+", skill_text[section_start:], flags=re.MULTILINE)
    section_end = section_start + next_h2.start() if next_h2 else len(skill_text)
    section_text = skill_text[section_start:section_end].lower()
    assert needle.lower() in section_text, (
        f"'When NOT to use scellrun' section must mention {needle!r}"
    )


@pytest.mark.parametrize(
    "subsection",
    [
        "Long-running remote jobs",
        "Disk space and run-dir naming",
        "Resume / inspect a half-finished run",
        "Secret handling",
    ],
)
def test_environment_hygiene_v112_subsections(skill_text: str, subsection: str) -> None:
    """v1.1.2 adds four operational H3 subsections under Environment hygiene; each must be present."""
    pattern = rf"^###\s+{re.escape(subsection)}"
    assert re.search(pattern, skill_text, flags=re.MULTILINE) is not None, (
        f"SKILL.md must contain an H3 subsection starting with {subsection!r}"
    )


def test_profile_specific_glossaries_section(skill_text: str) -> None:
    """v1.1.2 renames the chondrocyte glossary H2 to a profile-generic 'Profile-specific glossaries' container so future profiles slot in cleanly."""
    assert re.search(
        r"^##\s+Profile-specific glossaries\s*$",
        skill_text,
        flags=re.MULTILINE,
    ) is not None, "SKILL.md must contain an H2 section titled 'Profile-specific glossaries'"
    # And the chondrocyte panel must now sit under it as an H3.
    assert re.search(
        r"^###\s+`chondrocyte_markers`",
        skill_text,
        flags=re.MULTILINE,
    ) is not None, "Profile-specific glossaries must contain an H3 subsection for `chondrocyte_markers`"
