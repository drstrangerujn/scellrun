"""
Schema-level tests for README.md (v1.0.2).

These check structure and presence of differentiation anchors, not
phrasing — README content is allowed to evolve as long as the
operational hooks (decision log, joint-disease panel anchor, the field
quote anchoring "why this and not vanilla agent") stay grep-able.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
README_PATH = REPO_ROOT / "README.md"


@pytest.fixture(scope="module")
def readme_text() -> str:
    assert README_PATH.exists(), f"README.md missing at {README_PATH}"
    return README_PATH.read_text(encoding="utf-8")


@pytest.mark.parametrize(
    "needle",
    [
        # The decision log is the differentiation centerpiece — must be referenced.
        "00_decisions.jsonl",
        # Fan 2024 anchors the joint-disease panel claim.
        "Fan 2024",
        # The 20% reanalysis-divergence figure is the field quote that
        # justifies why this exists at all.
        "20%",
    ],
)
def test_readme_mentions_required_anchors(readme_text: str, needle: str) -> None:
    """The body must reference each differentiation anchor."""
    assert needle in readme_text, f"README.md must reference {needle!r}"


def test_readme_has_fenced_code_block(readme_text: str) -> None:
    """At least one fenced code block — the install / quick-start one-liner."""
    fences = re.findall(r"^```", readme_text, flags=re.MULTILINE)
    assert len(fences) >= 2, "README.md must contain at least one fenced code block"


@pytest.mark.parametrize(
    "heading",
    [
        "Why this exists",
        "Who this is for",
    ],
)
def test_readme_has_required_h2_sections(readme_text: str, heading: str) -> None:
    """README.md must carry the differentiation framing headings."""
    pattern = rf"^##\s+{re.escape(heading)}\s*$"
    assert re.search(pattern, readme_text, flags=re.MULTILINE) is not None, (
        f"README.md must have an H2 section titled {heading!r}"
    )
