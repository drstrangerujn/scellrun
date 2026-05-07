"""Tests for the `_pick_best_resolution` heuristic.

Closes the v0.7 ISSUES.md #005 trap and locks the v1.1.0+ tie-break
behavior — every resolution-selection branch (non-fragmented preferred,
fewest singletons fallback, no-≥2-cluster fallback) gets a guard test so
future refactors can't silently revert to the "always pick max
n_clusters" trap that picked the most fragmented resolution on real OA
data.
"""
from __future__ import annotations

from scellrun.analyze import _pick_best_resolution


def test_pick_best_resolution_empty_returns_fallback():
    """No quality table → return the fallback resolution."""
    res, reason = _pick_best_resolution(None, fallback=0.7)
    assert res == 0.7
    assert "fallback" in reason

    res, reason = _pick_best_resolution({}, fallback=0.5)
    assert res == 0.5
    assert "fallback" in reason


def test_pick_best_resolution_prefers_non_fragmented():
    """A clean resolution wins over fragmented ones."""
    quality = {
        0.3: {"n_clusters": 5, "n_singletons": 3, "largest_pct": 0.40},   # fragmented
        0.5: {"n_clusters": 7, "n_singletons": 1, "largest_pct": 0.30},   # non-fragmented
        0.8: {"n_clusters": 9, "n_singletons": 4, "largest_pct": 0.50},   # fragmented
    }
    res, reason = _pick_best_resolution(quality, fallback=0.5)
    assert res == 0.5
    assert "non-fragmented" in reason


def test_pick_best_resolution_picks_largest_n_among_clean():
    """Among multiple non-fragmented resolutions, pick the one with most clusters."""
    quality = {
        0.3: {"n_clusters": 5, "n_singletons": 0, "largest_pct": 0.40},
        0.5: {"n_clusters": 8, "n_singletons": 1, "largest_pct": 0.25},
        0.8: {"n_clusters": 12, "n_singletons": 1, "largest_pct": 0.20},
    }
    res, reason = _pick_best_resolution(quality)
    assert res == 0.8
    assert "largest n_clusters" in reason


def test_pick_best_resolution_ties_break_to_smaller_resolution():
    """When non-fragmented candidates have the same n_clusters, the smaller
    resolution number wins (more conservative clustering)."""
    quality = {
        0.3: {"n_clusters": 7, "n_singletons": 0, "largest_pct": 0.30},
        0.5: {"n_clusters": 7, "n_singletons": 0, "largest_pct": 0.28},
        0.8: {"n_clusters": 7, "n_singletons": 1, "largest_pct": 0.26},
    }
    res, _ = _pick_best_resolution(quality)
    assert res == 0.3


def test_pick_best_resolution_all_fragmented_picks_fewest_singletons():
    """v0.7 ISSUES.md #005: when every resolution has >1 singleton, the picker
    must NOT fall through to max n_clusters (which on real OA data picked the
    most fragmented resolution). Instead, prefer fewest singletons → smallest
    largest_pct."""
    quality = {
        0.3: {"n_clusters": 13, "n_singletons": 4, "largest_pct": 0.31},
        0.5: {"n_clusters": 16, "n_singletons": 3, "largest_pct": 0.27},  # fewest singletons → win
        0.8: {"n_clusters": 22, "n_singletons": 8, "largest_pct": 0.18},
    }
    res, reason = _pick_best_resolution(quality)
    assert res == 0.5
    assert "fewest singletons" in reason
    assert "fragmented" in reason


def test_pick_best_resolution_all_fragmented_ties_singletons_picks_balanced():
    """All-fragmented + same number of singletons → smallest largest_pct
    (most balanced, least one-cluster-dominated)."""
    quality = {
        0.3: {"n_clusters": 8, "n_singletons": 3, "largest_pct": 0.55},   # very dominated
        0.5: {"n_clusters": 10, "n_singletons": 3, "largest_pct": 0.32},  # most balanced → win
        0.8: {"n_clusters": 14, "n_singletons": 3, "largest_pct": 0.40},
    }
    res, _ = _pick_best_resolution(quality)
    assert res == 0.5


def test_pick_best_resolution_all_fragmented_full_tie_picks_smallest_res():
    """All-fragmented + tied singletons + tied largest_pct → smallest resolution."""
    quality = {
        0.3: {"n_clusters": 6, "n_singletons": 2, "largest_pct": 0.30},
        0.5: {"n_clusters": 9, "n_singletons": 2, "largest_pct": 0.30},
        0.8: {"n_clusters": 13, "n_singletons": 2, "largest_pct": 0.30},
    }
    res, _ = _pick_best_resolution(quality)
    assert res == 0.3


def test_pick_best_resolution_no_resolution_has_two_plus_clusters():
    """Tiny / failed runs may produce only single-cluster resolutions. Fall
    back instead of returning a meaningless 1-cluster pick."""
    quality = {
        0.1: {"n_clusters": 1, "n_singletons": 0, "largest_pct": 1.0},
        0.3: {"n_clusters": 1, "n_singletons": 0, "largest_pct": 1.0},
    }
    res, reason = _pick_best_resolution(quality, fallback=0.5)
    assert res == 0.5
    assert "fallback" in reason
    assert "≥2 clusters" in reason


def test_pick_best_resolution_v07_bml1_trap():
    """Regression guard for the exact ISSUES.md #005 trap on real OA data:
    every resolution had >1 singleton, max-n_clusters fell to the most
    fragmented one. v1.0.1+ behavior must NOT pick the highest-singleton row.
    """
    # Numbers shape-matching the v0.7 BML_1 dogfood quality table.
    quality = {
        0.05: {"n_clusters": 4, "n_singletons": 2, "largest_pct": 0.42},
        0.1: {"n_clusters": 7, "n_singletons": 3, "largest_pct": 0.36},
        0.3: {"n_clusters": 13, "n_singletons": 2, "largest_pct": 0.315},  # win
        0.5: {"n_clusters": 16, "n_singletons": 4, "largest_pct": 0.28},
        0.8: {"n_clusters": 22, "n_singletons": 8, "largest_pct": 0.22},
    }
    res, reason = _pick_best_resolution(quality)
    # Win goes to res=0.3 (fewest singletons + smaller largest_pct than 0.05).
    assert res == 0.3
    # And specifically NOT res=0.8 (max n_clusters but most fragmented).
    assert res != 0.8, "v0.7 #005 trap regressed: picked the most fragmented resolution"
