"""
Profiles package: each profile is a small module that exports per-assay
threshold dataclasses (and, for domain profiles, marker panels and other
opinionated structures).

Adding a profile = one PR adding a file under this package.

Profile names use dashes externally (`--profile joint-disease`) but the
underlying module uses underscores (`joint_disease`); resolution handles
the swap.
"""
from __future__ import annotations

from importlib import import_module
from typing import Any

_REGISTRY: list[str] = ["default", "joint-disease"]


def _module_name(profile: str) -> str:
    return profile.replace("-", "_")


def load(name: str) -> Any:
    """Load a profile module by external name. Raises ModuleNotFoundError if missing."""
    return import_module(f"scellrun.profiles.{_module_name(name)}")


def list_profiles() -> list[str]:
    return list(_REGISTRY)
