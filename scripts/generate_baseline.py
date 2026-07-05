"""Prepare/generate the future scientific numerical baseline.

This milestone intentionally does not generate numerical data in environments
without the scientific Python stack. The script checks prerequisites and prints a
clear next action instead of failing with import tracebacks.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

REQUIRED_MODULES = ["numpy", "scipy", "matplotlib", "nbformat", "nbclient"]
BASELINE_DIR = Path("regression/baseline")
PLANNED_DIRS = [BASELINE_DIR / "atomic", BASELINE_DIR / "imaging", BASELINE_DIR / "multishot"]


def missing_modules() -> list[str]:
    return [module for module in REQUIRED_MODULES if importlib.util.find_spec(module) is None]


def ensure_layout() -> None:
    for directory in PLANNED_DIRS:
        directory.mkdir(parents=True, exist_ok=True)


def main() -> int:
    ensure_layout()
    missing = missing_modules()
    if missing:
        print("Scientific baseline generation was not run.")
        print("Missing required dependencies: " + ", ".join(missing))
        print("Install the scientific Python stack, then rerun:")
        print("  python3 scripts/generate_baseline.py")
        print("See docs/baseline_specification.md for the planned output files.")
        return 2

    print("All required dependencies are available.")
    print("Baseline generation is intentionally gated until reviewed.")
    print("Next implementation step: execute the original notebook and save the arrays specified in docs/baseline_specification.md.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
