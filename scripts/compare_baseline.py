"""Compare future scientific baseline outputs.

Currently this script validates that the planned baseline layout exists and that
required scientific dependencies are available before numerical comparison is
implemented in the baseline-generation milestone.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

REQUIRED_FOR_ARRAY_COMPARE = ["numpy"]
BASELINE_DIR = Path("regression/baseline")


def main() -> int:
    missing = [module for module in REQUIRED_FOR_ARRAY_COMPARE if importlib.util.find_spec(module) is None]
    if missing:
        print("Scientific baseline comparison was not run.")
        print("Missing required dependencies: " + ", ".join(missing))
        return 2
    if not BASELINE_DIR.exists():
        print(f"Baseline directory is missing: {BASELINE_DIR}")
        return 1
    print("Baseline comparison infrastructure is present.")
    print("No array comparison was performed because scientific baseline arrays have not been generated yet.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
