"""Verify future scientific baseline layout and metadata presence."""

from __future__ import annotations

from pathlib import Path

BASELINE_DIR = Path("regression/baseline")
EXPECTED_DIRS = [BASELINE_DIR / "atomic", BASELINE_DIR / "imaging", BASELINE_DIR / "multishot"]


def main() -> int:
    missing_dirs = [str(path) for path in EXPECTED_DIRS if not path.exists()]
    if missing_dirs:
        print("Baseline layout is incomplete.")
        print("Missing directories: " + ", ".join(missing_dirs))
        return 1
    metadata = BASELINE_DIR / "metadata.json"
    if not metadata.exists():
        print("Baseline layout exists, but generated scientific metadata is not present yet.")
        print("This is expected before the scientific baseline generation milestone.")
        return 0
    print("Baseline metadata found. Full metadata validation will be implemented with baseline generation.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
