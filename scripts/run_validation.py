"""Run the lightweight repository validation commands.

This script is intentionally small. It sets ``PYTHONPATH=src`` for pytest so the
helper package can be imported without requiring packaging metadata, then runs
the notebook-section validator.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def run(command: list[str], *, env: dict[str, str] | None = None) -> None:
    print("+", " ".join(command))
    subprocess.run(command, cwd=REPO_ROOT, env=env, check=True)


def main() -> int:
    env = os.environ.copy()
    src_path = str(REPO_ROOT / "src")
    existing_pythonpath = env.get("PYTHONPATH")
    env["PYTHONPATH"] = src_path if not existing_pythonpath else f"{src_path}{os.pathsep}{existing_pythonpath}"

    try:
        run([sys.executable, "-m", "pytest", "-q"], env=env)
        run([sys.executable, "scripts/validate_notebook_sections.py"])
    except subprocess.CalledProcessError as exc:
        return exc.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
