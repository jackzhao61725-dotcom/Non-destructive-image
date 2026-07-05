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
    # Include both the helper source tree and repository root. The root is
    # needed so tests can import validation helpers from the top-level
    # ``scripts/`` directory when pytest is launched as a console script.
    pythonpath_entries = [str(REPO_ROOT / "src"), str(REPO_ROOT)]
    existing_pythonpath = env.get("PYTHONPATH")
    if existing_pythonpath:
        pythonpath_entries.append(existing_pythonpath)
    env["PYTHONPATH"] = os.pathsep.join(pythonpath_entries)

    try:
        run([sys.executable, "-m", "pytest", "-q"], env=env)
        run([sys.executable, "scripts/validate_notebook_sections.py"])
    except subprocess.CalledProcessError as exc:
        return exc.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
