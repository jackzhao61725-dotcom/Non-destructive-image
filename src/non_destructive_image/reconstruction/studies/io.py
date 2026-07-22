"""Small provenance and tabular I/O helpers for reconstruction studies."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

import numpy as np


def load_json(path: Path) -> dict[str, Any]:
    """Load a JSON object from ``path``."""

    with path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"expected a JSON object in {path}")
    return payload


def load_rows(path: Path) -> list[dict[str, str]]:
    """Load a CSV table without guessing numerical types."""

    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _json_value(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    return value


def write_json(path: Path, payload: dict[str, Any]) -> None:
    """Write JSON with deterministic key ordering and NumPy conversion."""

    with path.open("w", encoding="utf-8") as handle:
        json.dump(_json_value(payload), handle, indent=2, sort_keys=True)
        handle.write("\n")


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    """Write a non-empty row table with a stable declared column order."""

    if not rows:
        raise ValueError(f"cannot write an empty table to {path}")
    fieldnames = list(rows[0])
    if any(list(row) != fieldnames for row in rows):
        raise ValueError("every output row must use the same ordered fields")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def git_value(repository_root: Path, *arguments: str) -> str:
    """Return one Git value, using ``unknown`` when Git is unavailable."""

    try:
        return subprocess.check_output(
            ["git", *arguments],
            cwd=repository_root,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def file_sha256(path: Path) -> str:
    """Return the SHA-256 digest of a study input file."""

    return hashlib.sha256(path.read_bytes()).hexdigest()
