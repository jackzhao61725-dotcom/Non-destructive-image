"""Extract stored notebook outputs into a lightweight regression baseline.

This script does not execute the notebook. It records the outputs already stored
inside the reference notebook, including stream text and hashes of rich display
payloads, so future changes can detect accidental notebook-output drift.
"""

from __future__ import annotations

import base64
import hashlib
import json
from pathlib import Path
from typing import Any

NOTEBOOK_PATH = Path("1 calculations revised 2  multishot  6  extended.ipynb")
BASELINE_PATH = Path("regression/baseline/notebook_outputs.json")


def _normalise_text(value: Any) -> str:
    if isinstance(value, list):
        return "".join(str(part) for part in value)
    return str(value)


def _hash_payload(value: Any) -> dict[str, Any]:
    if isinstance(value, list):
        raw = "".join(value).encode()
    elif isinstance(value, str):
        try:
            raw = base64.b64decode(value, validate=True)
        except Exception:
            raw = value.encode()
    else:
        raw = json.dumps(value, sort_keys=True).encode()
    return {"sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}


def build_baseline() -> dict[str, Any]:
    notebook = json.loads(NOTEBOOK_PATH.read_text())
    cells: list[dict[str, Any]] = []
    for cell_index, cell in enumerate(notebook["cells"]):
        if cell.get("cell_type") != "code" or not cell.get("outputs"):
            continue
        cell_record: dict[str, Any] = {
            "cell": cell_index,
            "execution_count": cell.get("execution_count"),
            "first_source_line": next((line.strip() for line in cell.get("source", []) if line.strip()), ""),
            "outputs": [],
        }
        for output in cell["outputs"]:
            output_record: dict[str, Any] = {"output_type": output.get("output_type")}
            if "name" in output:
                output_record["name"] = output["name"]
            if "text" in output:
                output_record["text"] = _normalise_text(output["text"])
            if "ename" in output:
                output_record["ename"] = output["ename"]
                output_record["evalue"] = output.get("evalue")
            if "data" in output:
                output_record["data"] = {
                    mime_type: _hash_payload(payload)
                    for mime_type, payload in sorted(output["data"].items())
                }
            cell_record["outputs"].append(output_record)
        cells.append(cell_record)
    return {
        "notebook": str(NOTEBOOK_PATH),
        "baseline_kind": "stored_notebook_outputs",
        "note": "Extracted from outputs already stored in the reference notebook; notebook was not executed.",
        "code_cells_with_outputs": len(cells),
        "cells": cells,
    }


def main() -> None:
    BASELINE_PATH.parent.mkdir(parents=True, exist_ok=True)
    BASELINE_PATH.write_text(json.dumps(build_baseline(), indent=2, ensure_ascii=False) + "\n")
    print(f"Wrote {BASELINE_PATH}")


if __name__ == "__main__":
    main()
