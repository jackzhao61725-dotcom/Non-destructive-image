from __future__ import annotations

import json

from scripts.extract_notebook_output_baseline import build_baseline


def test_stored_notebook_outputs_match_baseline() -> None:
    baseline = json.loads(open("regression/baseline/notebook_outputs.json", encoding="utf-8").read())
    current = build_baseline()

    assert current == baseline
