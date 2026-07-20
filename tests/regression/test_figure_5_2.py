from __future__ import annotations

import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "results" / "dissertation_plots_v1" / "figure_5_2"


def _values() -> dict:
    return json.loads((OUTPUT / "figure_5_2_values.json").read_text(encoding="utf-8"))


def test_figure_5_2_reference_contract() -> None:
    values = _values()
    reference = values["canonical_gate"]["reference"]

    assert values["canonical_gate"]["passed"] is True
    assert reference["depletion_limited_frames"] == 10
    assert reference["usable_frames"] == {
        "Faraday dark-field": 0,
        "Faraday dual-port": 10,
    }
    assert reference["initial_snr"] == pytest.approx(
        {
            "Faraday dark-field": 6.227105171960334,
            "Faraday dual-port": 16.06033992848491,
        }
    )
    assert values["scan"]["usable_frame_range"] == {
        "Faraday dark-field": [0, 4],
        "Faraday dual-port": [0, 21],
    }


def test_figure_5_2_outputs_and_metadata_are_complete() -> None:
    required = (
        "figure_5_2.svg",
        "figure_5_2_data.csv",
        "figure_5_2_values.json",
        "figure_5_2_dual_port_heatmap.svg",
        "figure_5_2_operating_band.svg",
        "metadata.json",
    )
    assert all((OUTPUT / name).exists() and (OUTPUT / name).stat().st_size > 0 for name in required)

    metadata = json.loads((OUTPUT / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["model_config"].replace("\\", "/") == "configs/dissertation_v2_dcc3260m.json"
    assert metadata["criteria"]["minimum_snr"] == 9.0
    assert "uncalibrated" in metadata["scope"]["faraday"]
