from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.generate_figure_5_2 import DEFAULT_CONFIG, generate


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "results" / "dissertation_plots_v2_orca_fusion" / "figure_5_2"


def _values() -> dict:
    return json.loads((OUTPUT / "figure_5_2_values.json").read_text(encoding="utf-8"))


def test_figure_5_2_reference_contract() -> None:
    values = _values()
    reference = values["canonical_gate"]["reference"]

    assert values["canonical_gate"]["passed"] is True
    assert reference["depletion_limited_frames"] == 3
    assert reference["usable_frames"] == {
        "Faraday dark-field": 0,
        "Faraday dual-port": 3,
    }
    assert reference["initial_snr"] == pytest.approx(
        {
            "Faraday dark-field": 7.8000308699781264,
            "Faraday dual-port": 15.344634363699267,
        }
    )
    assert values["scan"]["usable_frame_range"] == {
        "Faraday dark-field": [0, 2],
        "Faraday dual-port": [0, 7],
    }
    assert values["criteria"]["central_block_shape"] == [5, 5]
    assert {
        label: condition["usable_frames"]["Faraday dual-port"]
        for label, condition in values["selected_conditions"].items()
    } == {"A": 3, "B": 6, "C": 3}


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
    assert metadata["model_config"].replace("\\", "/") == "configs/dissertation_v3_orca_fusion.json"
    assert metadata["criteria"]["minimum_snr"] == 8.0
    assert "kappa_F=-45/91" in metadata["scope"]["faraday"]
    assert "sigma_r=0.7 e-" in metadata["scope"]["detector"]
    assert "NA=0.130" in metadata["scope"]["optics"]
    assert metadata["result_status"] == {
        "status": "frozen_version_1_screening_output",
        "current_oxford_prediction": False,
        "canonical_regeneration_allowed": False,
        "regeneration_blocked_pending": "docs/multiframe_heating_model_optimisation.md",
    }

    correction = metadata["provenance_correction"]
    assert correction["generation_base_commit"] == "df00b42ae509921e663cccc0adc0a7a7a240c2a8"
    assert correction["generation_base_commit_status"] == "pre-rewrite local checkpoint identifier"
    assert correction["generation_worktree_dirty"] is True
    assert correction["original_dirty_diff_sealed"] is False
    assert correction["metadata_annotation_only"] is True
    assert len(correction["retrospective_source_snapshot_commit"]) == 40
    assert correction["retrospective_source_snapshot_status"] == (
        "post-generation preserved source/config snapshot; not a clean-generation claim"
    )


def test_figure_5_2_canonical_regeneration_is_blocked() -> None:
    with pytest.raises(RuntimeError, match="canonical Figure 5.2 is frozen"):
        generate(DEFAULT_CONFIG)
