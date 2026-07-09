from __future__ import annotations

import json
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

from scripts.recover_notebook_phase_stage import build_phase_stage, comparison_report


CONFIG_PATH = Path("configs/notebook_v1_defaults.json")


def _load_config() -> dict:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def test_notebook_phase_recovery_matches_helper_quantities() -> None:
    config = _load_config()
    stage = build_phase_stage(config)
    report = comparison_report(config, stage)

    assert stage["notebook_phase_map_rad"].shape == (1024, 1024)
    assert report["phase_peak_rad"]["absolute_difference"] == 0.0
    assert report["phase_map_comparison"]["max_absolute_difference"] == 0.0
    assert report["phase_map_comparison"]["max_relative_difference"] == 0.0
    assert report["dimensionless_detuning"]["absolute_difference"] == 0.0


def test_notebook_phase_recovery_stable_reference_values() -> None:
    config = _load_config()
    stage = build_phase_stage(config)
    phase_map = stage["notebook_phase_map_rad"]

    assert stage["detuning_hz"] == 1.5e9
    assert stage["notebook_phase_peak_rad"] == pytest.approx(0.20294165287929014, rel=1e-12)
    assert phase_map[512, 512] == pytest.approx(stage["notebook_phase_peak_rad"], rel=1e-12)
    assert np.unravel_index(np.argmax(phase_map), phase_map.shape) == (512, 512)
    assert np.min(phase_map) == 0.0
    assert np.isfinite(phase_map).all()
