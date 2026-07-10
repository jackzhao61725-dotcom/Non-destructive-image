from __future__ import annotations

import json
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

from scripts.recover_notebook_faraday_stage import build_faraday_stage, comparison_report


CONFIG_PATH = Path("configs/notebook_v1_defaults.json")


def _load_config() -> dict:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def test_notebook_faraday_recovery_matches_helper_quantities() -> None:
    config = _load_config()
    stage = build_faraday_stage(config)
    report = comparison_report(config, stage)

    assert stage["notebook_theta_f_map_rad"].shape == (1024, 1024)
    for key in [
        "theta_f_map",
        "sigma_plus_object_field",
        "sigma_minus_object_field",
        "sigma_plus_field",
        "sigma_minus_field",
        "output_ex_field",
        "output_ey_field",
        "dark_field_intensity",
        "dual_port_u_intensity",
        "dual_port_v_intensity",
        "dual_port_signal",
    ]:
        assert report[key]["max_absolute_difference"] == 0.0
        assert report[key]["max_relative_difference"] == 0.0


def test_notebook_faraday_recovery_stable_reference_values() -> None:
    config = _load_config()
    stage = build_faraday_stage(config)
    dark = stage["notebook_dark_field_intensity"]
    signal = stage["notebook_dual_port_signal"]

    assert stage["kappa_F"] == 1.0
    assert stage["theta_f_peak_rad"] == pytest.approx(0.20294165287929014, rel=1e-12)
    assert np.count_nonzero(stage["notebook_pupil"]) == 1245
    assert dark[512, 512] == pytest.approx(0.015956443479481053, rel=1e-12)
    assert np.max(dark) == pytest.approx(0.015956443479481053, rel=1e-12)
    assert np.mean(dark) == pytest.approx(8.457346602539056e-05, rel=1e-12)
    assert signal[512, 512] == pytest.approx(0.25116900542198173, rel=1e-12)
    assert np.max(signal) == pytest.approx(0.25116900542198173, rel=1e-12)
    assert np.min(signal) == pytest.approx(-0.04479708725007998, rel=1e-12)
    assert np.mean(signal) == pytest.approx(0.001879360319011973, rel=1e-12)
    assert np.unravel_index(np.argmax(dark), dark.shape) == (512, 512)
    assert np.unravel_index(np.argmax(signal), signal.shape) == (512, 512)
    assert np.isfinite(dark).all()
    assert np.isfinite(signal).all()
    assert np.all(dark >= 0)


def test_notebook_faraday_recovery_records_compatible_existing_baseline() -> None:
    config = _load_config()
    stage = build_faraday_stage(config)
    report = comparison_report(config, stage)
    baseline = report["baseline_comparison"]

    assert baseline["available"] is True
    assert baseline["metadata"]["baseline_name"] == "faraday_imaging_baseline_v1"
    assert baseline["metadata"]["kappa_F"] == 1.0
    assert baseline["metadata"]["microscopic_faraday_model"] is False
    assert baseline["arrays"]["theta_f_map_rad"]["max_absolute_difference"] < 3e-5
    assert baseline["arrays"]["dark_field_intensity"]["max_absolute_difference"] < 5e-6
    assert baseline["arrays"]["dual_port_signal"]["max_absolute_difference"] < 4e-5
