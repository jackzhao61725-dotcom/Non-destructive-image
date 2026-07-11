from __future__ import annotations

import json
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

from scripts.recover_notebook_faraday_camera_panel import build_faraday_camera_panel, comparison_report


CONFIG_PATH = Path("configs/notebook_v1_defaults.json")
OUTPUT_DIR = Path("results/notebook_aligned_recovery/faraday_camera_panel")


def _load_config() -> dict:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def test_faraday_camera_panel_outputs_exist_and_record_metadata() -> None:
    required = [
        "faraday_camera_panel.svg",
        "comparison_report.json",
        "faraday_camera_panel_summary.json",
        "metadata.json",
        "lineouts.csv",
        "frame_statistics.csv",
    ]
    for filename in required:
        assert (OUTPUT_DIR / filename).exists()

    metadata = json.loads((OUTPUT_DIR / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["notebook_counterpart"]["primary_cell"] == 51
    assert metadata["figure_type"] == "notebook-aligned recovery"
    assert metadata["physical_parameters"]["kappa_F"] == 1.0
    assert "phenomenological Version 1 placeholder" in metadata["physical_parameters"]["kappa_F_status"]
    assert metadata["comparison_status"]["exact_explicit_seed_helper_replay"] is True
    assert metadata["comparison_status"]["exact_arbitrary_notebook_global_rng_reproduction"] is False


def test_faraday_camera_panel_stage_has_expected_shapes_and_parameters() -> None:
    config = _load_config()
    stage = build_faraday_camera_panel(config)

    assert stage["dark"]["noisy_image"].shape == (68, 68)
    assert stage["dark"]["camera_image"].shape == (68, 68)
    assert stage["port_u"]["noisy_image"].shape == (68, 68)
    assert stage["port_v"]["noisy_image"].shape == (68, 68)
    assert stage["dual_port_noisy_signal"].shape == (68, 68)
    assert stage["dual_port_ideal_signal"].shape == (68, 68)

    assert stage["photon_scale"]["probe_power_mw"] == 5.0
    assert stage["photon_scale"]["photons_per_pixel"] == pytest.approx(3830.8091532666604, rel=1e-12)
    assert stage["read_noise_electrons"] == 7.0
    assert stage["rng_seed"] == 7
    assert np.isfinite(stage["dual_port_noisy_signal"]).all()


def test_faraday_camera_panel_comparison_report_is_stable() -> None:
    config = _load_config()
    stage = build_faraday_camera_panel(config)
    report = comparison_report(config, stage)

    assert report["notebook_counterpart"]["primary_cell"] == 51
    assert report["camera_binning"]["camera_shape"] == [68, 68]
    assert report["noise"]["exact_explicit_seed_replay"] is True
    assert report["dark_field_noisy"]["max_absolute_difference"] == 0.0
    assert report["dual_port_noisy_signal"]["max_absolute_difference"] == 0.0
    assert report["dark_field_noisy"]["stats"]["centre_value"] == pytest.approx(
        0.01617456315961192,
        rel=1e-12,
    )
    assert report["dual_port_noisy_signal"]["stats"]["centre_value"] == pytest.approx(
        0.23867929119031744,
        rel=1e-12,
    )
