from __future__ import annotations

import json
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

from scripts.recover_notebook_camera_stage import build_camera_stage, comparison_report


CONFIG_PATH = Path("configs/notebook_v1_defaults.json")


def _load_config() -> dict:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def test_notebook_camera_recovery_matches_helper_quantities() -> None:
    config = _load_config()
    stage = build_camera_stage(config)
    report = comparison_report(config, stage)

    assert stage["input_ideal_image"].shape == (1024, 1024)
    assert stage["notebook_binned_image"].shape == (68, 68)
    assert stage["notebook_camera_image"].shape == (68, 68)
    assert stage["trimmed_shape"] == [1020, 1020]

    for key in ["binned_image", "deterministic_counts", "normalised_camera_image"]:
        assert report[key]["max_absolute_difference"] == 0.0
        assert report[key]["max_relative_difference"] == 0.0

    assert report["normalisation_identity"]["max_absolute_difference"] < 2e-16
    assert np.isfinite(stage["notebook_camera_image"]).all()
    assert np.isfinite(stage["notebook_deterministic_counts"]).all()


def test_notebook_camera_recovery_stable_reference_values() -> None:
    config = _load_config()
    stage = build_camera_stage(config)
    camera = stage["notebook_camera_image"]
    counts = stage["notebook_deterministic_counts"]

    assert stage["input_stage"] == "pci_stage"
    assert stage["probe_power_mw"] == 2.0
    assert stage["bin_size"] == 15
    assert stage["default_exposure_s"] == pytest.approx(1e-4, rel=1e-12)
    assert stage["quantum_efficiency"] == 0.4
    assert stage["read_noise_electrons"] == 7.0
    assert stage["photon_scale"]["magnification"] == 2.0
    assert stage["photon_scale"]["object_pixel_m"] == pytest.approx(1.465e-6, rel=1e-12)
    assert stage["photon_scale"]["photons_per_pixel"] == pytest.approx(
        1532.3236613066642,
        rel=1e-12,
    )

    assert camera[34, 34] == pytest.approx(1.1339167724448815, rel=1e-12)
    assert np.max(camera) == pytest.approx(1.1339167724448815, rel=1e-12)
    assert np.min(camera) == pytest.approx(0.8662054798811654, rel=1e-12)
    assert np.mean(camera) == pytest.approx(0.9043924044803204, rel=1e-12)
    assert np.unravel_index(np.argmax(camera), camera.shape) == (34, 34)

    assert counts[34, 34] == pytest.approx(1737.5275003697766, rel=1e-12)
    assert np.max(counts) == pytest.approx(1737.5275003697766, rel=1e-12)
    assert np.min(counts) == pytest.approx(1327.3071523754033, rel=1e-12)
    assert np.mean(counts) == pytest.approx(1385.821880491222, rel=1e-12)


def test_notebook_camera_recovery_is_deterministic_noiseless_stage() -> None:
    config = _load_config()
    first = build_camera_stage(config)
    second = build_camera_stage(config)

    np.testing.assert_allclose(first["notebook_camera_image"], second["notebook_camera_image"])
    np.testing.assert_allclose(first["notebook_deterministic_counts"], second["notebook_deterministic_counts"])
    np.testing.assert_allclose(first["notebook_camera_image"], first["notebook_binned_image"])
