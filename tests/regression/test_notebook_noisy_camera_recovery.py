from __future__ import annotations

import json
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

from scripts.recover_notebook_noisy_camera_stage import build_noisy_camera_stage, comparison_report


CONFIG_PATH = Path("configs/notebook_v1_defaults.json")


def _load_config() -> dict:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def test_notebook_noisy_camera_recovery_matches_helper_for_seeded_rng() -> None:
    config = _load_config()
    stage = build_noisy_camera_stage(config)
    report = comparison_report(config, stage)

    assert stage["notebook_noisy_image"].shape == (68, 68)
    assert stage["notebook_noisy_counts"].shape == (68, 68)
    assert report["rng_policy"]["rng_seed"] == 7
    assert report["rng_policy"]["hidden_rng_introduced"] is False

    for key in ["binned_input", "noisy_counts", "noisy_camera_image"]:
        assert report[key]["max_absolute_difference"] == 0.0
        assert report[key]["max_relative_difference"] == 0.0


def test_notebook_noisy_camera_recovery_seed_policy_and_noise_statistics() -> None:
    config = _load_config()
    stage = build_noisy_camera_stage(config)
    noisy = stage["notebook_noisy_image"]
    deterministic = stage["notebook_deterministic_camera_image"]
    residual = stage["notebook_residual"]

    assert stage["rng_seed"] == 7
    assert stage["read_noise_electrons"] == 7.0
    assert stage["photons_per_pixel"] == pytest.approx(1532.3236613066642, rel=1e-12)
    assert np.isfinite(noisy).all()
    assert np.isfinite(stage["notebook_noisy_counts"]).all()

    np.testing.assert_allclose(noisy, build_noisy_camera_stage(config)["notebook_noisy_image"])
    assert np.max(np.abs(noisy - stage["different_seed_noisy_image"])) > 0.05
    assert abs(float(np.mean(noisy)) - float(np.mean(deterministic))) < 0.002
    assert np.std(residual) == pytest.approx(0.02517625607875563, rel=1e-12)
    assert np.std(residual) / np.mean(stage["notebook_expected_image_std"]) == pytest.approx(
        1.0184981488308726,
        rel=1e-12,
    )


def test_notebook_noisy_camera_recovery_stable_reference_values() -> None:
    config = _load_config()
    stage = build_noisy_camera_stage(config)
    noisy = stage["notebook_noisy_image"]
    counts = stage["notebook_noisy_counts"]

    assert noisy[34, 34] == pytest.approx(1.1547265453246764, rel=1e-12)
    assert np.max(noisy) == pytest.approx(1.1616523905537648, rel=1e-12)
    assert np.min(noisy) == pytest.approx(0.8196905139413712, rel=1e-12)
    assert np.mean(noisy) == pytest.approx(0.9047881903189995, rel=1e-12)
    assert np.unravel_index(np.argmax(noisy), noisy.shape) == (34, 38)

    assert counts[34, 34] == pytest.approx(1769.414807739904, rel=1e-12)
    assert np.max(counts) == pytest.approx(1780.027444258984, rel=1e-12)
    assert np.min(counts) == pytest.approx(1256.0311694609832, rel=1e-12)
    assert np.mean(counts) == pytest.approx(1386.4283524966404, rel=1e-12)
