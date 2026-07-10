from __future__ import annotations

import json
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

from scripts.recover_notebook_dgi_stage import build_dgi_stage, comparison_report


CONFIG_PATH = Path("configs/notebook_v1_defaults.json")


def _load_config() -> dict:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def test_notebook_dgi_recovery_matches_helper_quantities() -> None:
    config = _load_config()
    stage = build_dgi_stage(config)
    report = comparison_report(config, stage)

    assert stage["notebook_dgi_image_intensity"].shape == (1024, 1024)
    assert report["object_field"]["max_absolute_difference"] == 0.0
    assert report["object_field"]["max_relative_difference"] == 0.0
    assert report["dgi_stop"]["reference_field_absolute_difference"] == 0.0
    assert report["dgi_image_intensity"]["max_absolute_difference"] == 0.0
    assert report["dgi_image_intensity"]["max_relative_difference"] == 0.0
    assert report["propagated_scattered_field"]["max_absolute_difference"] < 1e-18
    assert report["propagated_scattered_field"]["max_relative_difference"] < 1e-10


def test_notebook_dgi_recovery_stable_reference_values() -> None:
    config = _load_config()
    stage = build_dgi_stage(config)
    intensity = stage["notebook_dgi_image_intensity"]

    assert stage["stop_optical_depth"] == 4.0
    assert stage["notebook_reference_field"] == pytest.approx(0.01, rel=1e-15)
    assert stage["dgi_reference_intensity"] == pytest.approx(0.0001, rel=1e-15)
    assert stage["pupil_stage"]["numerical_aperture"] == pytest.approx(0.08, rel=1e-15)
    assert np.count_nonzero(stage["notebook_pupil"]) == 1245
    assert intensity[512, 512] == pytest.approx(0.01595651906354294, rel=1e-12)
    assert np.max(intensity) == pytest.approx(0.01595651906354294, rel=1e-12)
    assert np.min(intensity) == pytest.approx(9.854763408159256e-05, rel=1e-12)
    assert np.mean(intensity) == pytest.approx(0.00018381101889794308, rel=1e-12)
    assert np.unravel_index(np.argmax(intensity), intensity.shape) == (512, 512)
    assert np.isfinite(intensity).all()
    assert np.all(intensity >= 0)
