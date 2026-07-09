from __future__ import annotations

import json
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

from scripts.recover_notebook_pci_stage import build_pci_stage, comparison_report


CONFIG_PATH = Path("configs/notebook_v1_defaults.json")


def _load_config() -> dict:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def test_notebook_pci_recovery_matches_helper_quantities() -> None:
    config = _load_config()
    stage = build_pci_stage(config)
    report = comparison_report(config, stage)

    assert stage["notebook_pci_image_intensity"].shape == (1024, 1024)
    assert report["object_field"]["max_absolute_difference"] == 0.0
    assert report["object_field"]["max_relative_difference"] == 0.0
    assert report["phase_plate"]["reference_field_absolute_difference"] == 0.0
    assert report["pci_image_intensity"]["max_absolute_difference"] == 0.0
    assert report["pci_image_intensity"]["max_relative_difference"] == 0.0
    assert report["propagated_scattered_field"]["max_absolute_difference"] < 2e-16
    assert report["propagated_scattered_field"]["max_relative_difference"] < 5e-9


def test_notebook_pci_recovery_stable_reference_values() -> None:
    config = _load_config()
    stage = build_pci_stage(config)
    intensity = stage["notebook_pci_image_intensity"]

    assert stage["phase_plate_transmittance"] == 0.95
    assert stage["phase_plate_phase_rad"] == pytest.approx(np.pi / 2, rel=1e-15)
    assert stage["plate_background_intensity"] == pytest.approx(0.9025, rel=1e-15)
    assert stage["pupil_stage"]["numerical_aperture"] == pytest.approx(0.08, rel=1e-15)
    assert np.count_nonzero(stage["notebook_pupil"]) == 1245
    assert intensity[512, 512] == pytest.approx(1.1585677695076864, rel=1e-12)
    assert np.max(intensity) == pytest.approx(1.1585677695076864, rel=1e-12)
    assert np.min(intensity) == pytest.approx(0.8603460757846358, rel=1e-12)
    assert np.mean(intensity) == pytest.approx(0.904372478370771, rel=1e-12)
    assert np.unravel_index(np.argmax(intensity), intensity.shape) == (512, 512)
    assert np.isfinite(intensity).all()
    assert np.all(intensity >= 0)
