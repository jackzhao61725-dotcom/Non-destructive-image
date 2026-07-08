from __future__ import annotations

import json
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

from non_destructive_image import simulate_pci_image


BASELINE_PATH = Path("regression/baseline/imaging/pci_dgi_imaging_baseline_v1.npz")


def test_simulate_pci_image_matches_pci_dgi_baseline() -> None:
    with np.load(BASELINE_PATH) as baseline:
        metadata = json.loads(str(baseline["metadata_json"]))
        result = simulate_pci_image(
            baseline["scalar_phase_map_rad"],
            baseline["pupil"],
            phase_plate_transmittance=metadata["t_p"],
            phase_plate_phase=metadata["theta_rad"],
            return_intermediates=True,
        )

        np.testing.assert_allclose(result["object_field"], baseline["object_field"], rtol=1e-10, atol=1e-12)
        np.testing.assert_allclose(result["scattered_field"], baseline["scattered_field"], rtol=1e-10, atol=1e-12)
        np.testing.assert_allclose(
            result["propagated_scattered_field"],
            baseline["propagated_scattered_field"],
            rtol=1e-10,
            atol=1e-12,
        )
        np.testing.assert_allclose(result["pci_reference_field"], baseline["pci_reference_field"], rtol=1e-10, atol=1e-12)
        np.testing.assert_allclose(result["pci_image_intensity"], baseline["pci_image_intensity"], rtol=1e-10, atol=1e-12)


def test_simulate_pci_image_returns_valid_intensity() -> None:
    with np.load(BASELINE_PATH) as baseline:
        metadata = json.loads(str(baseline["metadata_json"]))
        intensity = simulate_pci_image(
            baseline["scalar_phase_map_rad"],
            baseline["pupil"],
            phase_plate_transmittance=metadata["t_p"],
            phase_plate_phase=metadata["theta_rad"],
        )

    assert intensity.shape == (1024, 1024)
    assert np.isfinite(intensity).all()
    assert np.all(intensity >= 0)

