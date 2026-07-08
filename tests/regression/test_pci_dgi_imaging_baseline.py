from __future__ import annotations

import json
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

from scripts.generate_pci_dgi_imaging_baseline import build_baseline_arrays


BASELINE_PATH = Path("regression/baseline/imaging/pci_dgi_imaging_baseline_v1.npz")

REQUIRED_KEYS = {
    "grid_axis_m",
    "grid_a_m",
    "grid_b_m",
    "spatial_frequency_axis_m_inv",
    "spatial_frequency_x_m_inv",
    "spatial_frequency_y_m_inv",
    "radii_m",
    "column_density_m2",
    "phase_peak_rad",
    "thomas_fermi_profile",
    "scalar_phase_map_rad",
    "object_field",
    "pupil",
    "scattered_field",
    "propagated_scattered_field",
    "pci_reference_field",
    "pci_image_intensity",
    "dgi_reference_field",
    "dgi_image_intensity",
    "metadata_json",
}

ARRAY_2D_KEYS = {
    "grid_a_m",
    "grid_b_m",
    "spatial_frequency_x_m_inv",
    "spatial_frequency_y_m_inv",
    "thomas_fermi_profile",
    "scalar_phase_map_rad",
    "object_field",
    "pupil",
    "scattered_field",
    "propagated_scattered_field",
    "pci_image_intensity",
    "dgi_image_intensity",
}


def _load_baseline() -> np.lib.npyio.NpzFile:
    assert BASELINE_PATH.exists()
    return np.load(BASELINE_PATH)


def test_pci_dgi_imaging_baseline_file_and_keys_exist() -> None:
    with _load_baseline() as baseline:
        assert REQUIRED_KEYS.issubset(set(baseline.files))


def test_pci_dgi_imaging_baseline_shapes_and_values() -> None:
    with _load_baseline() as baseline:
        for key in ARRAY_2D_KEYS:
            assert baseline[key].shape == (1024, 1024)
            assert np.isfinite(baseline[key]).all()

        assert baseline["grid_axis_m"].shape == (1024,)
        assert baseline["spatial_frequency_axis_m_inv"].shape == (1024,)
        assert baseline["radii_m"].shape == (3,)
        assert baseline["column_density_m2"].shape == (3,)
        assert np.isfinite(baseline["phase_peak_rad"])
        assert np.isfinite(baseline["pci_reference_field"])
        assert np.isfinite(baseline["dgi_reference_field"])
        assert np.isfinite(baseline["radii_m"]).all()
        assert np.isfinite(baseline["column_density_m2"]).all()
        assert np.all(baseline["pci_image_intensity"] >= 0)
        assert np.all(baseline["dgi_image_intensity"] >= 0)


def test_pci_dgi_imaging_baseline_metadata_is_present() -> None:
    with _load_baseline() as baseline:
        metadata = json.loads(str(baseline["metadata_json"]))

    assert metadata["baseline_name"] == "pci_dgi_imaging_baseline_v1"
    assert metadata["full_notebook_execution"] is False
    assert metadata["fft_convention"] == "np.fft.ifft2(np.fft.fft2(scattered_field) * pupil)"
    assert metadata["imaging_axis"] == "x"
    assert metadata["transverse_plane"] == ["y", "z"]


def test_pci_dgi_imaging_baseline_regenerates_equal_arrays() -> None:
    regenerated = build_baseline_arrays()
    with _load_baseline() as baseline:
        for key in REQUIRED_KEYS:
            if key == "metadata_json":
                stored_metadata = json.loads(str(baseline[key]))
                regenerated_metadata = json.loads(str(regenerated[key]))
                for stable_key in [
                    "baseline_name",
                    "baseline_type",
                    "source_notebook_sha256",
                    "notebook_reference",
                    "full_notebook_execution",
                    "imaging_axis",
                    "transverse_plane",
                    "detuning_hz",
                    "dgi_od",
                    "ngrid",
                    "fov_m",
                    "phase_peak_rad",
                    "pci_reference_field_real",
                    "pci_reference_field_imag",
                    "dgi_reference_field",
                    "t_p",
                    "theta_rad",
                    "numerical_aperture",
                    "fft_convention",
                ]:
                    assert regenerated_metadata[stable_key] == stored_metadata[stable_key]
                continue

            np.testing.assert_allclose(regenerated[key], baseline[key], rtol=1e-10, atol=1e-12)

