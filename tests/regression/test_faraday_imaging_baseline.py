from __future__ import annotations

import json
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

from scripts.generate_faraday_imaging_baseline import build_baseline_arrays


BASELINE_PATH = Path("regression/baseline/imaging/faraday_imaging_baseline_v1.npz")

REQUIRED_KEYS = {
    "grid_axis_m",
    "grid_a_m",
    "grid_b_m",
    "spatial_frequency_axis_m_inv",
    "spatial_frequency_x_m_inv",
    "spatial_frequency_y_m_inv",
    "radii_m",
    "column_density_m2",
    "scalar_phase_peak_rad",
    "theta_f_peak_rad",
    "thomas_fermi_profile",
    "scalar_phase_map_rad",
    "theta_f_map_rad",
    "pupil",
    "sigma_plus_object_field",
    "sigma_minus_object_field",
    "sigma_plus_scattered_field",
    "sigma_minus_scattered_field",
    "sigma_plus_propagated_scattered_field",
    "sigma_minus_propagated_scattered_field",
    "sigma_plus_field",
    "sigma_minus_field",
    "output_ex_field",
    "output_ey_field",
    "dark_field_intensity",
    "dual_port_u_intensity",
    "dual_port_v_intensity",
    "dual_port_signal",
    "metadata_json",
}

ARRAY_2D_KEYS = {
    "grid_a_m",
    "grid_b_m",
    "spatial_frequency_x_m_inv",
    "spatial_frequency_y_m_inv",
    "thomas_fermi_profile",
    "scalar_phase_map_rad",
    "theta_f_map_rad",
    "pupil",
    "sigma_plus_object_field",
    "sigma_minus_object_field",
    "sigma_plus_scattered_field",
    "sigma_minus_scattered_field",
    "sigma_plus_propagated_scattered_field",
    "sigma_minus_propagated_scattered_field",
    "sigma_plus_field",
    "sigma_minus_field",
    "output_ex_field",
    "output_ey_field",
    "dark_field_intensity",
    "dual_port_u_intensity",
    "dual_port_v_intensity",
    "dual_port_signal",
}

INTENSITY_KEYS = {
    "dark_field_intensity",
    "dual_port_u_intensity",
    "dual_port_v_intensity",
}


def _load_baseline() -> np.lib.npyio.NpzFile:
    assert BASELINE_PATH.exists()
    return np.load(BASELINE_PATH)


def test_faraday_imaging_baseline_file_and_keys_exist() -> None:
    with _load_baseline() as baseline:
        assert REQUIRED_KEYS.issubset(set(baseline.files))


def test_faraday_imaging_baseline_shapes_and_values() -> None:
    with _load_baseline() as baseline:
        for key in ARRAY_2D_KEYS:
            assert baseline[key].shape == (1024, 1024)
            assert np.isfinite(baseline[key]).all()

        assert baseline["grid_axis_m"].shape == (1024,)
        assert baseline["spatial_frequency_axis_m_inv"].shape == (1024,)
        assert baseline["radii_m"].shape == (3,)
        assert baseline["column_density_m2"].shape == (3,)
        assert np.isfinite(baseline["theta_f_peak_rad"])
        assert np.isfinite(baseline["theta_f_map_rad"]).all()
        assert np.isfinite(baseline["radii_m"]).all()
        assert np.isfinite(baseline["column_density_m2"]).all()

        for key in INTENSITY_KEYS:
            assert np.all(baseline[key] >= 0)


def test_faraday_imaging_baseline_metadata_is_present() -> None:
    with _load_baseline() as baseline:
        metadata = json.loads(str(baseline["metadata_json"]))

    assert metadata["baseline_name"] == "faraday_imaging_baseline_v1"
    assert metadata["full_notebook_execution"] is False
    assert metadata["microscopic_faraday_model"] is False
    assert metadata["kappa_F"] == 1.0
    assert metadata["linear_recombination_convention"] == "Ex=(Pp+Pm)/2, Ey=1j*(Pp-Pm)/2"
    assert metadata["imaging_axis"] == "x"
    assert metadata["transverse_plane"] == ["y", "z"]


def test_faraday_imaging_baseline_regenerates_equal_arrays() -> None:
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
                    "microscopic_faraday_model",
                    "imaging_axis",
                    "transverse_plane",
                    "detuning_hz",
                    "kappa_F",
                    "ngrid",
                    "fov_m",
                    "scalar_phase_peak_rad",
                    "theta_f_peak_rad",
                    "numerical_aperture",
                    "circular_component_convention",
                    "propagation_convention",
                    "linear_recombination_convention",
                    "dark_field_convention",
                    "dual_port_convention",
                ]:
                    assert regenerated_metadata[stable_key] == stored_metadata[stable_key]
                continue

            np.testing.assert_allclose(regenerated[key], baseline[key], rtol=1e-10, atol=1e-12)

