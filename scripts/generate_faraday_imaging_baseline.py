from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path
from typing import Any

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK_PATH = "1 calculations revised 2  multishot  6  extended.ipynb"
OUTPUT_PATH = REPO_ROOT / "regression" / "baseline" / "imaging" / "faraday_imaging_baseline_v1.npz"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_baseline_arrays() -> dict[str, Any]:
    """Build notebook-equivalent deterministic Faraday imaging arrays.

    This mirrors notebook section 17.2 / cell 51. It preserves the current
    phenomenological Faraday model ``theta_F = kappa_F * phi_peak`` with
    ``kappa_F = 1.0`` and the circular-component propagation / linear
    recombination conventions used by ``sim_faraday_fields``.
    """

    hbar = 1.054571817e-34
    amu = 1.66053907e-27
    a0 = 5.29177211e-11

    atomic_mass = 166 * amu
    wavelength_m = 401.0e-9
    gamma_rad_per_s = 2 * np.pi * 29.5e6
    atom_number = 2.5e4
    scattering_length_m = 72 * a0
    trap_frequencies_hz = np.array([293.0, 14.0, 233.0])

    probe_diameter_m = 24.0e-3
    f1_m = 150e-3
    detuning_hz = 1.5e9
    kappa_f = 1.0
    imaging_axis = 0
    ngrid = 1024
    fov_m = 100e-6

    sigma0_m2 = 3 * wavelength_m**2 / (2 * np.pi)
    trap_angular_frequencies = 2 * np.pi * trap_frequencies_hz
    omega_bar = trap_angular_frequencies.prod() ** (1 / 3)
    harmonic_oscillator_length_m = np.sqrt(hbar / (atomic_mass * omega_bar))
    chemical_potential_j = (
        0.5
        * (15 * atom_number * scattering_length_m / harmonic_oscillator_length_m) ** (2 / 5)
        * hbar
        * omega_bar
    )
    peak_density_m3 = chemical_potential_j * atomic_mass / (4 * np.pi * hbar**2 * scattering_length_m)
    radii_m = np.sqrt(2 * chemical_potential_j / (atomic_mass * trap_angular_frequencies**2))
    column_density_m2 = (4 / 3) * peak_density_m3 * radii_m

    delta = 2 * detuning_hz * 2 * np.pi / gamma_rad_per_s
    scalar_phase_peak_rad = sigma0_m2 * column_density_m2[imaging_axis] * delta / (2 * (1 + delta**2))
    theta_f_peak_rad = kappa_f * scalar_phase_peak_rad

    dgrid_m = fov_m / ngrid
    grid_axis_m = (np.arange(ngrid) - ngrid // 2) * dgrid_m
    grid_a_m, grid_b_m = np.meshgrid(grid_axis_m, grid_axis_m)
    plane = [index for index in range(3) if index != imaging_axis]
    thomas_fermi_profile = np.maximum(
        0,
        1 - grid_a_m**2 / radii_m[plane[0]] ** 2 - grid_b_m**2 / radii_m[plane[1]] ** 2,
    ) ** 1.5

    scalar_phase_map_rad = scalar_phase_peak_rad * thomas_fermi_profile
    theta_f_map_rad = theta_f_peak_rad * thomas_fermi_profile

    numerical_aperture = (probe_diameter_m / 2) / f1_m
    spatial_frequency_axis_m_inv = np.fft.fftfreq(ngrid, dgrid_m)
    spatial_frequency_x_m_inv, spatial_frequency_y_m_inv = np.meshgrid(
        spatial_frequency_axis_m_inv,
        spatial_frequency_axis_m_inv,
    )
    pupil = (
        np.sqrt(spatial_frequency_x_m_inv**2 + spatial_frequency_y_m_inv**2)
        <= numerical_aperture / wavelength_m
    ).astype(float)

    sigma_plus_object_field = np.exp(1j * theta_f_map_rad)
    sigma_minus_object_field = np.exp(-1j * theta_f_map_rad)
    sigma_plus_scattered_field = sigma_plus_object_field - 1
    sigma_minus_scattered_field = sigma_minus_object_field - 1
    sigma_plus_propagated_scattered_field = np.fft.ifft2(np.fft.fft2(sigma_plus_scattered_field) * pupil)
    sigma_minus_propagated_scattered_field = np.fft.ifft2(np.fft.fft2(sigma_minus_scattered_field) * pupil)
    sigma_plus_field = 1 + sigma_plus_propagated_scattered_field
    sigma_minus_field = 1 + sigma_minus_propagated_scattered_field

    output_ex_field = (sigma_plus_field + sigma_minus_field) / 2
    output_ey_field = 1j * (sigma_plus_field - sigma_minus_field) / 2
    dark_field_intensity = np.abs(output_ey_field) ** 2
    dual_port_u_intensity = np.abs(output_ex_field + output_ey_field) ** 2 / 2
    dual_port_v_intensity = np.abs(output_ex_field - output_ey_field) ** 2 / 2
    dual_port_signal = (dual_port_v_intensity - dual_port_u_intensity) / (
        dual_port_v_intensity + dual_port_u_intensity
    )

    metadata = {
        "baseline_name": "faraday_imaging_baseline_v1",
        "baseline_type": "representative deterministic notebook-equivalent Faraday imaging arrays",
        "source_notebook": NOTEBOOK_PATH,
        "source_notebook_sha256": _sha256(REPO_ROOT / NOTEBOOK_PATH),
        "notebook_reference": "Section 17.2 / cell 51 sim_faraday_fields and faraday_maps implementation",
        "full_notebook_execution": False,
        "microscopic_faraday_model": False,
        "notes": [
            "Generated from notebook-equivalent formulas and parameters, not by executing the full notebook.",
            "No stochastic camera noise is included.",
            "The current phenomenological kappa_F placeholder is preserved exactly.",
        ],
        "imaging_axis": "x",
        "transverse_plane": ["y", "z"],
        "detuning_hz": detuning_hz,
        "kappa_F": kappa_f,
        "ngrid": ngrid,
        "fov_m": fov_m,
        "scalar_phase_peak_rad": float(scalar_phase_peak_rad),
        "theta_f_peak_rad": float(theta_f_peak_rad),
        "numerical_aperture": numerical_aperture,
        "circular_component_convention": "sigma_plus=exp(+1j*theta_F_map), sigma_minus=exp(-1j*theta_F_map)",
        "propagation_convention": "1 + np.fft.ifft2(np.fft.fft2(exp(+-1j*theta_F_map)-1) * pupil)",
        "linear_recombination_convention": "Ex=(Pp+Pm)/2, Ey=1j*(Pp-Pm)/2",
        "dark_field_convention": "I_dark=abs(Ey)**2",
        "dual_port_convention": "I_u=abs(Ex+Ey)**2/2, I_v=abs(Ex-Ey)**2/2, S=(I_v-I_u)/(I_v+I_u)",
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
    }

    return {
        "grid_axis_m": grid_axis_m,
        "grid_a_m": grid_a_m,
        "grid_b_m": grid_b_m,
        "spatial_frequency_axis_m_inv": spatial_frequency_axis_m_inv,
        "spatial_frequency_x_m_inv": spatial_frequency_x_m_inv,
        "spatial_frequency_y_m_inv": spatial_frequency_y_m_inv,
        "radii_m": radii_m,
        "column_density_m2": column_density_m2,
        "scalar_phase_peak_rad": np.array(scalar_phase_peak_rad),
        "theta_f_peak_rad": np.array(theta_f_peak_rad),
        "thomas_fermi_profile": thomas_fermi_profile,
        "scalar_phase_map_rad": scalar_phase_map_rad,
        "theta_f_map_rad": theta_f_map_rad,
        "pupil": pupil,
        "sigma_plus_object_field": sigma_plus_object_field,
        "sigma_minus_object_field": sigma_minus_object_field,
        "sigma_plus_scattered_field": sigma_plus_scattered_field,
        "sigma_minus_scattered_field": sigma_minus_scattered_field,
        "sigma_plus_propagated_scattered_field": sigma_plus_propagated_scattered_field,
        "sigma_minus_propagated_scattered_field": sigma_minus_propagated_scattered_field,
        "sigma_plus_field": sigma_plus_field,
        "sigma_minus_field": sigma_minus_field,
        "output_ex_field": output_ex_field,
        "output_ey_field": output_ey_field,
        "dark_field_intensity": dark_field_intensity,
        "dual_port_u_intensity": dual_port_u_intensity,
        "dual_port_v_intensity": dual_port_v_intensity,
        "dual_port_signal": dual_port_signal,
        "metadata_json": np.array(json.dumps(metadata, sort_keys=True)),
    }


def write_baseline(output_path: Path = OUTPUT_PATH) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_path, **build_baseline_arrays())
    return output_path


def main() -> None:
    output_path = write_baseline()
    print(output_path.relative_to(REPO_ROOT).as_posix())


if __name__ == "__main__":
    main()

