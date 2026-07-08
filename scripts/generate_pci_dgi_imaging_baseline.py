from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path
from typing import Any

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK_PATH = "1 calculations revised 2  multishot  6  extended.ipynb"
OUTPUT_PATH = REPO_ROOT / "regression" / "baseline" / "imaging" / "pci_dgi_imaging_baseline_v1.npz"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_baseline_arrays() -> dict[str, Any]:
    """Build the notebook-equivalent deterministic PCI/DGI imaging arrays.

    This intentionally mirrors notebook section 7.3, especially cell 18:
    ``object_field = exp(1j * phase_map)``, scattered field
    ``object_field - 1``, propagation ``ifft2(fft2(scattered) * pupil)``,
    PCI reference ``t_p * exp(1j * theta)``, and DGI reference
    ``10**(-OD/2)``.
    """

    hbar = 1.054571817e-34
    h = 2 * np.pi * hbar
    c = 2.99792458e8
    k_b = 1.380649e-23
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
    t_p = 0.95
    theta_rad = np.pi / 2
    detuning_hz = 1.5e9
    dgi_od = 4.0
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
    phase_peak_rad = sigma0_m2 * column_density_m2[imaging_axis] * delta / (2 * (1 + delta**2))

    dgrid_m = fov_m / ngrid
    grid_axis_m = (np.arange(ngrid) - ngrid // 2) * dgrid_m
    grid_a_m, grid_b_m = np.meshgrid(grid_axis_m, grid_axis_m)
    plane = [index for index in range(3) if index != imaging_axis]
    thomas_fermi_profile = np.maximum(
        0,
        1 - grid_a_m**2 / radii_m[plane[0]] ** 2 - grid_b_m**2 / radii_m[plane[1]] ** 2,
    ) ** 1.5

    scalar_phase_map_rad = phase_peak_rad * thomas_fermi_profile
    object_field = np.exp(1j * scalar_phase_map_rad)
    scattered_field = object_field - 1

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

    propagated_scattered_field = np.fft.ifft2(np.fft.fft2(scattered_field) * pupil)
    pci_reference_field = np.array(t_p * np.exp(1j * theta_rad))
    dgi_reference_field = np.array(10 ** (-dgi_od / 2))
    pci_image_intensity = np.abs(pci_reference_field + propagated_scattered_field) ** 2
    dgi_image_intensity = np.abs(dgi_reference_field + propagated_scattered_field) ** 2

    metadata = {
        "baseline_name": "pci_dgi_imaging_baseline_v1",
        "baseline_type": "representative deterministic notebook-equivalent PCI/DGI imaging arrays",
        "source_notebook": NOTEBOOK_PATH,
        "source_notebook_sha256": _sha256(REPO_ROOT / NOTEBOOK_PATH),
        "notebook_reference": "Section 7.3 / cell 18 sim_image scalar PCI-DGI Fourier-optics implementation",
        "full_notebook_execution": False,
        "notes": [
            "Generated from notebook-equivalent formulas and parameters, not by executing the full notebook.",
            "No stochastic camera noise is included.",
            "This baseline is intended to lock the current scalar PCI/DGI Fourier-imaging convention before orchestration migration.",
        ],
        "imaging_axis": "x",
        "transverse_plane": ["y", "z"],
        "detuning_hz": detuning_hz,
        "dgi_od": dgi_od,
        "ngrid": ngrid,
        "fov_m": fov_m,
        "phase_peak_rad": float(phase_peak_rad),
        "pci_reference_field_real": float(np.real(pci_reference_field)),
        "pci_reference_field_imag": float(np.imag(pci_reference_field)),
        "dgi_reference_field": float(dgi_reference_field),
        "t_p": t_p,
        "theta_rad": theta_rad,
        "numerical_aperture": numerical_aperture,
        "fft_convention": "np.fft.ifft2(np.fft.fft2(scattered_field) * pupil)",
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
        "phase_peak_rad": np.array(phase_peak_rad),
        "thomas_fermi_profile": thomas_fermi_profile,
        "scalar_phase_map_rad": scalar_phase_map_rad,
        "object_field": object_field,
        "pupil": pupil,
        "scattered_field": scattered_field,
        "propagated_scattered_field": propagated_scattered_field,
        "pci_reference_field": pci_reference_field,
        "pci_image_intensity": pci_image_intensity,
        "dgi_reference_field": dgi_reference_field,
        "dgi_image_intensity": dgi_image_intensity,
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

