"""Audit and regenerate thesis-facing numerical conclusions from explicit parameter sets."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
from scipy.special import zeta

from scripts.recover_notebook_multishot_stage import (
    _basic_constants,
    build_blur_axis,
    dimensionless_detuning,
    intensity_at_atoms_notebook,
    pci_snr_pixel_for_phase,
    photons_per_camera_pixel,
    reabsorption_fraction,
    scalar_phase_peak,
    scattered_photons_per_atom,
    self_consistent_total_atoms,
)
from scripts.recover_notebook_pci_stage import build_pupil


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = REPO_ROOT / "configs" / "thesis_numerical_contract_v1.json"
OUTPUT_DIR = REPO_ROOT / "results" / "thesis_numerical_consistency_v1"
REPORT_PATH = REPO_ROOT / "docs" / "thesis_numerical_consistency_correction_report.md"


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def git_value(*args: str) -> str:
    commands = (["git", *args], [r"C:\Program Files\Git\cmd\git.exe", *args])
    for command in commands:
        try:
            return subprocess.check_output(command, cwd=REPO_ROOT, text=True, stderr=subprocess.DEVNULL).strip()
        except (OSError, subprocess.CalledProcessError):
            continue
    return "unknown"


def photons_per_object_pixel(
    notebook: dict[str, Any], power_mw: float, exposure_s: float, quantum_efficiency: float
) -> float:
    constants = notebook["constants"]
    atom = notebook["atom"]
    geometry = notebook["imaging_geometry"]
    h_planck = 2 * np.pi * float(constants["hbar"])
    photon_energy = h_planck * float(constants["speed_of_light"]) / float(atom["transition_wavelength_m"])
    magnification = float(geometry["focal_length_2_m"]) / float(geometry["focal_length_1_m"])
    object_pixel_m = float(geometry["camera_pixel_m"]) / magnification
    return float(
        intensity_at_atoms_notebook(notebook, power_mw)
        * object_pixel_m**2
        * exposure_s
        * quantum_efficiency
        / photon_energy
    )


def clean_loss_budget(scattered: float, loss_fraction: float, eta_coll: float) -> float:
    return float(-np.log1p(-loss_fraction) / (eta_coll * scattered))


def heating_budget(
    notebook: dict[str, Any],
    constants: dict[str, Any],
    thermal: dict[str, float],
    scattered: float,
    loss_fraction: float,
    reabsorption: float,
) -> float:
    temperature = float(notebook["condensate"]["temperature_k"])
    critical_temperature = thermal["critical_temperature_k"]
    initial_fraction = 1 - (temperature / critical_temperature) ** 3
    target_temperature = critical_temperature * (1 - (1 - loss_fraction) * initial_fraction) ** (1 / 3)
    energy_coefficient = (
        3
        * (float(zeta(4)) / float(zeta(3)))
        * constants["boltzmann_constant"]
        / critical_temperature**3
    )
    deposited_energy = scattered * (1 + reabsorption) * constants["e_rec"]
    return float(energy_coefficient * (target_temperature**4 - temperature**4) / deposited_energy)


def accumulated_snr_row(
    *,
    notebook: dict[str, Any],
    constants: dict[str, Any],
    detuning_ghz: float,
    power_mw: float,
    exposure_us: float,
    quantum_efficiency: float,
    read_noise: float,
    loss_fraction: float,
    eta_coll: float,
    t_p: float,
    parameter_set: str,
) -> dict[str, Any]:
    detuning_hz = detuning_ghz * 1e9
    exposure_s = exposure_us * 1e-6
    phase = scalar_phase_peak(detuning_hz, float(constants["column_density"][0]), constants)
    scattered = scattered_photons_per_atom(notebook, constants, detuning_hz, power_mw, exposure_s)
    n_max = clean_loss_budget(scattered, loss_fraction, eta_coll)
    photons = photons_per_object_pixel(notebook, power_mw, exposure_s, quantum_efficiency)
    signal = 2 * t_p * abs(phase) * photons
    shot_noise = t_p * np.sqrt(photons)
    shot_plus_read_noise = np.sqrt(t_p**2 * photons + read_noise**2)
    snr_shot_only = float(signal / shot_noise)
    snr_shot_plus_read = float(signal / shot_plus_read_noise)
    return {
        "parameter_set": parameter_set,
        "Delta_over_2pi_GHz": detuning_ghz,
        "power_mW": power_mw,
        "exposure_us": exposure_us,
        "QE": quantum_efficiency,
        "read_noise_e_rms": read_noise,
        "spatial_observable": "analytical peak object-space camera pixel; pre-NA/PSF",
        "destruction_model": "continuous optimistic clean loss; 30% condensate loss",
        "phase_rad": phase,
        "scattered_photons_per_atom": scattered,
        "N_max_clean_continuous": n_max,
        "PCI_SNR_shot_only_per_frame": snr_shot_only,
        "PCI_SNR_shot_plus_read_per_frame": snr_shot_plus_read,
        "PCI_SNR_total_shot_only": snr_shot_only * np.sqrt(n_max),
        "PCI_SNR_total_shot_plus_read": snr_shot_plus_read * np.sqrt(n_max),
    }


def legacy_table_row(
    *,
    notebook: dict[str, Any],
    constants: dict[str, Any],
    legacy: dict[str, Any],
    detuning_ghz: float,
    snr_exposure_us: float,
    nmax_exposure_us: float,
) -> dict[str, Any]:
    power_mw = float(legacy["probe_power_mw"])
    detuning_hz = detuning_ghz * 1e9
    t_p = float(notebook["pci_recovery"]["phase_plate_transmittance"])
    phase = scalar_phase_peak(detuning_hz, float(constants["column_density"][0]), constants)
    photons = photons_per_object_pixel(notebook, power_mw, snr_exposure_us * 1e-6, 1.0)
    snr = 2 * t_p * phase * np.sqrt(photons)
    scattered = scattered_photons_per_atom(
        notebook, constants, detuning_hz, power_mw, nmax_exposure_us * 1e-6
    )
    n_max = clean_loss_budget(scattered, 0.3, float(notebook["multishot_recovery"]["eta_coll"]))
    return {
        "Delta_over_2pi_GHz": detuning_ghz,
        "power_mW": power_mw,
        "SNR_exposure_us": snr_exposure_us,
        "N_max_exposure_us": nmax_exposure_us,
        "QE": 1.0,
        "read_noise_e_rms": 0.0,
        "spatial_observable": "legacy ideal analytical peak object-space camera pixel",
        "PCI_SNR_per_frame": float(snr),
        "N_max_clean_continuous": n_max,
        "PCI_SNR_total": float(snr * np.sqrt(n_max)),
    }


def budget_row(
    *,
    notebook: dict[str, Any],
    constants: dict[str, Any],
    thermal: dict[str, float],
    parameter_set: str,
    power_mw: float,
    exposure_us: float,
) -> dict[str, Any]:
    detuning_ghz = 1.5
    detuning_hz = detuning_ghz * 1e9
    loss_fraction = 0.3
    scattered = scattered_photons_per_atom(
        notebook, constants, detuning_hz, power_mw, exposure_us * 1e-6
    )
    reabs = reabsorption_fraction(detuning_hz, constants)
    clean = clean_loss_budget(scattered, loss_fraction, float(notebook["multishot_recovery"]["eta_coll"]))
    heating = heating_budget(notebook, constants, thermal, scattered, loss_fraction, 0.0)
    heating_reabs = heating_budget(notebook, constants, thermal, scattered, loss_fraction, reabs)
    return {
        "parameter_set": parameter_set,
        "Delta_over_2pi_GHz": detuning_ghz,
        "power_mW": power_mw,
        "exposure_us": exposure_us,
        "loss_threshold_fraction": loss_fraction,
        "N_max_clean_continuous": clean,
        "N_stop_heating_continuous": heating,
        "N_stop_heating_reabs_continuous": heating_reabs,
        "strict_integer_clean_frames": int(np.floor(clean)),
        "strict_integer_heating_reabs_frames": int(np.floor(heating_reabs)),
        "reabsorption_fraction": reabs,
    }


def operating_point_snr_row(
    notebook: dict[str, Any],
    constants: dict[str, Any],
    pupil: np.ndarray,
    blur_axis: dict[int, float],
    detuning_ghz: float,
    axis: int,
    power_mw: float = 2.0,
    exposure_us: float = 40.0,
) -> dict[str, Any]:
    phase = scalar_phase_peak(detuning_ghz * 1e9, float(constants["column_density"][axis]), constants)
    exposure_s = exposure_us * 1e-6
    ideal_photons = photons_per_object_pixel(notebook, power_mw, exposure_s, 1.0)
    t_p = float(notebook["pci_recovery"]["phase_plate_transmittance"])
    theta = float(notebook["pci_recovery"]["phase_plate_phase_rad"])
    ideal_snr = 2 * t_p * phase * np.sqrt(ideal_photons)
    pixel_snr = pci_snr_pixel_for_phase(
        notebook, phase, power_mw, axis, exposure_s, blur_axis
    )

    ngrid = int(notebook["grid"]["ngrid"])
    fov = float(notebook["grid"]["field_of_view_m"])
    grid_axis = (np.arange(ngrid) - ngrid // 2) * (fov / ngrid)
    ga, gb = np.meshgrid(grid_axis, grid_axis)
    plane = [index for index in range(3) if index != axis]
    profile = np.maximum(
        0,
        1 - ga**2 / constants["radii"][plane[0]] ** 2 - gb**2 / constants["radii"][plane[1]] ** 2,
    ) ** 1.5
    propagated = np.fft.ifft2(np.fft.fft2(np.exp(1j * phase * profile) - 1) * pupil)
    image = np.abs(t_p * np.exp(1j * theta) + propagated) ** 2
    bin_size = int(notebook["camera_recovery"]["bin_size"])
    usable = (ngrid // bin_size) * bin_size
    binned = image[:usable, :usable].reshape(
        usable // bin_size, bin_size, usable // bin_size, bin_size
    ).mean(axis=(1, 3))
    photons = photons_per_camera_pixel(notebook, power_mw, exposure_s)
    read_noise = float(notebook["camera_recovery"]["read_noise_electrons"])
    background = t_p**2
    middle = binned.shape[0] // 2
    block = binned[middle - 1:middle + 2, middle - 1:middle + 2]
    signal = (block - background).sum() * photons
    noise = np.sqrt(block.sum() * photons + block.size * read_noise**2)
    return {
        "Delta_over_2pi_GHz": detuning_ghz,
        "axis": "xyz"[axis],
        "power_mW": power_mw,
        "exposure_us": exposure_us,
        "ideal_QE1_peak_pixel_SNR": float(ideal_snr),
        "realistic_QE0p4_read7_peak_pixel_SNR": float(pixel_snr),
        "realistic_QE0p4_read7_resolution_element_SNR": float(signal / noise),
        "NA_PSF_included_in_realistic_values": True,
    }


def build_audit(contract: dict[str, Any]) -> dict[str, Any]:
    notebook_path = REPO_ROOT / contract["notebook_defaults_config"]
    notebook = load_json(notebook_path)
    constants = _basic_constants(notebook)
    thermal = self_consistent_total_atoms(notebook, constants)
    pupil = build_pupil(notebook)["pupil"]
    blur_axis = build_blur_axis(notebook, constants)
    canonical = contract["canonical_accumulated_snr"]
    accumulated_rows = [
        accumulated_snr_row(
            notebook=notebook,
            constants=constants,
            detuning_ghz=float(detuning),
            power_mw=float(canonical["probe_power_mw"]),
            exposure_us=float(canonical["exposure_time_us"]),
            quantum_efficiency=float(canonical["quantum_efficiency"]),
            read_noise=float(canonical["read_noise_electrons_rms"]),
            loss_fraction=float(canonical["destruction_threshold_fraction"]),
            eta_coll=float(notebook["multishot_recovery"]["eta_coll"]),
            t_p=float(canonical["phase_plate_amplitude_transmittance"]),
            parameter_set="canonical_accumulated_snr",
        )
        for detuning in canonical["detuning_over_2pi_ghz"]
    ]

    legacy_rows = [
        legacy_table_row(
            notebook=notebook,
            constants=constants,
            legacy=contract["legacy_table_5_1_audit"],
            detuning_ghz=float(detuning),
            snr_exposure_us=float(contract["legacy_table_5_1_audit"]["actual_snr_exposure_time_us"]),
            nmax_exposure_us=float(contract["legacy_table_5_1_audit"]["actual_nmax_exposure_time_us"]),
        )
        for detuning in canonical["detuning_over_2pi_ghz"]
    ]
    corrected_legacy_rows = [
        legacy_table_row(
            notebook=notebook,
            constants=constants,
            legacy=contract["legacy_table_5_1_audit"],
            detuning_ghz=float(detuning),
            snr_exposure_us=float(contract["legacy_table_5_1_audit"]["labelled_exposure_time_us"]),
            nmax_exposure_us=float(contract["legacy_table_5_1_audit"]["labelled_exposure_time_us"]),
        )
        for detuning in canonical["detuning_over_2pi_ghz"]
    ]

    budgets = [
        budget_row(
            notebook=notebook,
            constants=constants,
            thermal=thermal,
            parameter_set="legacy_label_claim_P2_tau15",
            power_mw=2.0,
            exposure_us=15.0,
        ),
        budget_row(
            notebook=notebook,
            constants=constants,
            thermal=thermal,
            parameter_set="actual_legacy_calls_P2_tau40",
            power_mw=2.0,
            exposure_us=40.0,
        ),
        budget_row(
            notebook=notebook,
            constants=constants,
            thermal=thermal,
            parameter_set="canonical_route_a_P3p5_tau40",
            power_mw=3.5,
            exposure_us=40.0,
        ),
    ]

    along_rows = []
    along = contract["along_cigar_reference"]
    for detuning_ghz in along["detuning_over_2pi_ghz"]:
        detuning_hz = float(detuning_ghz) * 1e9
        phase = scalar_phase_peak(detuning_hz, float(constants["column_density"][1]), constants)
        scattered = scattered_photons_per_atom(
            notebook,
            constants,
            detuning_hz,
            float(along["probe_power_mw"]),
            float(along["exposure_time_us"]) * 1e-6,
        )
        along_rows.append({
            "Delta_over_2pi_GHz": float(detuning_ghz),
            "power_mW": float(along["probe_power_mw"]),
            "exposure_us": float(along["exposure_time_us"]),
            "axis": "y (along cigar)",
            "peak_phase_rad": phase,
            "phase_regime": "phase-wrapped" if phase >= np.pi else ("linear" if phase < 0.5 else "nonlinear"),
            "N_max_clean_continuous": clean_loss_budget(
                scattered, float(along["destruction_threshold_fraction"]), float(notebook["multishot_recovery"]["eta_coll"])
            ),
        })

    wavelength = float(notebook["atom"]["transition_wavelength_m"])
    resolution_rows = [
        {"NA": na, "Rayleigh_resolution_um": 0.61 * wavelength / na * 1e6}
        for na in (0.08, 0.40)
    ]
    delta_15 = dimensionless_detuning(1.5e9, constants)
    legacy_15 = next(row for row in legacy_rows if row["Delta_over_2pi_GHz"] == 1.5)
    corrected_15 = next(row for row in corrected_legacy_rows if row["Delta_over_2pi_GHz"] == 1.5)
    canonical_15 = next(row for row in accumulated_rows if row["Delta_over_2pi_GHz"] == 1.5)

    corrections = [
        {
            "legacy_claim": "Table 5.1 accumulated SNR = 171.7 at a 40 us operating point",
            "status": "DEPRECATED_MIXED_PARAMETERS",
            "diagnosis": "single-frame SNR used 100 us while N_max used 40 us",
            "legacy_value": legacy_15["PCI_SNR_total"],
            "same_legacy_model_with_consistent_40us": corrected_15["PCI_SNR_total"],
            "canonical_fig3p2_shot_only": canonical_15["PCI_SNR_total_shot_only"],
            "canonical_fig3p2_shot_plus_read": canonical_15["PCI_SNR_total_shot_plus_read"],
        },
        {
            "legacy_claim": "52/25/24 frames at P=2 mW and tau=15 us",
            "status": "DEPRECATED_MISLABELLED_PARAMETERS",
            "diagnosis": "functions used the global 40 us default",
            "actual_parameter_set": "P=2 mW, tau=40 us",
            "correct_15us_parameter_set": "P=2 mW, tau=15 us",
        },
        {
            "legacy_claim": "Route A gives 14 realistic and 30 optimistic usable frames",
            "status": "RETAIN_WITH_COUNTING_QUALIFIER",
            "diagnosis": "14 and 30 are rounded continuous budgets; strict accepted-frame counts are 13 and 29",
        },
    ]

    operating_point_rows = [
        operating_point_snr_row(notebook, constants, pupil, blur_axis, 1.5, 0),
        operating_point_snr_row(notebook, constants, pupil, blur_axis, 13.0, 1),
    ]
    corrections.append({
        "legacy_claim": "operating_point_report prints tau=40 us beside SNR values evaluated at the 100 us camera default",
        "status": "DEPRECATED_MIXED_PARAMETERS",
        "diagnosis": "tau_s=None reaches N_phot_pix default_exposure=100 us while the report label falls back to global tau=40 us",
        "corrected_values_source": "operating_point_snr_rows",
    })

    full_model_path = (
        REPO_ROOT
        / "results"
        / "dissertation_plots_v1"
        / "full_multishot_accumulated_snr"
        / "full_multishot_accumulated_snr_data.csv"
    )
    with full_model_path.open(encoding="utf-8", newline="") as handle:
        full_model_rows = [
            row for row in csv.DictReader(handle) if float(row["detuning_ghz"]) == 1.5
        ]
    if len(full_model_rows) != 4:
        raise ValueError("Expected four full-model PCI/DGI rows at 1.5 GHz")
    expected_frames = budgets[2]["strict_integer_heating_reabs_frames"]
    if any(int(float(row["n_frames_full"])) != expected_frames for row in full_model_rows):
        raise ValueError("Full multishot output does not match the canonical Route A frame count")

    return {
        "contract": contract,
        "constants": {
            "dimensionless_detuning_at_1p5GHz": delta_15,
            "recommended_thesis_rendering": "delta approximately 101.7",
        },
        "accumulated_snr_rows": accumulated_rows,
        "legacy_table_rows": legacy_rows,
        "corrected_legacy_rows": corrected_legacy_rows,
        "budget_rows": budgets,
        "along_cigar_rows": along_rows,
        "resolution_rows": resolution_rows,
        "full_model_rows": full_model_rows,
        "operating_point_snr_rows": operating_point_rows,
        "corrections": corrections,
    }


def fmt(value: float, digits: int = 3) -> str:
    return f"{value:.{digits}f}"


def write_report(path: Path, audit: dict[str, Any]) -> None:
    canonical = audit["accumulated_snr_rows"]
    budgets = {row["parameter_set"]: row for row in audit["budget_rows"]}
    along = audit["along_cigar_rows"]
    legacy_15 = next(row for row in audit["legacy_table_rows"] if row["Delta_over_2pi_GHz"] == 1.5)
    corrected_15 = next(row for row in audit["corrected_legacy_rows"] if row["Delta_over_2pi_GHz"] == 1.5)
    canonical_15 = next(row for row in canonical if row["Delta_over_2pi_GHz"] == 1.5)
    p2_15 = budgets["legacy_label_claim_P2_tau15"]
    p2_40 = budgets["actual_legacy_calls_P2_tau40"]
    route = budgets["canonical_route_a_P3p5_tau40"]
    full_rows = sorted(audit["full_model_rows"], key=lambda row: (row["mode"], row["noise_model"]))
    operating_rows = audit["operating_point_snr_rows"]
    lines = [
        "# Thesis Numerical Consistency: Confirmation and Correction Report",
        "",
        "## Verdict",
        "",
        "Two legacy notebook conclusions mixed parameters and must not be quoted as thesis results:",
        "",
        "1. `171.7` combined a 100 us single-frame SNR with a 40 us clean-loss budget.",
        "2. `52/25/24` were labelled as 15 us results but were evaluated with the global 40 us default.",
        "",
        "The canonical thesis-facing parameter contract is `configs/thesis_numerical_contract_v1.json`. "
        "Every retained number below is regenerated by `scripts/audit_thesis_numerical_consistency.py`.",
        "",
        "## Canonical Parameter Sets",
        "",
        "### A. Fig. 3.2 and supporting accumulated-SNR table",
        "",
        "- `P = 3.5 mW`; `tau = 40 us`; `QE = 0.40`; read noise `7 e- rms`.",
        "- Across-cigar imaging axis (`x`).",
        "- Analytical peak object-space camera pixel before NA/PSF blur; no spatial summation.",
        "- Continuous optimistic clean-loss budget at 30% condensate loss.",
        "- Identical-frame accumulation `SNR_total = SNR_shot sqrt(N_max_clean)`.",
        "",
        "| $\\lvert\\Delta\\rvert/2\\pi$ (GHz) | PCI SNR/frame, shot only | PCI SNR/frame, shot+read | $N_{max}$ clean | Total, shot only | Total, shot+read |",
        "| ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in canonical:
        lines.append(
            f"| {row['Delta_over_2pi_GHz']:.2f} | {row['PCI_SNR_shot_only_per_frame']:.2f} | "
            f"{row['PCI_SNR_shot_plus_read_per_frame']:.2f} | {row['N_max_clean_continuous']:.2f} | "
            f"{row['PCI_SNR_total_shot_only']:.2f} | {row['PCI_SNR_total_shot_plus_read']:.2f} |"
        )
    lines += [
        "",
        "These rows are the only values that should accompany the current idealised Fig. 3.2.",
        "",
        "### B. Quantitative Route A multishot result",
        "",
        "- `|Delta|/2pi = 1.5 GHz`; `P = 3.5 mW`; `tau = 40 us`; 30% condensate-loss threshold.",
        f"- Continuous clean-loss upper bound: `{route['N_max_clean_continuous']:.2f}` pulses.",
        f"- Continuous heating plus reabsorption crossing: `{route['N_stop_heating_reabs_continuous']:.2f}` pulses.",
        f"- Strict integer accepted frames: `{route['strict_integer_heating_reabs_frames']}` realistic, "
        f"`{route['strict_integer_clean_frames']}` clean-loss upper bound.",
        "- Therefore `about 14/30` is allowed only when explicitly called a rounded continuous budget, not an integer frame count.",
        "",
        "The current full evolving calculation uses the same power, exposure, QE, and read noise, but changes both the spatial observable and destruction model:",
        "",
        "| Mode | Noise model | Spatial observable | Strict frames | Full accumulated SNR |",
        "| --- | --- | --- | ---: | ---: |",
    ]
    for row in full_rows:
        lines.append(
            f"| {row['mode']} | {row['noise_model']} | fixed-ROI matched filter | "
            f"{int(float(row['n_frames_full']))} | {float(row['snr_total_full']):.2f} |"
        )
    lines += [
        "",
        "These full-model values support the quantitative evolving-sequence discussion. They are not numerically interchangeable with the peak-pixel Fig. 3.2 curves; the different spatial observable must be stated.",
        "",
        "## Corrections",
        "",
        "### Table 5.1 accumulated SNR",
        "",
        f"- Legacy mixed result at 1.5 GHz: `{legacy_15['PCI_SNR_total']:.2f}` (the printed `171.7`).",
        f"- Same legacy ideal model with 40 us used consistently: `{corrected_15['PCI_SNR_total']:.2f}`.",
        f"- Canonical Fig. 3.2 result: `{canonical_15['PCI_SNR_total_shot_only']:.2f}` shot only and "
        f"`{canonical_15['PCI_SNR_total_shot_plus_read']:.2f}` with read noise.",
        "- Both legacy and canonical analytical quantities are per peak object-space pixel, not per resolution element.",
        "- Power is not the source of the ideal accumulated-SNR difference because its dependence cancels against the clean-loss budget.",
        "",
        "### The 52/25/24 budget statement",
        "",
        f"At the parameters actually used by the calls (`P=2 mW`, `tau=40 us`): clean `{p2_40['N_max_clean_continuous']:.2f}`, "
        f"heating `{p2_40['N_stop_heating_continuous']:.2f}`, heating+reabsorption `{p2_40['N_stop_heating_reabs_continuous']:.2f}`.",
        "",
        f"At the parameters stated by the old label (`P=2 mW`, `tau=15 us`): clean `{p2_15['N_max_clean_continuous']:.2f}`, "
        f"heating `{p2_15['N_stop_heating_continuous']:.2f}`, heating+reabsorption `{p2_15['N_stop_heating_reabs_continuous']:.2f}`.",
        "",
        "Hence `52/25/24 at 15 us` is invalid. Either change the label to 40 us or use the corrected 15 us values.",
        "",
        "### Operating-point SNR printout",
        "",
        "The legacy report printed `tau=40 us` but evaluated its SNR functions at the 100 us camera default. Explicit 40 us results are:",
        "",
        "| $\\lvert\\Delta\\rvert/2\\pi$ (GHz) | Axis | Ideal QE=1 peak pixel | QE=0.4 + read, peak pixel with NA/PSF | QE=0.4 + read, resolution element |",
        "| ---: | --- | ---: | ---: | ---: |",
    ]
    for row in operating_rows:
        lines.append(
            f"| {row['Delta_over_2pi_GHz']:.1f} | {row['axis']} | "
            f"{row['ideal_QE1_peak_pixel_SNR']:.2f} | "
            f"{row['realistic_QE0p4_read7_peak_pixel_SNR']:.2f} | "
            f"{row['realistic_QE0p4_read7_resolution_element_SNR']:.2f} |"
        )
    lines += [
        "",
        "The old `23.86/9.32/18.88` and `57.63/9.70/13.78` rows are therefore not valid 40 us SNR rows.",
        "",
        "## Confirmed Values",
        "",
        f"- Dimensionless detuning at `|Delta|/2pi = 1.5 GHz`: `{audit['constants']['dimensionless_detuning_at_1p5GHz']:.6f}`; write `$\\delta \\approx 101.7$`.",
        f"- Along-cigar peak phase at 1.5 GHz: `{along[0]['peak_phase_rad']:.3f} rad` (phase-wrapped).",
        f"- Along-cigar peak phase at 13 GHz: `{along[1]['peak_phase_rad']:.3f} rad` (linear-limit edge).",
        f"- Along-cigar 13 GHz clean-loss continuous budget at `P=2 mW`, `tau=40 us`: `{along[1]['N_max_clean_continuous']:.0f}` pulses.",
        f"- Rayleigh resolution at NA=0.080: `{audit['resolution_rows'][0]['Rayleigh_resolution_um']:.2f} um`.",
        f"- Rayleigh resolution at NA=0.40: `{audit['resolution_rows'][1]['Rayleigh_resolution_um']:.2f} um`.",
        "",
        "## Mandatory Thesis Labelling Rules",
        "",
        "1. Write GHz detuning as `|Delta|/2pi = X GHz`.",
        "2. Every thesis-figure generation must emit a synchronized provenance table containing quantity, value, detuning, power, exposure, imaging axis, normalisation, N_max model, QE/read noise, and repository path.",
        "3. Attach a parameter-set identifier to every table, figure, and quoted result.",
        "4. State power and exposure together; never rely on a function default in thesis-generating code.",
        "5. State `per pixel`, `per resolution element`, or `matched ROI`; these are not interchangeable.",
        "6. State QE, read noise, and whether NA/PSF filtering is included for every SNR.",
        "7. State `continuous budget` or `strict integer accepted frames` for every frame count.",
        "8. State `clean loss` or `heating plus reabsorption` for every destruction result.",
        "9. Do not cite `171.7` or `52/25/24 at 15 us` as corrected thesis results.",
        "",
        "## Scope Boundary",
        "",
        "The dissertation source containing Sections 3.3, 3.5, 5.1-5.4 and Tables 5.1-5.3 is not present in this repository. "
        "This audit therefore provides the corrected source numbers and labelling contract, but cannot certify the final manuscript wording or table headers until that source is supplied.",
        "",
        "No simulator physics, helper API, notebook, or regression baseline is changed by this audit.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate(contract_path: Path) -> dict[str, Path]:
    contract = load_json(contract_path)
    audit = build_audit(contract)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    outputs = {
        "canonical_table": OUTPUT_DIR / "canonical_accumulated_snr_table.csv",
        "legacy_table": OUTPUT_DIR / "legacy_table_5_1_audit.csv",
        "corrected_legacy_table": OUTPUT_DIR / "corrected_legacy_table_5_1.csv",
        "budgets": OUTPUT_DIR / "destruction_budget_comparison.csv",
        "along_cigar": OUTPUT_DIR / "along_cigar_reference.csv",
        "full_model": OUTPUT_DIR / "full_multishot_reference_at_1p5GHz.csv",
        "operating_point_snr": OUTPUT_DIR / "operating_point_snr_corrections.csv",
        "summary": OUTPUT_DIR / "corrected_numbers.json",
        "metadata": OUTPUT_DIR / "metadata.json",
        "report": REPORT_PATH,
    }
    write_csv(outputs["canonical_table"], audit["accumulated_snr_rows"])
    write_csv(outputs["legacy_table"], audit["legacy_table_rows"])
    write_csv(outputs["corrected_legacy_table"], audit["corrected_legacy_rows"])
    write_csv(outputs["budgets"], audit["budget_rows"])
    write_csv(outputs["along_cigar"], audit["along_cigar_rows"])
    write_csv(outputs["full_model"], audit["full_model_rows"])
    write_csv(outputs["operating_point_snr"], audit["operating_point_snr_rows"])
    write_json(outputs["summary"], {
        "constants": audit["constants"],
        "corrections": audit["corrections"],
        "resolution_rows": audit["resolution_rows"],
        "status": "PASS_WITH_LEGACY_CORRECTIONS",
    })
    write_json(outputs["metadata"], {
        "audit": "thesis numerical parameter consistency",
        "calibration_status": contract["calibration_status"],
        "contract": str(contract_path.resolve().relative_to(REPO_ROOT)),
        "git_branch": git_value("branch", "--show-current"),
        "git_commit": git_value("rev-parse", "HEAD"),
        "notebook_modified": False,
        "simulator_code_modified": False,
        "status": "PASS_WITH_LEGACY_CORRECTIONS",
    })
    write_report(outputs["report"], audit)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONTRACT)
    args = parser.parse_args()
    outputs = generate(args.config.resolve())
    print("Thesis numerical consistency audit passed with legacy corrections.")
    for name, path in outputs.items():
        print(f"{name}: {path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
