"""Generate the fixed-destruction-budget accumulated-SNR comparison."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from non_destructive_image import (
    build_thomas_fermi_state,
    scalar_phase_shift,
    scattered_photons_per_atom,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_NOTEBOOK_CONFIG = REPO_ROOT / "configs" / "notebook_v1_defaults.json"


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _git_commit() -> str:
    for command in (["git", "rev-parse", "HEAD"], [r"C:\Program Files\Git\cmd\git.exe", "rev-parse", "HEAD"]):
        try:
            return subprocess.check_output(command, cwd=REPO_ROOT, text=True, stderr=subprocess.DEVNULL).strip()
        except (OSError, subprocess.CalledProcessError):
            continue
    return "unknown"


def _log_slope(x: np.ndarray, y: np.ndarray, mask: np.ndarray | None = None) -> float:
    valid = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if mask is not None:
        valid &= mask
    if np.count_nonzero(valid) < 3:
        raise ValueError("At least three positive finite points are required for a scaling fit")
    return float(np.polyfit(np.log(x[valid]), np.log(y[valid]), 1)[0])


def _relative_range(values: np.ndarray) -> float:
    mean = float(np.mean(values))
    return float(np.ptp(values) / abs(mean)) if mean else float("inf")


def build_accumulated_snr_data(
    notebook_config: dict[str, Any],
    plot_config: dict[str, Any],
) -> tuple[dict[str, np.ndarray], dict[str, Any], dict[str, Any]]:
    cfg = plot_config["accumulated_snr_invariance"]
    if cfg["destruction_model"] != "clean_loss":
        raise ValueError("This figure currently supports the notebook clean-loss budget only")
    if cfg["spacing"] != "log":
        raise ValueError("accumulated-SNR detuning spacing must be logarithmic")

    constants = notebook_config["constants"]
    atom = notebook_config["atom"]
    condensate = notebook_config["condensate"]
    geometry = notebook_config["imaging_geometry"]
    multishot = notebook_config["multishot_recovery"]

    atomic_mass = atom["mass_number"] * constants["atomic_mass_unit"]
    scattering_length = condensate["scattering_length_bohr"] * constants["bohr_radius_m"]
    state = build_thomas_fermi_state(
        condensate["atom_number"],
        scattering_length,
        condensate["trap_frequencies_hz"],
        atomic_mass,
        constants["hbar"],
        constants["boltzmann_constant"],
    )

    axis = int(cfg["imaging_axis"])
    detuning_ghz = np.geomspace(
        float(cfg["detuning_min_ghz"]),
        float(cfg["detuning_max_ghz"]),
        int(cfg["num_points"]),
    )
    detuning_hz = detuning_ghz * 1e9
    power_mw = float(cfg["probe_power_mw"])
    exposure_s = float(cfg["exposure_time_us"]) * 1e-6
    qe = float(cfg["quantum_efficiency"])
    read_noise = float(cfg["read_noise_electrons"])
    t_p = float(cfg["phase_plate_amplitude_transmittance"])
    loss_fraction = float(cfg["destruction_budget_fraction"])
    eta_coll = float(multishot["eta_coll"])

    h_planck = 2 * np.pi * constants["hbar"]
    saturation_intensity = (
        np.pi * h_planck * constants["speed_of_light"] * atom["natural_linewidth_rad_s"]
        / (3 * atom["transition_wavelength_m"] ** 3)
    )
    probe_intensity = 2 * power_mw * 1e-3 / (np.pi * (geometry["probe_diameter_m"] / 2) ** 2)
    magnification = geometry["focal_length_2_m"] / geometry["focal_length_1_m"]
    object_pixel_m = geometry["camera_pixel_m"] / magnification
    photon_energy = h_planck * constants["speed_of_light"] / atom["transition_wavelength_m"]
    detected_photons = probe_intensity * object_pixel_m**2 * exposure_s * qe / photon_energy

    phase = np.asarray([
        scalar_phase_shift(
            float(value),
            float(state.column_density[axis]),
            atom["resonant_cross_section_m2"],
            atom["natural_linewidth_rad_s"],
        )
        for value in detuning_hz
    ])
    scattered = np.asarray([
        scattered_photons_per_atom(
            float(value),
            power_mw,
            exposure_s,
            saturation_intensity,
            atom["natural_linewidth_rad_s"],
            geometry["probe_diameter_m"],
            use_peak_intensity=True,
        )
        for value in detuning_hz
    ])
    n_max = -np.log1p(-loss_fraction) / (eta_coll * scattered)

    # Notebook section 8.2: linear PCI contrast against the phase-plate background.
    pci_signal_e = 2 * t_p * np.abs(phase) * detected_photons
    pci_noise_e = np.sqrt(t_p**2 * detected_photons + read_noise**2)
    pci_snr_realistic = pci_signal_e / pci_noise_e
    pci_snr_shot_only = pci_signal_e / np.sqrt(t_p**2 * detected_photons)

    # Notebook DGI construction: the passed field is the scattered field. At its
    # peak, |exp(i phi)-1|^2 = 4 sin^2(phi/2), retaining the exact phase response.
    dgi_intensity = 4 * np.sin(phase / 2) ** 2
    dgi_signal_e = dgi_intensity * detected_photons
    dgi_snr_realistic = dgi_signal_e / np.sqrt(dgi_signal_e + read_noise**2)
    dgi_snr_shot_only = np.sqrt(dgi_signal_e)

    sqrt_n_max = np.sqrt(n_max)
    pci_total_realistic = pci_snr_realistic * sqrt_n_max
    pci_total_shot_only = pci_snr_shot_only * sqrt_n_max
    dgi_total_realistic = dgi_snr_realistic * sqrt_n_max
    dgi_total_shot_only = dgi_snr_shot_only * sqrt_n_max

    read_dominated = dgi_signal_e <= 0.25 * read_noise**2
    ideal_mode_ratio = pci_total_shot_only / dgi_total_shot_only
    checks = {
        "quadratic_realistic_snr_shot_high_detuning_log_slope": _log_slope(
            detuning_ghz, dgi_snr_realistic, read_dominated
        ),
        "quadratic_shot_noise_only_snr_shot_log_slope": _log_slope(detuning_ghz, dgi_snr_shot_only),
        "n_max_log_slope": _log_slope(detuning_ghz, n_max),
        "pci_realistic_total_relative_range": _relative_range(pci_total_realistic),
        "pci_shot_noise_total_relative_range": _relative_range(pci_total_shot_only),
        "quadratic_shot_noise_total_relative_range": _relative_range(dgi_total_shot_only),
        "quadratic_realistic_total_log_slope": _log_slope(detuning_ghz, dgi_total_realistic),
        "quadratic_realistic_total_high_detuning_log_slope": _log_slope(
            detuning_ghz, dgi_total_realistic, read_dominated
        ),
        "n_max_mode_independent": True,
        "pci_realistic_le_shot_noise_limit": bool(np.all(pci_total_realistic <= pci_total_shot_only)),
        "dgi_realistic_le_shot_noise_limit": bool(np.all(dgi_total_realistic <= dgi_total_shot_only)),
        "pci_to_dgi_shot_noise_total_ratio_median": float(np.median(ideal_mode_ratio)),
        "pci_to_dgi_shot_noise_total_ratio_range": [
            float(np.min(ideal_mode_ratio)),
            float(np.max(ideal_mode_ratio)),
        ],
        "pci_to_dgi_shot_noise_ratio_relative_range": _relative_range(ideal_mode_ratio),
        "pci_and_dgi_shot_noise_curves_coincide": bool(
            np.allclose(pci_total_shot_only, dgi_total_shot_only, rtol=0.03, atol=0.0)
        ),
        "ideal_cross_mode_offset_interpretation": (
            "Expected mode-dependent prefactor: PCI uses the carrier as a homodyne reference, "
            "giving approximately 2|phi|sqrt(N_ph), while DGI gives approximately "
            "|phi|sqrt(N_ph). This is not a hidden normalisation correction."
        ),
        "read_dominated_points": int(np.count_nonzero(read_dominated)),
    }
    checks["passed"] = bool(
        -2.2 <= checks["quadratic_realistic_snr_shot_high_detuning_log_slope"] <= -1.7
        and -1.1 <= checks["quadratic_shot_noise_only_snr_shot_log_slope"] <= -0.9
        and 1.9 <= checks["n_max_log_slope"] <= 2.1
        and checks["pci_realistic_total_relative_range"] <= 0.01
        and checks["pci_shot_noise_total_relative_range"] <= 0.01
        and checks["quadratic_shot_noise_total_relative_range"] <= 0.03
        and checks["quadratic_realistic_total_log_slope"] < -0.2
        and -1.1 <= checks["quadratic_realistic_total_high_detuning_log_slope"] <= -0.7
        and checks["pci_realistic_le_shot_noise_limit"]
        and checks["dgi_realistic_le_shot_noise_limit"]
        and 1.9 <= checks["pci_to_dgi_shot_noise_total_ratio_median"] <= 2.1
        and checks["pci_to_dgi_shot_noise_ratio_relative_range"] <= 0.04
    )

    data = {
        "detuning_hz": detuning_hz,
        "detuning_ghz": detuning_ghz,
        "phase_rad": phase,
        "scattered_photons_per_atom": scattered,
        "n_max": n_max,
        "sqrt_n_max": sqrt_n_max,
        "pci_realistic_snr_shot": pci_snr_realistic,
        "pci_realistic_snr_total": pci_total_realistic,
        "pci_shot_only_snr_shot": pci_snr_shot_only,
        "pci_shot_only_snr_total": pci_total_shot_only,
        "dgi_realistic_snr_shot": dgi_snr_realistic,
        "dgi_realistic_snr_total": dgi_total_realistic,
        "dgi_shot_only_snr_shot": dgi_snr_shot_only,
        "dgi_shot_only_snr_total": dgi_total_shot_only,
    }
    params = {
        "linear_mode": cfg["linear_mode"],
        "quadratic_mode": cfg["quadratic_mode"],
        "destruction_budget_fraction": loss_fraction,
        "destruction_model": cfg["destruction_model"],
        "probe_power_mw": power_mw,
        "exposure_time_us": float(cfg["exposure_time_us"]),
        "quantum_efficiency": qe,
        "read_noise_electrons": read_noise,
        "detected_photons_per_pixel_at_unit_intensity": float(detected_photons),
        "eta_coll": eta_coll,
        "imaging_axis": axis,
    }
    return data, params, checks


def _write_csv(path: Path, data: dict[str, np.ndarray]) -> None:
    rows = []
    series = [
        ("PCI", "realistic", data["pci_realistic_snr_shot"], data["pci_realistic_snr_total"]),
        ("PCI", "shot_noise_only", data["pci_shot_only_snr_shot"], data["pci_shot_only_snr_total"]),
        ("DGI", "realistic", data["dgi_realistic_snr_shot"], data["dgi_realistic_snr_total"]),
        ("DGI", "shot_noise_only", data["dgi_shot_only_snr_shot"], data["dgi_shot_only_snr_total"]),
    ]
    for mode, noise_model, snr_shot, snr_total in series:
        for index in range(len(data["detuning_hz"])):
            rows.append([
                data["detuning_hz"][index], data["detuning_ghz"][index], mode, noise_model,
                snr_shot[index], data["n_max"][index], data["sqrt_n_max"][index], snr_total[index],
            ])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Delta_Hz", "Delta_over_2pi_GHz", "mode", "noise_model", "SNR_shot", "N_max", "sqrt_N_max", "SNR_total"])
        writer.writerows(rows)


def _plot(path: Path, data: dict[str, np.ndarray], cfg: dict[str, Any]) -> None:
    plt.style.use("default")
    fig, ax = plt.subplots(figsize=tuple(cfg["figure_size_inches"]))
    ax.plot(data["detuning_ghz"], data["pci_realistic_snr_total"], color="#1769aa", lw=2.2, label="PCI (realistic)")
    ax.plot(data["detuning_ghz"], data["pci_shot_only_snr_total"], color="#1769aa", lw=2.0, ls="--", label="PCI (shot-noise limit)")
    ax.plot(data["detuning_ghz"], data["dgi_realistic_snr_total"], color="#3f8c4d", lw=2.2, label="DGI (realistic)")
    ax.plot(data["detuning_ghz"], data["dgi_shot_only_snr_total"], color="#3f8c4d", lw=2.0, ls="--", label="DGI (shot-noise limit)")
    ax.set_xscale(cfg["x_axis_scale"])
    ax.set_yscale(cfg["y_axis_scale"])
    ax.set_xlabel(r"$|\Delta|/2\pi$ (GHz)")
    ax.set_ylabel(r"$\mathrm{SNR}_{\mathrm{total}}=\mathrm{SNR}_{\mathrm{shot}}\sqrt{N_{\max}}$")
    ax.grid(True, which="major", alpha=0.22)
    ax.legend(frameon=False, loc="best")
    ax.margins(x=0.02, y=0.08)
    fig.tight_layout()
    fig.savefig(path, format="svg", metadata={"Date": None})
    plt.close(fig)


def generate(config_path: Path) -> dict[str, Path]:
    plot_config = _load_json(config_path)
    notebook_path = REPO_ROOT / plot_config.get("notebook_defaults_config", str(DEFAULT_NOTEBOOK_CONFIG.relative_to(REPO_ROOT)))
    notebook_config = _load_json(notebook_path)
    data, params, checks = build_accumulated_snr_data(notebook_config, plot_config)
    if not checks["passed"]:
        raise RuntimeError(f"Accumulated-SNR scaling checks failed: {checks}")

    cfg = plot_config["accumulated_snr_invariance"]
    output_dir = REPO_ROOT / cfg["output_directory"]
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "figure": output_dir / cfg["figure_filename"],
        "data": output_dir / "accumulated_snr_invariance_data.csv",
        "summary": output_dir / "accumulated_snr_invariance_summary.json",
        "metadata": output_dir / "metadata.json",
    }
    _write_csv(outputs["data"], data)
    _plot(outputs["figure"], data, cfg)

    summary = {
        "parameters": params,
        "scaling_checks": checks,
        "curve_ranges": {
            "pci_realistic_snr_total": [float(np.min(data["pci_realistic_snr_total"])), float(np.max(data["pci_realistic_snr_total"]))],
            "pci_shot_noise_only_snr_total": [float(np.min(data["pci_shot_only_snr_total"])), float(np.max(data["pci_shot_only_snr_total"]))],
            "dgi_realistic_snr_total": [float(np.min(data["dgi_realistic_snr_total"])), float(np.max(data["dgi_realistic_snr_total"]))],
            "dgi_shot_noise_only_snr_total": [float(np.min(data["dgi_shot_only_snr_total"])), float(np.max(data["dgi_shot_only_snr_total"]))],
            "n_max": [float(np.min(data["n_max"])), float(np.max(data["n_max"]))],
        },
    }
    metadata = {
        "figure_type": "fixed-destruction-budget accumulated-SNR simulation",
        "git_branch": "work/accumulated-snr-invariance-plot",
        "git_commit": _git_commit(),
        "script": str(Path(__file__).resolve().relative_to(REPO_ROOT)),
        "config": str(config_path.resolve().relative_to(REPO_ROOT)),
        "notebook_defaults_config": str(notebook_path.resolve().relative_to(REPO_ROOT)),
        "mode_choices": {"linear": "PCI", "quadratic": "DGI"},
        "noise_assumptions": {
            "realistic": "Poisson photon noise plus 7 e- rms read noise in quadrature",
            "shot_noise_only": "Poisson photon noise with read noise disabled",
        },
        "destruction_budget": "30% condensate loss using the notebook clean-loss model",
        "detuning_range_ghz": [float(data["detuning_ghz"][0]), float(data["detuning_ghz"][-1])],
        "axis_scales": {"x": cfg["x_axis_scale"], "y": cfg["y_axis_scale"]},
        "physical_interpretation": "N_max is mode-independent and grows approximately as detuning squared. Both PCI curves and the DGI shot-noise-limit curve are invariant; realistic DGI falls at large detuning because read noise dominates its quadratic signal.",
        "ideal_cross_mode_comparison": {
            "coincide": checks["pci_and_dgi_shot_noise_curves_coincide"],
            "pci_to_dgi_ratio_median": checks["pci_to_dgi_shot_noise_total_ratio_median"],
            "ratio_range": checks["pci_to_dgi_shot_noise_total_ratio_range"],
            "interpretation": checks["ideal_cross_mode_offset_interpretation"],
        },
        "same_mode_noise_checks": {
            "pci_realistic_le_shot_noise_limit": checks["pci_realistic_le_shot_noise_limit"],
            "dgi_realistic_le_shot_noise_limit": checks["dgi_realistic_le_shot_noise_limit"],
        },
        "caption": "Accumulated signal-to-noise at a fixed 30% destruction budget. Shot-noise-limit curves provide the ideal reference for each mode, while solid curves include camera read noise. PCI uses a carrier reference and therefore retains an expected factor-of-two ideal prefactor relative to DGI under the present SNR definitions. Both ideal curves are detuning-invariant; the realistic DGI curve falls at large detuning as read noise acts on its quadratic signal.",
        "scaling_checks": checks,
        "caveats": [
            "Version 1 representative and uncalibrated.",
            "No experimental absorption/RAI calibration has been applied.",
            "This is not a final operating-point prediction.",
            "The figure uses fixed-destruction-budget simulation logic.",
        ],
    }
    _write_json(outputs["summary"], summary)
    _write_json(outputs["metadata"], metadata)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=REPO_ROOT / "configs" / "dissertation_plots_v1.json")
    args = parser.parse_args()
    outputs = generate(args.config if args.config.is_absolute() else REPO_ROOT / args.config)
    for name, path in outputs.items():
        print(f"- {name}: {path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
