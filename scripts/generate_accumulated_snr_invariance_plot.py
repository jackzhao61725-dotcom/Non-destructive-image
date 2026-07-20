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


def _git_branch() -> str:
    for command in (
        ["git", "branch", "--show-current"],
        [r"C:\Program Files\Git\cmd\git.exe", "branch", "--show-current"],
    ):
        try:
            return subprocess.check_output(
                command,
                cwd=REPO_ROOT,
                text=True,
                stderr=subprocess.DEVNULL,
            ).strip()
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
    detuning_ghz = np.unique(np.concatenate([
        np.geomspace(
            float(cfg["detuning_min_ghz"]),
            float(cfg["detuning_max_ghz"]),
            int(cfg["num_points"]),
        ),
        np.asarray([1.5]),
    ]))
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
    pci_snr_shot_plus_read = pci_signal_e / pci_noise_e
    pci_snr_shot_only = pci_signal_e / np.sqrt(t_p**2 * detected_photons)

    # Notebook DGI construction: the passed field is the scattered field. At its
    # peak, |exp(i phi)-1|^2 = 4 sin^2(phi/2), retaining the exact phase response.
    dgi_intensity = 4 * np.sin(phase / 2) ** 2
    dgi_signal_e = dgi_intensity * detected_photons
    dgi_snr_shot_plus_read = dgi_signal_e / np.sqrt(dgi_signal_e + read_noise**2)
    dgi_snr_shot_only = np.sqrt(dgi_signal_e)

    sqrt_n_max = np.sqrt(n_max)
    pci_total_shot_plus_read = pci_snr_shot_plus_read * sqrt_n_max
    pci_total_shot_only = pci_snr_shot_only * sqrt_n_max
    dgi_total_shot_plus_read = dgi_snr_shot_plus_read * sqrt_n_max
    dgi_total_shot_only = dgi_snr_shot_only * sqrt_n_max

    read_dominated = dgi_signal_e <= 0.25 * read_noise**2
    high_detuning = np.arange(detuning_ghz.size) >= int(0.75 * detuning_ghz.size)
    asymptotic_read_noise_regime = int(np.count_nonzero(read_dominated)) >= 3
    slope_mask = read_dominated if asymptotic_read_noise_regime else high_detuning
    ideal_mode_ratio = pci_total_shot_only / dgi_total_shot_only
    checks = {
        "quadratic_shot_plus_read_snr_high_detuning_log_slope": _log_slope(
            detuning_ghz, dgi_snr_shot_plus_read, slope_mask
        ),
        "quadratic_shot_noise_only_snr_shot_log_slope": _log_slope(detuning_ghz, dgi_snr_shot_only),
        "n_max_log_slope": _log_slope(detuning_ghz, n_max),
        "pci_shot_plus_read_total_relative_range": _relative_range(pci_total_shot_plus_read),
        "pci_shot_noise_total_relative_range": _relative_range(pci_total_shot_only),
        "quadratic_shot_noise_total_relative_range": _relative_range(dgi_total_shot_only),
        "quadratic_shot_plus_read_total_log_slope": _log_slope(detuning_ghz, dgi_total_shot_plus_read),
        "quadratic_shot_plus_read_total_high_detuning_log_slope": _log_slope(
            detuning_ghz, dgi_total_shot_plus_read, slope_mask
        ),
        "n_max_mode_independent": True,
        "pci_shot_plus_read_le_shot_noise_limit": bool(np.all(pci_total_shot_plus_read <= pci_total_shot_only)),
        "dgi_shot_plus_read_le_shot_noise_limit": bool(np.all(dgi_total_shot_plus_read <= dgi_total_shot_only)),
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
        "high_detuning_fit_basis": (
            "read-dominated points" if asymptotic_read_noise_regime else "highest-detuning quartile; exact read-noise asymptote not reached"
        ),
        "maximum_abs_phase_rad": float(np.max(np.abs(phase))),
        "notebook_linear_phase_regime_satisfied": bool(np.max(np.abs(phase)) < 0.5),
        "full_fourier_camera_pipeline_used": False,
    }
    checks["passed"] = bool(
        -2.2 <= checks["quadratic_shot_plus_read_snr_high_detuning_log_slope"] <= -0.9
        and -1.1 <= checks["quadratic_shot_noise_only_snr_shot_log_slope"] <= -0.9
        and 1.9 <= checks["n_max_log_slope"] <= 2.1
        and checks["pci_shot_plus_read_total_relative_range"] <= 0.01
        and checks["pci_shot_noise_total_relative_range"] <= 0.01
        and checks["quadratic_shot_noise_total_relative_range"] <= 0.03
        and checks["quadratic_shot_plus_read_total_log_slope"] < -0.05
        and -1.2 <= checks["quadratic_shot_plus_read_total_high_detuning_log_slope"] < 0.0
        and checks["pci_shot_plus_read_le_shot_noise_limit"]
        and checks["dgi_shot_plus_read_le_shot_noise_limit"]
        and 1.9 <= checks["pci_to_dgi_shot_noise_total_ratio_median"] <= 2.1
        and checks["pci_to_dgi_shot_noise_ratio_relative_range"] <= 0.04
        and checks["notebook_linear_phase_regime_satisfied"]
    )

    data = {
        "detuning_hz": detuning_hz,
        "detuning_ghz": detuning_ghz,
        "phase_rad": phase,
        "scattered_photons_per_atom": scattered,
        "n_max": n_max,
        "sqrt_n_max": sqrt_n_max,
        "pci_shot_plus_read_snr_shot": pci_snr_shot_plus_read,
        "pci_shot_plus_read_snr_total": pci_total_shot_plus_read,
        "pci_shot_only_snr_shot": pci_snr_shot_only,
        "pci_shot_only_snr_total": pci_total_shot_only,
        "dgi_shot_plus_read_snr_shot": dgi_snr_shot_plus_read,
        "dgi_shot_plus_read_snr_total": dgi_total_shot_plus_read,
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
        "quantum_efficiency_status": cfg.get("quantum_efficiency_status", "unspecified"),
        "read_noise_electrons": read_noise,
        "read_noise_status": cfg.get("read_noise_status", "unspecified"),
        "camera_model": cfg.get("camera_model", "unspecified"),
        "effective_magnification": float(magnification),
        "object_plane_pixel_um": float(object_pixel_m * 1e6),
        "detected_photons_per_pixel_at_unit_intensity": float(detected_photons),
        "eta_coll": eta_coll,
        "imaging_axis": axis,
    }
    return data, params, checks


def _write_csv(
    path: Path,
    data: dict[str, np.ndarray],
    params: dict[str, Any],
    repo_path: str,
) -> None:
    rows = []
    series = [
        ("PCI", "shot_plus_read_noise", data["pci_shot_plus_read_snr_shot"], data["pci_shot_plus_read_snr_total"]),
        ("PCI", "shot_noise_only", data["pci_shot_only_snr_shot"], data["pci_shot_only_snr_total"]),
        ("DGI", "shot_plus_read_noise", data["dgi_shot_plus_read_snr_shot"], data["dgi_shot_plus_read_snr_total"]),
        ("DGI", "shot_noise_only", data["dgi_shot_only_snr_shot"], data["dgi_shot_only_snr_total"]),
    ]
    for mode, noise_model, snr_shot, snr_total in series:
        for index in range(len(data["detuning_hz"])):
            rows.append([
                data["detuning_hz"][index], data["detuning_ghz"][index], mode, noise_model,
                snr_shot[index], data["n_max"][index], data["sqrt_n_max"][index], snr_total[index],
                params["probe_power_mw"], params["exposure_time_us"], params["imaging_axis"],
                "absolute dimensionless SNR; analytical peak object-space camera pixel; pre-NA/PSF; no spatial sum",
                "30% continuous optimistic clean loss; identical-frame accumulation",
                (
                    f"QE={params['quantum_efficiency']}; read=0 e- rms (shot-noise limit)"
                    if noise_model == "shot_noise_only"
                    else f"QE={params['quantum_efficiency']}; read={params['read_noise_electrons']} e- rms screening reference"
                ),
                repo_path,
            ])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "Delta_Hz", "Delta_over_2pi_GHz", "mode", "noise_model", "SNR_shot", "N_max",
            "sqrt_N_max", "SNR_acc", "P_mW", "tau_us", "imaging_axis", "normalisation",
            "N_max_model", "QE_read", "repo_path",
        ])
        writer.writerows(rows)


def _write_parameter_register(
    path: Path,
    data: dict[str, np.ndarray],
    params: dict[str, Any],
    repo_path: str,
) -> None:
    matches = np.flatnonzero(np.isclose(data["detuning_ghz"], 1.5, rtol=0.0, atol=1e-12))
    if len(matches) != 1:
        raise RuntimeError("Fig. 3.2 data must contain exactly one 1.5 GHz reference point")
    index = int(matches[0])
    normalisation = "absolute dimensionless SNR; analytical peak object-space camera pixel; pre-NA/PSF; no spatial sum"
    nmax_model = "30% continuous optimistic clean loss; identical-frame accumulation"
    common = [
        1.5,
        params["probe_power_mw"],
        params["exposure_time_us"],
        f"x ({params['imaging_axis']})",
        normalisation,
        nmax_model,
    ]
    rows = [
        ["PCI accumulated SNR (shot + read)", data["pci_shot_plus_read_snr_total"][index], *common,
         f"QE={params['quantum_efficiency']}; read={params['read_noise_electrons']} e- rms screening reference", repo_path],
        ["PCI accumulated SNR (shot only)", data["pci_shot_only_snr_total"][index], *common,
         f"QE={params['quantum_efficiency']}; read=0 e- rms", repo_path],
        ["DGI accumulated SNR (shot + read)", data["dgi_shot_plus_read_snr_total"][index], *common,
         f"QE={params['quantum_efficiency']}; read={params['read_noise_electrons']} e- rms screening reference", repo_path],
        ["DGI accumulated SNR (shot only)", data["dgi_shot_only_snr_total"][index], *common,
         f"QE={params['quantum_efficiency']}; read=0 e- rms", repo_path],
        ["N_max clean continuous", data["n_max"][index], 1.5, params["probe_power_mw"],
         params["exposure_time_us"], f"x ({params['imaging_axis']})",
         "absolute continuous pulse budget; no integer rounding", nmax_model, "not applicable", repo_path],
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "量", "值", "Δ/2π (GHz)", "P (mW)", "τ (µs)", "imaging axis", "归一化",
            "N_max模型", "QE/read", "repo 路径",
        ])
        writer.writerows(rows)


def _plot(
    path: Path,
    png_path: Path,
    pdf_path: Path,
    data: dict[str, np.ndarray],
    cfg: dict[str, Any],
) -> None:
    plt.style.use("default")
    matplotlib.rcParams["svg.hashsalt"] = "figure-3-2-accumulated-snr-invariance"
    matplotlib.rcParams.update({
        "font.size": 9.5,
        "axes.labelsize": 10.5,
        "legend.fontsize": 8.8,
        "xtick.labelsize": 9.0,
        "ytick.labelsize": 9.0,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    fig, ax = plt.subplots(figsize=tuple(cfg["figure_size_inches"]))
    pci_color = "#0072B2"
    dgi_color = "#D55E00"
    ax.plot(data["detuning_ghz"], data["pci_shot_plus_read_snr_total"], color=pci_color, lw=2.2, label="PCI: shot + read noise")
    ax.plot(data["detuning_ghz"], data["pci_shot_only_snr_total"], color=pci_color, lw=2.0, ls="--", label="PCI: shot-noise limit")
    ax.plot(data["detuning_ghz"], data["dgi_shot_plus_read_snr_total"], color=dgi_color, lw=2.2, label="DGI: shot + read noise")
    ax.plot(data["detuning_ghz"], data["dgi_shot_only_snr_total"], color=dgi_color, lw=2.0, ls="--", label="DGI: shot-noise limit")

    reference_detuning = 1.5
    reference_index = int(np.flatnonzero(np.isclose(data["detuning_ghz"], reference_detuning, rtol=0.0, atol=1e-12))[0])
    reference_series = (
        (data["pci_shot_plus_read_snr_total"], pci_color, -13),
        (data["pci_shot_only_snr_total"], pci_color, 8),
        (data["dgi_shot_plus_read_snr_total"], dgi_color, -13),
        (data["dgi_shot_only_snr_total"], dgi_color, 8),
    )
    ax.axvline(reference_detuning, color="#6B7280", lw=1.0, ls=":", zorder=0)
    for values, color, y_offset in reference_series:
        value = float(values[reference_index])
        ax.scatter(
            [reference_detuning],
            [value],
            s=28,
            marker="o",
            facecolor="white",
            edgecolor=color,
            linewidth=1.3,
            zorder=4,
        )
        ax.annotate(
            f"{value:.1f}",
            (reference_detuning, value),
            xytext=(5, y_offset),
            textcoords="offset points",
            color=color,
            fontsize=8.0,
            ha="left",
            va="center",
        )
    ax.text(
        reference_detuning,
        0.65 * max(
            float(np.max(data["pci_shot_plus_read_snr_total"])),
            float(np.max(data["pci_shot_only_snr_total"])),
        ),
        "1.5 GHz reference",
        color="#4B5563",
        fontsize=7.8,
        rotation=90,
        ha="left",
        va="bottom",
    )

    ax.set_xscale(cfg["x_axis_scale"])
    ax.set_yscale(cfg["y_axis_scale"])
    ax.set_xlabel(r"$|\Delta|/2\pi$ (GHz)")
    ax.set_ylabel(r"Accumulated peak-pixel SNR, $\mathrm{SNR}_{\mathrm{acc}}$")
    ax.set_xlim(float(data["detuning_ghz"][0]), float(data["detuning_ghz"][-1]))
    y_max = max(
        float(np.max(data["pci_shot_plus_read_snr_total"])),
        float(np.max(data["pci_shot_only_snr_total"])),
        float(np.max(data["dgi_shot_plus_read_snr_total"])),
        float(np.max(data["dgi_shot_only_snr_total"])),
    )
    ax.set_ylim(0.0, 1.10 * y_max)
    x_ticks = [0.75, 1.0, 1.5, 2.0, 3.0, 5.0]
    ax.set_xticks(x_ticks, labels=["0.75", "1.0", "1.5", "2.0", "3.0", "5.0"])
    ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.grid(True, axis="y", which="major", color="#D1D5DB", lw=0.7, alpha=0.65)
    ax.legend(frameon=False, loc="lower left", ncols=1, handlelength=2.8)
    fig.tight_layout()
    fig.savefig(path, format="svg", metadata={"Date": None})
    fig.savefig(png_path, format="png", dpi=300, metadata={"Software": "matplotlib"})
    fig.savefig(pdf_path, format="pdf", metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)
    svg_text = path.read_text(encoding="utf-8")
    path.write_text("\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n", encoding="utf-8")


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
        "figure_png": output_dir / "figure_3_2.png",
        "figure_pdf": output_dir / "figure_3_2.pdf",
        "data": output_dir / "accumulated_snr_invariance_data.csv",
        "parameter_register": output_dir / "figure_3_2_parameter_register.csv",
        "summary": output_dir / "accumulated_snr_invariance_summary.json",
        "metadata": output_dir / "metadata.json",
    }
    data_repo_path = str(outputs["data"].relative_to(REPO_ROOT)).replace("\\", "/")
    _write_csv(outputs["data"], data, params, data_repo_path)
    _write_parameter_register(outputs["parameter_register"], data, params, data_repo_path)
    _plot(outputs["figure"], outputs["figure_png"], outputs["figure_pdf"], data, cfg)

    reference_index = int(np.flatnonzero(np.isclose(data["detuning_ghz"], 1.5, rtol=0.0, atol=1e-12))[0])

    summary = {
        "parameters": params,
        "scaling_checks": checks,
        "curve_ranges": {
            "pci_shot_plus_read_snr_total": [float(np.min(data["pci_shot_plus_read_snr_total"])), float(np.max(data["pci_shot_plus_read_snr_total"]))],
            "pci_shot_noise_only_snr_total": [float(np.min(data["pci_shot_only_snr_total"])), float(np.max(data["pci_shot_only_snr_total"]))],
            "dgi_shot_plus_read_snr_total": [float(np.min(data["dgi_shot_plus_read_snr_total"])), float(np.max(data["dgi_shot_plus_read_snr_total"]))],
            "dgi_shot_noise_only_snr_total": [float(np.min(data["dgi_shot_only_snr_total"])), float(np.max(data["dgi_shot_only_snr_total"]))],
            "n_max": [float(np.min(data["n_max"])), float(np.max(data["n_max"]))],
        },
        "reference_point_1p5GHz": {
            "Delta_over_2pi_GHz": 1.5,
            "pci_shot_plus_read_snr_total": float(data["pci_shot_plus_read_snr_total"][reference_index]),
            "pci_shot_noise_only_snr_total": float(data["pci_shot_only_snr_total"][reference_index]),
            "dgi_shot_plus_read_snr_total": float(data["dgi_shot_plus_read_snr_total"][reference_index]),
            "dgi_shot_noise_only_snr_total": float(data["dgi_shot_only_snr_total"][reference_index]),
            "n_max_clean_continuous": float(data["n_max"][reference_index]),
        },
    }
    metadata = {
        "figure_type": "fixed-destruction-budget accumulated-SNR simulation",
        "git_branch": _git_branch(),
        "git_commit": _git_commit(),
        "script": str(Path(__file__).resolve().relative_to(REPO_ROOT)),
        "config": str(config_path.resolve().relative_to(REPO_ROOT)),
        "notebook_defaults_config": str(notebook_path.resolve().relative_to(REPO_ROOT)),
        "parameter_register": str(outputs["parameter_register"].relative_to(REPO_ROOT)),
        "mode_choices": {"linear": "PCI", "quadratic": "DGI"},
        "noise_assumptions": {
            "shot_plus_read_noise": (
                f"Poisson photon noise plus {params['read_noise_electrons']:.6g} e- rms read noise in quadrature; "
                "the numerical read noise is a plausible dissertation screening reference, not a measured value; "
                "this label does not imply the full Fourier/pupil/binning camera pipeline"
            ),
            "shot_noise_only": "Poisson photon noise with read noise disabled",
        },
        "destruction_budget": "30% condensate loss using the notebook clean-loss model",
        "detuning_range_ghz": [float(data["detuning_ghz"][0]), float(data["detuning_ghz"][-1])],
        "axis_scales": {"x": cfg["x_axis_scale"], "y": cfg["y_axis_scale"]},
        "physical_interpretation": "Within the notebook small-phase analytical scaling model, N_max is mode-independent and grows approximately as detuning squared. Both PCI curves and the DGI shot-noise-limit curve are invariant; the DGI curve including the reference read noise decreases with detuning as its quadratic signal approaches the detector-noise floor.",
        "model_scope": {
            "classification": "analytical scaling figure",
            "pci_signal": "notebook linear small-phase contrast 2*t_p*|phi|",
            "dgi_signal": "opaque-stop limit 4*sin(phi/2)^2",
            "camera_terms": "detected photon scale, Poisson variance, and Gaussian read-noise variance only",
            "not_included": [
                "finite-NA Fourier propagation and axis-dependent blur",
                "camera binning and spatial integration",
                "finite OD=4 stop leakage",
                "condensate depletion during the accumulated sequence",
                "heating plus reabsorption destruction model"
            ],
            "destruction_boundary": "N_max uses the notebook optimistic clean-loss model, not its later realistic heating-plus-reabsorption model"
        },
        "ideal_cross_mode_comparison": {
            "coincide": checks["pci_and_dgi_shot_noise_curves_coincide"],
            "pci_to_dgi_ratio_median": checks["pci_to_dgi_shot_noise_total_ratio_median"],
            "ratio_range": checks["pci_to_dgi_shot_noise_total_ratio_range"],
            "interpretation": checks["ideal_cross_mode_offset_interpretation"],
        },
        "same_mode_noise_checks": {
            "pci_shot_plus_read_le_shot_noise_limit": checks["pci_shot_plus_read_le_shot_noise_limit"],
            "dgi_shot_plus_read_le_shot_noise_limit": checks["dgi_shot_plus_read_le_shot_noise_limit"],
        },
        "caption": (
            "Accumulated signal-to-noise in the small-phase analytical model at a fixed 30% clean-loss budget. "
            f"Dashed curves are photon-shot-noise limits; solid curves use the provisional dissertation reference of {params['read_noise_electrons']:.6g} e- rms. "
            "The dotted vertical line and open markers identify the 1.5 GHz reference point. PCI uses a carrier "
            "reference and therefore has an expected factor-of-two ideal prefactor relative to DGI under these SNR "
            "definitions. Both ideal curves are approximately detuning-invariant, while the DGI curve including read "
            "noise decreases as its quadratic signal approaches the detector-noise floor. This scaling figure does not "
            "include the full Fourier/pupil/binning pipeline or the later heating-plus-reabsorption destruction model."
        ),
        "scaling_checks": checks,
        "caveats": [
            "Version 1 representative and uncalibrated.",
            "No experimental absorption/RAI calibration has been applied.",
            "This is not a final operating-point prediction.",
            "The figure uses fixed-destruction-budget simulation logic.",
            "The clean-loss budget is an optimistic notebook bound; it is not the later realistic heating-plus-reabsorption model.",
            "The curves are analytical scaling references, not outputs of the full spatial imaging and camera pipeline.",
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
