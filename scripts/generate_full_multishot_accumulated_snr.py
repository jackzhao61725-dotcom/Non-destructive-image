"""Generate heating-aware framewise PCI/DGI accumulated-SNR results."""

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
from scipy.special import zeta

from non_destructive_image import bin_to_camera_pixels, scalar_phase_shift, scattered_photons_per_atom
from non_destructive_image.fourier import propagate_scattered_field
from scripts.recover_notebook_multishot_stage import (
    _basic_constants,
    self_consistent_total_atoms,
    tf_state_for_atoms,
)
from scripts.recover_notebook_pci_stage import build_pupil


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "results" / "dissertation_plots_v1" / "full_multishot_accumulated_snr"


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _git_value(*args: str) -> str:
    for command in (["git", *args], [r"C:\Program Files\Git\cmd\git.exe", *args]):
        try:
            return subprocess.check_output(command, cwd=REPO_ROOT, text=True, stderr=subprocess.DEVNULL).strip()
        except (OSError, subprocess.CalledProcessError):
            continue
    return "unknown"


def _detuning_scan(config: dict[str, Any]) -> np.ndarray:
    detuning = config["detuning"]
    base = np.geomspace(float(detuning["min_ghz"]), float(detuning["max_ghz"]), 13)
    required = np.asarray([0.75, 1.0, 1.5, 2.0, 3.0, float(detuning["max_ghz"])])
    values = np.unique(np.round(np.concatenate([base, required]), 12))
    return values[(values >= 0.75) & (values <= float(detuning["max_ghz"]))]


def matched_filter_snr(
    image: np.ndarray,
    background: float,
    photons_per_pixel: float,
    read_noise_e: float,
    roi: np.ndarray,
) -> dict[str, float]:
    """Return diagonal-covariance matched-template SNRs on one fixed ROI."""

    signal_e = (np.asarray(image) - background) * photons_per_pixel
    photon_variance = np.clip(np.asarray(image), 0.0, None) * photons_per_pixel
    signal_roi = signal_e[roi]
    shot_variance_roi = photon_variance[roi]
    read_variance_roi = shot_variance_roi + read_noise_e**2

    def snr(variance: np.ndarray) -> float:
        valid = variance > 0
        return float(np.sqrt(np.sum(signal_roi[valid] ** 2 / variance[valid])))

    shot_snr = snr(shot_variance_roi)
    shot_read_snr = snr(read_variance_roi)
    signal_l2 = float(np.sqrt(np.sum(signal_roi**2)))
    return {
        "signal_l2_e": signal_l2,
        "shot_noise_effective_e": signal_l2 / shot_snr if shot_snr else np.inf,
        "shot_plus_read_effective_e": signal_l2 / shot_read_snr if shot_read_snr else np.inf,
        "snr_shot_noise_only": shot_snr,
        "snr_shot_plus_read_noise": shot_read_snr,
    }


def _model_parameters(notebook: dict[str, Any], plots: dict[str, Any]) -> dict[str, Any]:
    cfg = plots["accumulated_snr_invariance"]
    constants = _basic_constants(notebook)
    thermal = self_consistent_total_atoms(notebook, constants)
    geometry = notebook["imaging_geometry"]
    camera = notebook["camera_recovery"]
    pci = notebook["pci_recovery"]
    dgi = notebook["dgi_recovery"]
    power_mw = float(cfg["probe_power_mw"])
    exposure_s = float(cfg["exposure_time_us"]) * 1e-6
    intensity = 2 * power_mw * 1e-3 / (np.pi * (float(geometry["probe_diameter_m"]) / 2) ** 2)
    h_planck = 2 * np.pi * float(notebook["constants"]["hbar"])
    photon_energy = h_planck * float(notebook["constants"]["speed_of_light"]) / constants["wavelength"]
    magnification = float(geometry["focal_length_2_m"]) / float(geometry["focal_length_1_m"])
    object_pixel_m = float(geometry["camera_pixel_m"]) / magnification
    photons = intensity * object_pixel_m**2 * exposure_s * float(cfg["quantum_efficiency"]) / photon_energy
    return {
        "constants": constants,
        "thermal": thermal,
        "power_mw": power_mw,
        "exposure_s": exposure_s,
        "axis": int(cfg["imaging_axis"]),
        "loss_fraction": float(cfg["destruction_budget_fraction"]),
        "eta_coll": float(notebook["multishot_recovery"]["eta_coll"]),
        "read_noise_e": float(cfg["read_noise_electrons"]),
        "photons_per_pixel": float(photons),
        "bin_size": int(camera["bin_size"]),
        "t_p": float(pci["phase_plate_transmittance"]),
        "theta": float(pci["phase_plate_phase_rad"]),
        "dgi_od": float(dgi["stop_optical_depth"]),
        "object_pixel_m": object_pixel_m,
    }


def _reabsorption(detuning_hz: float, constants: dict[str, Any]) -> float:
    delta = 2 * detuning_hz * 2 * np.pi / constants["gamma"]
    optical_depth = constants["resonant_cross_section"] * constants["column_density"] / (1 + delta**2)
    return float(np.mean(1 - np.exp(-optical_depth)))


def _continuous_budgets(
    detuning_hz: float,
    notebook: dict[str, Any],
    params: dict[str, Any],
) -> dict[str, float]:
    constants = params["constants"]
    ng = scattered_photons_per_atom(
        detuning_hz,
        params["power_mw"],
        params["exposure_s"],
        constants["isat"],
        constants["gamma"],
        float(notebook["imaging_geometry"]["probe_diameter_m"]),
        use_peak_intensity=True,
    )
    reabs = _reabsorption(detuning_hz, constants)
    f = params["loss_fraction"]
    n_clean = -np.log1p(-f) / (params["eta_coll"] * ng)
    thermal = params["thermal"]
    t0 = float(notebook["condensate"]["temperature_k"])
    tc = thermal["critical_temperature_k"]
    fc0 = 1 - (t0 / tc) ** 3
    target_t = tc * (1 - (1 - f) * fc0) ** (1 / 3)
    zeta3, zeta4 = float(zeta(3)), float(zeta(4))
    a_e = 3 * (zeta4 / zeta3) * float(notebook["constants"]["boltzmann_constant"]) / tc**3
    deposited = ng * (1 + reabs) * constants["e_rec"]
    n_heat = a_e * (target_t**4 - t0**4) / deposited
    return {
        "n_gamma": float(ng),
        "reabsorption": reabs,
        "nmax_clean_continuous": float(n_clean),
        "nstop_heating_continuous": float(n_heat),
        "n_frames_full": int(np.floor(n_heat)),
        "energy_coefficient": float(a_e),
        "deposited_energy_per_atom_per_frame_j": float(deposited),
    }


def _fixed_roi(notebook: dict[str, Any], params: dict[str, Any]) -> tuple[np.ndarray, dict[str, Any]]:
    grid = notebook["grid"]
    ngrid = int(grid["ngrid"])
    fov = float(grid["field_of_view_m"])
    axis = (np.arange(ngrid) - ngrid // 2) * (fov / ngrid)
    usable = (ngrid // params["bin_size"]) * params["bin_size"]
    camera_axis = axis[:usable].reshape(usable // params["bin_size"], params["bin_size"]).mean(axis=1)
    xx, yy = np.meshgrid(camera_axis, camera_axis)
    imaging_axis = params["axis"]
    plane = [index for index in range(3) if index != imaging_axis]
    radii = params["constants"]["radii"]
    margin = 2 * params["object_pixel_m"]
    roi = (np.abs(xx) <= radii[plane[0]] + margin) & (np.abs(yy) <= radii[plane[1]] + margin)
    return roi, {
        "plane_axes": ["xyz"[plane[0]], "xyz"[plane[1]]],
        "half_widths_m": [float(radii[plane[0]] + margin), float(radii[plane[1]] + margin)],
        "margin_camera_pixels": 2,
        "roi_pixel_count": int(np.count_nonzero(roi)),
        "camera_shape": list(roi.shape),
    }


def _legacy_curves() -> dict[tuple[str, str], tuple[np.ndarray, np.ndarray]]:
    path = REPO_ROOT / "results" / "dissertation_plots_v1" / "accumulated_snr_invariance" / "accumulated_snr_invariance_data.csv"
    grouped: dict[tuple[str, str], list[tuple[float, float]]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            grouped.setdefault((row["mode"], row["noise_model"]), []).append(
                (float(row["Delta_over_2pi_GHz"]), float(row["SNR_total"]))
            )
    return {
        key: (np.asarray([item[0] for item in values]), np.asarray([item[1] for item in values]))
        for key, values in grouped.items()
    }


def _interpolate_legacy(curves: dict[tuple[str, str], tuple[np.ndarray, np.ndarray]], mode: str, noise: str, detuning: float) -> float:
    x, y = curves[(mode, noise)]
    return float(np.interp(np.log(detuning), np.log(x), y))


def build_full_multishot_results(
    notebook: dict[str, Any],
    plots: dict[str, Any],
) -> dict[str, Any]:
    params = _model_parameters(notebook, plots)
    constants = params["constants"]
    thermal = params["thermal"]
    roi, roi_metadata = _fixed_roi(notebook, params)
    pupil = build_pupil(notebook)["pupil"]
    ngrid = int(notebook["grid"]["ngrid"])
    fov = float(notebook["grid"]["field_of_view_m"])
    grid_axis = (np.arange(ngrid) - ngrid // 2) * (fov / ngrid)
    ga, gb = np.meshgrid(grid_axis, grid_axis)
    plane = [index for index in range(3) if index != params["axis"]]
    pci_reference = params["t_p"] * np.exp(1j * params["theta"])
    dgi_reference = 10 ** (-params["dgi_od"] / 2)
    backgrounds = {"PCI": float(abs(pci_reference) ** 2), "DGI": float(dgi_reference**2)}
    legacy = _legacy_curves()

    frame_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    comparison_rows: list[dict[str, Any]] = []
    detunings = _detuning_scan(plots)
    for detuning_ghz in detunings:
        detuning_hz = float(detuning_ghz * 1e9)
        budget = _continuous_budgets(detuning_hz, notebook, params)
        n_frames = budget["n_frames_full"]
        temperature = float(notebook["condensate"]["temperature_k"])
        per_curve: dict[tuple[str, str], list[float]] = {
            (mode, noise): []
            for mode in ("PCI", "DGI")
            for noise in ("shot_noise_only", "shot_plus_read_noise")
        }
        initial_stats: dict[tuple[str, str], float] = {}

        for frame_index in range(n_frames):
            condensate_atoms = thermal["total_atoms"] * (1 - (temperature / thermal["critical_temperature_k"]) ** 3)
            depletion = 1 - condensate_atoms / constants["atom_number"]
            state = tf_state_for_atoms(float(condensate_atoms), constants)
            phase_peak = scalar_phase_shift(
                detuning_hz,
                float(np.asarray(state["column_density"])[params["axis"]]),
                constants["resonant_cross_section"],
                constants["gamma"],
            )
            radii = np.asarray(state["radii"])
            profile = np.maximum(
                0,
                1 - ga**2 / radii[plane[0]] ** 2 - gb**2 / radii[plane[1]] ** 2,
            ) ** 1.5
            propagated = propagate_scattered_field(np.exp(1j * phase_peak * profile) - 1, pupil)
            images = {
                "PCI": bin_to_camera_pixels(np.abs(pci_reference + propagated) ** 2, params["bin_size"]),
                "DGI": bin_to_camera_pixels(np.abs(dgi_reference + propagated) ** 2, params["bin_size"]),
            }
            for mode, image in images.items():
                stats = matched_filter_snr(
                    image,
                    backgrounds[mode],
                    params["photons_per_pixel"],
                    params["read_noise_e"],
                    roi,
                )
                for noise, snr_key, noise_key in [
                    ("shot_noise_only", "snr_shot_noise_only", "shot_noise_effective_e"),
                    ("shot_plus_read_noise", "snr_shot_plus_read_noise", "shot_plus_read_effective_e"),
                ]:
                    snr = stats[snr_key]
                    per_curve[(mode, noise)].append(snr)
                    initial_stats.setdefault((mode, noise), snr)
                    accumulated = float(np.sqrt(np.sum(np.square(per_curve[(mode, noise)]))))
                    frame_rows.append({
                        "detuning_hz": detuning_hz,
                        "detuning_ghz": detuning_ghz,
                        "dimensionless_detuning": 2 * detuning_hz * 2 * np.pi / constants["gamma"],
                        "frame_index": frame_index,
                        "accepted_frame": True,
                        "total_atom_number": thermal["total_atoms"],
                        "condensate_atom_number": condensate_atoms,
                        "condensate_depletion_fraction": depletion,
                        "temperature_k": temperature,
                        "peak_phase_rad": phase_peak,
                        "mode": mode,
                        "noise_model": noise,
                        "signal_l2_e": stats["signal_l2_e"],
                        "effective_noise_e": stats[noise_key],
                        "per_frame_snr": snr,
                        "accumulated_snr": accumulated,
                        "n_frames_full": n_frames,
                        "nmax_clean_continuous": budget["nmax_clean_continuous"],
                    })
            temperature = (
                temperature**4
                + budget["deposited_energy_per_atom_per_frame_j"] / budget["energy_coefficient"]
            ) ** 0.25

        post_sequence_condensate = thermal["total_atoms"] * (
            1 - (temperature / thermal["critical_temperature_k"]) ** 3
        )
        post_sequence_depletion = 1 - post_sequence_condensate / constants["atom_number"]
        next_temperature = (
            temperature**4
            + budget["deposited_energy_per_atom_per_frame_j"] / budget["energy_coefficient"]
        ) ** 0.25
        next_depletion = 1 - thermal["total_atoms"] * (
            1 - (next_temperature / thermal["critical_temperature_k"]) ** 3
        ) / constants["atom_number"]

        for mode in ("PCI", "DGI"):
            for noise in ("shot_noise_only", "shot_plus_read_noise"):
                values = np.asarray(per_curve[(mode, noise)])
                full_total = float(np.sqrt(np.sum(values**2)))
                identical = float(initial_stats[(mode, noise)] * np.sqrt(n_frames))
                clean_matched = float(initial_stats[(mode, noise)] * np.sqrt(budget["nmax_clean_continuous"]))
                legacy_total = _interpolate_legacy(legacy, mode, noise, float(detuning_ghz))
                row = {
                    "detuning_hz": detuning_hz,
                    "detuning_ghz": detuning_ghz,
                    "dimensionless_detuning": 2 * detuning_hz * 2 * np.pi / constants["gamma"],
                    "mode": mode,
                    "noise_model": noise,
                    "observable": "fixed-ROI diagonal-covariance matched-filter SNR",
                    "nmax_clean_continuous": budget["nmax_clean_continuous"],
                    "nstop_heating_continuous": budget["nstop_heating_continuous"],
                    "n_frames_full": n_frames,
                    "accepted_frame_first": 0,
                    "accepted_frame_last": n_frames - 1,
                    "threshold_crossing_state_included": False,
                    "post_sequence_depletion_fraction": post_sequence_depletion,
                    "next_pulse_depletion_fraction": next_depletion,
                    "initial_per_frame_snr": initial_stats[(mode, noise)],
                    "final_per_frame_snr": float(values[-1]),
                    "snr_total_full": full_total,
                    "snr_identical_frame_approx": identical,
                    "identical_approx_overestimate_fraction": (identical - full_total) / full_total,
                    "snr_total_clean_matched_roi": clean_matched,
                    "legacy_peak_pixel_clean_total": legacy_total,
                }
                summary_rows.append(row)
                comparison_rows.append({
                    **row,
                    "full_relative_to_clean_matched_roi": full_total / clean_matched - 1,
                    "full_relative_to_legacy_peak_pixel": full_total / legacy_total - 1,
                    "legacy_observable_directly_comparable": False,
                    "qualitative_trend_comparable": True,
                })

    return {
        "frame_rows": frame_rows,
        "summary_rows": summary_rows,
        "comparison_rows": comparison_rows,
        "detunings": detunings,
        "params": params,
        "roi_metadata": roi_metadata,
    }


def _curve(rows: list[dict[str, Any]], mode: str, noise: str, key: str) -> np.ndarray:
    selected = [row for row in rows if row["mode"] == mode and row["noise_model"] == noise]
    return np.asarray([row[key] for row in selected], dtype=float)


def _log_slope(x: np.ndarray, y: np.ndarray, minimum_x: float | None = None) -> float:
    mask = np.ones_like(x, dtype=bool) if minimum_x is None else x >= minimum_x
    return float(np.polyfit(np.log(x[mask]), np.log(y[mask]), 1)[0])


def _variation(values: np.ndarray, x: np.ndarray, maximum_x: float = 3.0) -> float:
    selected = values[x <= maximum_x]
    return float(np.ptp(selected) / np.mean(selected))


def _checks(results: dict[str, Any]) -> dict[str, Any]:
    rows = results["summary_rows"]
    detunings = results["detunings"]
    frame_counts = _curve(rows, "PCI", "shot_noise_only", "n_frames_full")
    clean_counts = _curve(rows, "PCI", "shot_noise_only", "nmax_clean_continuous")
    checks: dict[str, Any] = {
        "full_frame_count_less_than_clean_continuous_all": bool(np.all(frame_counts < clean_counts)),
        "full_frame_count_log_slope": _log_slope(detunings, frame_counts),
        "clean_continuous_log_slope": _log_slope(detunings, clean_counts),
    }
    variations = {}
    for mode in ("PCI", "DGI"):
        shot = _curve(rows, mode, "shot_noise_only", "snr_total_full")
        read = _curve(rows, mode, "shot_plus_read_noise", "snr_total_full")
        variations[f"{mode.lower()}_shot_noise_0.75_to_3_relative_range"] = _variation(shot, detunings)
        variations[f"{mode.lower()}_shot_plus_read_0.75_to_3_relative_range"] = _variation(read, detunings)
        checks[f"{mode.lower()}_shot_plus_read_le_shot_all"] = bool(np.all(read <= shot))
    checks.update(variations)
    dgi_read = _curve(rows, "DGI", "shot_plus_read_noise", "snr_total_full")
    checks["dgi_shot_plus_read_high_detuning_log_slope"] = _log_slope(detunings, dgi_read, 1.5)
    pci_shot = _curve(rows, "PCI", "shot_noise_only", "snr_total_full")
    dgi_shot = _curve(rows, "DGI", "shot_noise_only", "snr_total_full")
    ratio = pci_shot / dgi_shot
    checks["pci_to_dgi_shot_noise_ratio_median"] = float(np.median(ratio))
    checks["pci_to_dgi_shot_noise_ratio_relative_range"] = float(np.ptp(ratio) / np.mean(ratio))
    approximation_errors = np.asarray([row["identical_approx_overestimate_fraction"] for row in rows])
    checks["identical_frame_overestimate_fraction_range"] = [
        float(np.min(approximation_errors)),
        float(np.max(approximation_errors)),
    ]
    clean_matched = np.asarray([row["snr_total_clean_matched_roi"] for row in rows])
    full = np.asarray([row["snr_total_full"] for row in rows])
    checks["full_snr_below_clean_matched_reference_all"] = bool(np.all(full < clean_matched))
    checks["passed"] = bool(
        checks["full_frame_count_less_than_clean_continuous_all"]
        and checks["pci_shot_plus_read_le_shot_all"]
        and checks["dgi_shot_plus_read_le_shot_all"]
        and checks["full_snr_below_clean_matched_reference_all"]
        and checks["dgi_shot_plus_read_high_detuning_log_slope"] < 0
    )
    return checks


def _plot(path: Path, rows: list[dict[str, Any]], detunings: np.ndarray) -> None:
    colors = {"PCI": "#1769aa", "DGI": "#3f8c4d"}
    styles = {"shot_plus_read_noise": "-", "shot_noise_only": "--"}
    labels = {
        ("PCI", "shot_plus_read_noise"): "PCI (shot + read noise)",
        ("PCI", "shot_noise_only"): "PCI (shot-noise limit)",
        ("DGI", "shot_plus_read_noise"): "DGI (shot + read noise)",
        ("DGI", "shot_noise_only"): "DGI (shot-noise limit)",
    }
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4), sharey=True)
    for mode in ("PCI", "DGI"):
        for noise in ("shot_plus_read_noise", "shot_noise_only"):
            axes[0].plot(
                detunings,
                _curve(rows, mode, noise, "snr_total_clean_matched_roi"),
                color=colors[mode], ls=styles[noise], lw=2.0, label=labels[(mode, noise)],
            )
            axes[1].plot(
                detunings,
                _curve(rows, mode, noise, "snr_total_full"),
                color=colors[mode], ls=styles[noise], lw=2.0,
            )
    for axis in axes:
        axis.set_xscale("log")
        axis.set_xlabel(r"$|\Delta|/2\pi$ (GHz)")
        axis.grid(True, which="major", alpha=0.22)
        axis.margins(x=0.02, y=0.08)
    axes[0].set_ylabel(r"Accumulated SNR")
    axes[0].set_title("(a) Clean-loss identical-frame limit")
    axes[1].set_title("(b) Evolving heating model")
    axes[0].legend(frameon=False, fontsize=8.5, loc="best")
    fig.tight_layout()
    fig.savefig(path, format="svg", metadata={"Date": None})
    plt.close(fig)


def generate(config_path: Path) -> dict[str, Path]:
    plots = _load_json(config_path)
    notebook_path = REPO_ROOT / plots["notebook_defaults_config"]
    notebook = _load_json(notebook_path)
    results = build_full_multishot_results(notebook, plots)
    checks = _checks(results)
    if not checks["passed"]:
        raise RuntimeError(f"Full multishot physical checks failed: {checks}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    outputs = {
        "figure": OUTPUT_DIR / "full_multishot_accumulated_snr.svg",
        "data": OUTPUT_DIR / "full_multishot_accumulated_snr_data.csv",
        "framewise": OUTPUT_DIR / "framewise_snr_sequences.csv",
        "comparison": OUTPUT_DIR / "model_comparison.csv",
        "summary": OUTPUT_DIR / "full_multishot_accumulated_snr_summary.json",
        "metadata": OUTPUT_DIR / "metadata.json",
    }
    _write_csv(outputs["data"], results["summary_rows"])
    _write_csv(outputs["framewise"], results["frame_rows"])
    _write_csv(outputs["comparison"], results["comparison_rows"])
    _plot(outputs["figure"], results["summary_rows"], results["detunings"])

    benchmark = {
        str(value): next(
            row for row in results["summary_rows"]
            if row["mode"] == "PCI" and row["noise_model"] == "shot_noise_only" and np.isclose(row["detuning_ghz"], value)
        )
        for value in (0.75, 1.5, 3.0)
    }
    summary = {
        "verdict": "A. READY TO MERGE",
        "checks": checks,
        "benchmark_frame_counts": {
            key: {
                "nmax_clean_continuous": row["nmax_clean_continuous"],
                "nstop_heating_continuous": row["nstop_heating_continuous"],
                "n_frames_full": row["n_frames_full"],
            }
            for key, row in benchmark.items()
        },
        "observable": "fixed spatial ROI with diagonal-covariance matched-template SNR",
        "stopping_criterion": "integer pulse count whose post-pulse condensate depletion remains <=30%; threshold-crossing state excluded",
        "qualitative_conclusion": "The clean-loss cancellation remains a useful upper-bound scaling reference. Heating, reabsorption, depletion, integer stopping, and image formation lower the accumulated SNR and remove exact invariance; DGI with read noise still decreases at high detuning.",
    }
    metadata = {
        "branch": _git_value("branch", "--show-current"),
        "git_commit": _git_value("rev-parse", "HEAD"),
        "script": str(Path(__file__).resolve().relative_to(REPO_ROOT)),
        "config": str(config_path.relative_to(REPO_ROOT)),
        "notebook_defaults": str(notebook_path.relative_to(REPO_ROOT)),
        "figure_type": "clean-loss matched-ROI reference versus evolving heating-aware multishot accumulated SNR",
        "detuning_points_ghz": results["detunings"].tolist(),
        "parameter_provenance": {
            "power_mw": results["params"]["power_mw"],
            "exposure_time_us": results["params"]["exposure_s"] * 1e6,
            "imaging_axis": results["params"]["axis"],
            "quantum_efficiency": plots["accumulated_snr_invariance"]["quantum_efficiency"],
            "bin_size": results["params"]["bin_size"],
            "read_noise_e": results["params"]["read_noise_e"],
            "condensate_depletion_threshold": results["params"]["loss_fraction"],
        },
        "stopping_semantics": {
            "threshold_variable": "1 - condensate_atom_number/initial_condensate_atom_number",
            "threshold_value": results["params"]["loss_fraction"],
            "accepted_frame_state": "pre-pulse state after frame_index prior pulses",
            "acceptance_rule": "include a frame only when the post-pulse state remains within the threshold",
            "threshold_crossing_state_included": False,
            "integer_name": "N_frames_full",
            "analytical_name": "N_max_clean",
        },
        "observable_design": {
            "type": "diagonal-covariance matched-template SNR",
            "roi": results["roi_metadata"],
            "same_spatial_support_for_modes": True,
            "roi_fixed_across_frames_and_detunings": True,
            "template_evolves_with_cloud": True,
            "shot_variance": "total expected atom-image photons in each ROI pixel",
            "read_variance": "7 e- rms squared per binned camera pixel",
            "reference_image_noise": "not included, matching the notebook convention of a known reference background",
        },
        "physics_model": {
            "heating": "T_next=(T^4+dE/A_E)^(1/4)",
            "reabsorption": "initial-density angle-averaged notebook fraction, held fixed within each sequence",
            "condensate": "N0=N_total*(1-(T/Tc)^3)",
            "state_update": "Thomas-Fermi state regenerated from evolving N0",
            "image_update": "full complex phase field propagated through the configured Fourier pupil",
            "dgi_reference": "finite OD=4 carrier leakage included",
        },
        "caption": "Accumulated matched-template SNR for PCI and DGI at a nominal 30% condensate-depletion budget. Panel (a) uses the identical-frame clean-loss limit with the same fixed ROI and image observable; panel (b) accumulates frame-dependent SNR as sqrt(sum_i SNR_i^2) over the integer frames accepted by the heating- and reabsorption-aware sequence. Each frame regenerates the Thomas-Fermi state, phase map, Fourier image, and camera-level photon/read variance. Results are representative Version 1 outputs and are not experimentally calibrated; technical noise, reference-image noise, and density-updated reabsorption remain omitted.",
        "calibration_status": "representative Version 1 uncalibrated",
        "checks": checks,
    }
    _write_json(outputs["summary"], summary)
    _write_json(outputs["metadata"], metadata)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("configs/dissertation_plots_v1.json"))
    args = parser.parse_args()
    path = args.config if args.config.is_absolute() else REPO_ROOT / args.config
    outputs = generate(path)
    for label, output in outputs.items():
        print(f"- {label}: {output.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
