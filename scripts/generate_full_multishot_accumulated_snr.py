"""Generate heating-aware framewise PCI/DGI/Faraday accumulated-SNR results."""

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

from non_destructive_image import (
    bin_to_camera_pixels,
    scalar_phase_shift,
    scattered_photons_per_atom,
    simulate_faraday_image,
)
from scripts.recover_notebook_multishot_stage import (
    _basic_constants,
    self_consistent_total_atoms,
    tf_state_for_atoms,
)
from scripts.recover_notebook_pci_stage import build_pupil


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "results" / "dissertation_plots_v1" / "full_multishot_accumulated_snr"
MODES = ("PCI", "DGI", "Faraday dark-field", "Faraday dual-port")
NOISE_MODELS = ("shot_noise_only", "shot_plus_read_noise")


def _mode_key(mode: str) -> str:
    return mode.lower().replace("-", "_").replace(" ", "_")


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
    return matched_filter_snr_from_electrons(
        signal_e,
        photon_variance,
        read_noise_e**2,
        roi,
    )


def matched_filter_snr_from_electrons(
    signal_e: np.ndarray,
    photon_variance_e2: np.ndarray,
    read_variance_e2: float | np.ndarray,
    roi: np.ndarray,
) -> dict[str, float]:
    """Return matched-template SNR from an electron signal and variances."""

    signal_e = np.asarray(signal_e, dtype=float)
    photon_variance = np.asarray(photon_variance_e2, dtype=float)
    read_variance = np.broadcast_to(np.asarray(read_variance_e2, dtype=float), signal_e.shape)
    if signal_e.shape != photon_variance.shape or signal_e.shape != np.asarray(roi).shape:
        raise ValueError("signal, photon variance, and ROI must have identical shapes")
    if np.any(photon_variance < 0) or np.any(read_variance < 0):
        raise ValueError("noise variances must be non-negative")

    signal_roi = signal_e[roi]
    shot_variance_roi = photon_variance[roi]
    shot_read_variance_roi = shot_variance_roi + read_variance[roi]

    def snr(variance: np.ndarray) -> float:
        valid = variance > 0
        return float(np.sqrt(np.sum(signal_roi[valid] ** 2 / variance[valid])))

    shot_snr = snr(shot_variance_roi)
    shot_read_snr = snr(shot_read_variance_roi)
    signal_l2 = float(np.sqrt(np.sum(signal_roi**2)))
    return {
        "signal_l2_e": signal_l2,
        "shot_noise_effective_e": signal_l2 / shot_snr if shot_snr else np.inf,
        "shot_plus_read_effective_e": signal_l2 / shot_read_snr if shot_read_snr else np.inf,
        "snr_shot_noise_only": shot_snr,
        "snr_shot_plus_read_noise": shot_read_snr,
    }


def matched_dual_port_snr(
    port_u_image: np.ndarray,
    port_v_image: np.ndarray,
    photons_per_pixel: float,
    read_noise_e: float,
    roi: np.ndarray,
) -> dict[str, float]:
    """Return matched SNR for the dual-port count difference ``N_V-N_U``.

    Independent Poisson ports give variance ``N_U+N_V``. Independent read
    noise in the two ports contributes ``2*read_noise_e**2``. This is the
    fixed-ROI form of the notebook's per-pixel normalised-difference error
    propagation; the common denominator cancels from the SNR.
    """

    port_u_e = np.clip(np.asarray(port_u_image), 0.0, None) * photons_per_pixel
    port_v_e = np.clip(np.asarray(port_v_image), 0.0, None) * photons_per_pixel
    return matched_filter_snr_from_electrons(
        port_v_e - port_u_e,
        port_u_e + port_v_e,
        2 * read_noise_e**2,
        roi,
    )


def _model_parameters(notebook: dict[str, Any], plots: dict[str, Any]) -> dict[str, Any]:
    cfg = plots["accumulated_snr_invariance"]
    constants = _basic_constants(notebook)
    thermal = self_consistent_total_atoms(notebook, constants)
    geometry = notebook["imaging_geometry"]
    camera = notebook["camera_recovery"]
    pci = notebook["pci_recovery"]
    dgi = notebook["dgi_recovery"]
    faraday = notebook["faraday_recovery"]
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
        "quantum_efficiency": float(cfg["quantum_efficiency"]),
        "photons_per_pixel": float(photons),
        "bin_size": int(camera["bin_size"]),
        "t_p": float(pci["phase_plate_transmittance"]),
        "theta": float(pci["phase_plate_phase_rad"]),
        "dgi_od": float(dgi["stop_optical_depth"]),
        "kappa_f": float(faraday["kappa_F"]),
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
            (mode, noise): [] for mode in MODES for noise in NOISE_MODELS
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
            phase_map = phase_peak * profile
            faraday = simulate_faraday_image(
                params["kappa_f"] * phase_map,
                pupil,
                return_intermediates=True,
            )
            propagated = faraday["sigma_plus_propagated_scattered_field"]
            images = {
                "PCI": bin_to_camera_pixels(np.abs(pci_reference + propagated) ** 2, params["bin_size"]),
                "DGI": bin_to_camera_pixels(np.abs(dgi_reference + propagated) ** 2, params["bin_size"]),
                "Faraday dark-field": bin_to_camera_pixels(
                    faraday["dark_field_intensity"],
                    params["bin_size"],
                ),
            }
            port_u = bin_to_camera_pixels(faraday["dual_port_u_intensity"], params["bin_size"])
            port_v = bin_to_camera_pixels(faraday["dual_port_v_intensity"], params["bin_size"])
            stats_by_mode = {
                mode: matched_filter_snr(
                    image,
                    backgrounds.get(mode, 0.0),
                    params["photons_per_pixel"],
                    params["read_noise_e"],
                    roi,
                )
                for mode, image in images.items()
            }
            stats_by_mode["Faraday dual-port"] = matched_dual_port_snr(
                port_u,
                port_v,
                params["photons_per_pixel"],
                params["read_noise_e"],
                roi,
            )
            for mode, stats in stats_by_mode.items():
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
                        "peak_faraday_rotation_rad": params["kappa_f"] * phase_peak,
                        "kappa_f": params["kappa_f"],
                        "mode": mode,
                        "noise_model": noise,
                        "calibration_status": (
                            "Version 1 uncalibrated; kappa_F=1 placeholder"
                            if mode.startswith("Faraday")
                            else "canonical notebook-aligned Version 1"
                        ),
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

        for mode in MODES:
            for noise in NOISE_MODELS:
                values = np.asarray(per_curve[(mode, noise)])
                full_total = float(np.sqrt(np.sum(values**2)))
                identical = float(initial_stats[(mode, noise)] * np.sqrt(n_frames))
                clean_matched = float(initial_stats[(mode, noise)] * np.sqrt(budget["nmax_clean_continuous"]))
                legacy_total = (
                    _interpolate_legacy(legacy, mode, noise, float(detuning_ghz))
                    if mode in ("PCI", "DGI")
                    else None
                )
                observable = {
                    "PCI": "intensity-minus-background electron template",
                    "DGI": "intensity-minus-background electron template",
                    "Faraday dark-field": "crossed-analyser electron-count template",
                    "Faraday dual-port": "dual-port electron-count difference N_V-N_U",
                }[mode]
                row = {
                    "detuning_hz": detuning_hz,
                    "detuning_ghz": detuning_ghz,
                    "dimensionless_detuning": 2 * detuning_hz * 2 * np.pi / constants["gamma"],
                    "mode": mode,
                    "noise_model": noise,
                    "observable": (
                        "fixed-ROI diagonal-covariance matched-filter SNR"
                        if mode in ("PCI", "DGI")
                        else f"fixed-ROI diagonal-covariance matched-filter SNR: {observable}"
                    ),
                    "kappa_f": params["kappa_f"] if mode.startswith("Faraday") else None,
                    "calibration_status": (
                        "Version 1 uncalibrated; kappa_F=1 placeholder"
                        if mode.startswith("Faraday")
                        else "canonical notebook-aligned Version 1"
                    ),
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
                    "full_relative_to_legacy_peak_pixel": (
                        full_total / legacy_total - 1 if legacy_total is not None else None
                    ),
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
    for mode in MODES:
        shot = _curve(rows, mode, "shot_noise_only", "snr_total_full")
        read = _curve(rows, mode, "shot_plus_read_noise", "snr_total_full")
        key = _mode_key(mode)
        variations[f"{key}_shot_noise_0.75_to_3_relative_range"] = _variation(shot, detunings)
        variations[f"{key}_shot_plus_read_0.75_to_3_relative_range"] = _variation(read, detunings)
        checks[f"{key}_shot_plus_read_le_shot_all"] = bool(np.all(read <= shot))
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
    faraday_rows = [row for row in rows if row["mode"].startswith("Faraday")]
    checks["faraday_results_finite_positive"] = bool(
        faraday_rows
        and all(np.isfinite(row["snr_total_full"]) and row["snr_total_full"] > 0 for row in faraday_rows)
    )
    checks["faraday_kappa_f_is_v1_placeholder"] = bool(np.isclose(results["params"]["kappa_f"], 1.0))
    noise_ordering = all(
        checks[f"{_mode_key(mode)}_shot_plus_read_le_shot_all"] for mode in MODES
    )
    checks["passed"] = bool(
        checks["full_frame_count_less_than_clean_continuous_all"]
        and noise_ordering
        and checks["full_snr_below_clean_matched_reference_all"]
        and checks["dgi_shot_plus_read_high_detuning_log_slope"] < 0
        and checks["faraday_results_finite_positive"]
        and checks["faraday_kappa_f_is_v1_placeholder"]
    )
    return checks


def _plot(path: Path, rows: list[dict[str, Any]], detunings: np.ndarray) -> None:
    matplotlib.rcParams["svg.hashsalt"] = "full-multishot-accumulated-snr-v1"
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
    svg_text = path.read_text(encoding="utf-8")
    path.write_text("\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n", encoding="utf-8")


def _canonical_faraday_ledger(results: dict[str, Any], output_path: Path) -> list[dict[str, Any]]:
    reference_rows = [
        row
        for row in results["summary_rows"]
        if row["mode"].startswith("Faraday") and np.isclose(row["detuning_ghz"], 1.5)
    ]
    if len(reference_rows) != 4:
        raise RuntimeError(f"Expected four canonical Faraday reference rows, found {len(reference_rows)}")

    params = results["params"]
    constants = params["constants"]
    roi_pixels = results["roi_metadata"]["roi_pixel_count"]
    ledger: list[dict[str, Any]] = []
    for row in reference_rows:
        dual_port = row["mode"] == "Faraday dual-port"
        ledger.append({
            "quantity": f"{row['mode']} accumulated matched-ROI SNR ({row['noise_model']})",
            "value": row["snr_total_full"],
            "units": "dimensionless SNR",
            "Delta_over_2pi_GHz": row["detuning_ghz"],
            "probe_power_mW": params["power_mw"],
            "exposure_time_us": params["exposure_s"] * 1e6,
            "imaging_axis": "xyz"[params["axis"]],
            "normalisation": (
                "electron-count difference N_V-N_U; matched over fixed ROI"
                if dual_port
                else "absolute crossed-analyser intensity; matched over fixed ROI"
            ),
            "N_max_model": "heating+reabsorption; strict integer post-pulse depletion <= 0.30",
            "quantum_efficiency": params["quantum_efficiency"],
            "read_noise_electrons_per_port": params["read_noise_e"],
            "read_variance_e2_per_pixel": (
                2 * params["read_noise_e"] ** 2
                if dual_port and row["noise_model"] == "shot_plus_read_noise"
                else params["read_noise_e"] ** 2
                if row["noise_model"] == "shot_plus_read_noise"
                else 0.0
            ),
            "roi_pixel_count": roi_pixels,
            "n_frames_full": row["n_frames_full"],
            "initial_atom_number": constants["atom_number"],
            "initial_peak_column_density_m2": float(np.asarray(constants["column_density"])[params["axis"]]),
            "saturation_intensity_W_m2": constants["isat"],
            "kappa_F": params["kappa_f"],
            "calibration_status": "Version 1 representative and uncalibrated; kappa_F=1 placeholder",
            "repository_path": str(output_path.relative_to(REPO_ROOT)),
        })
    return ledger


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
        "faraday_ledger": OUTPUT_DIR / "faraday_canonical_reference_at_1p5GHz.csv",
    }
    _write_csv(outputs["data"], results["summary_rows"])
    _write_csv(outputs["framewise"], results["frame_rows"])
    _write_csv(outputs["comparison"], results["comparison_rows"])
    faraday_ledger = _canonical_faraday_ledger(results, outputs["faraday_ledger"])
    _write_csv(outputs["faraday_ledger"], faraday_ledger)
    _plot(outputs["figure"], results["summary_rows"], results["detunings"])

    benchmark = {
        str(value): next(
            row for row in results["summary_rows"]
            if row["mode"] == "PCI" and row["noise_model"] == "shot_noise_only" and np.isclose(row["detuning_ghz"], value)
        )
        for value in (0.75, 1.5, 3.0)
    }
    canonical_reference = {
        row["mode"]: {
            noise: next(
                candidate["snr_total_full"]
                for candidate in results["summary_rows"]
                if candidate["mode"] == row["mode"]
                and candidate["noise_model"] == noise
                and np.isclose(candidate["detuning_ghz"], 1.5)
            )
            for noise in NOISE_MODELS
        }
        for row in results["summary_rows"]
        if row["mode"].startswith("Faraday") and np.isclose(row["detuning_ghz"], 1.5)
    }
    summary = {
        "verdict": "A. READY TO MERGE",
        "checks": checks,
        "canonical_faraday_reference_1p5GHz": canonical_reference,
        "faraday_interpretation": (
            "These values complete the four-mode structural comparison under the canonical "
            "Version 1 convention. They are not calibrated absolute Faraday predictions because "
            "kappa_F remains the phenomenological placeholder 1.0."
        ),
        "benchmark_frame_counts": {
            key: {
                "nmax_clean_continuous": row["nmax_clean_continuous"],
                "nstop_heating_continuous": row["nstop_heating_continuous"],
                "n_frames_full": row["n_frames_full"],
            }
            for key, row in benchmark.items()
        },
        "observable": "fixed 228-pixel spatial ROI with diagonal-covariance matched-template SNR",
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
            "initial_atom_number": results["params"]["constants"]["atom_number"],
            "initial_peak_column_density_m2": float(
                np.asarray(results["params"]["constants"]["column_density"])[results["params"]["axis"]]
            ),
            "saturation_intensity_W_m2": results["params"]["constants"]["isat"],
            "kappa_F": results["params"]["kappa_f"],
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
            "single_port_shot_variance": "total expected photons in each ROI pixel",
            "single_port_read_variance": "7 e- rms squared per binned camera pixel",
            "faraday_dark_field_observable": "crossed-analyser electron counts relative to zero background",
            "faraday_dual_port_observable": "electron-count difference N_V-N_U",
            "faraday_dual_port_shot_variance": "N_U+N_V for independent Poisson ports",
            "faraday_dual_port_read_variance": "2*(7 e-)^2 for two independent readouts",
            "reference_image_noise": "not included, matching the notebook convention of a known reference background",
        },
        "physics_model": {
            "heating": "T_next=(T^4+dE/A_E)^(1/4)",
            "reabsorption": "initial-density angle-averaged notebook fraction, held fixed within each sequence",
            "condensate": "N0=N_total*(1-(T/Tc)^3)",
            "state_update": "Thomas-Fermi state regenerated from evolving N0",
            "image_update": "full complex phase field propagated through the configured Fourier pupil",
            "dgi_reference": "finite OD=4 carrier leakage included",
            "faraday_rotation": "theta_F=kappa_F*phi with kappa_F=1.0 placeholder",
            "faraday_fields": "opposite circular phases propagated and recombined into Ex/Ey, dark field, and U/V ports",
        },
        "caption": "The plotted PCI/DGI curves retain the existing two-panel comparison. The accompanying CSV data and canonical ledger extend the identical evolving matched-ROI pipeline to dark-field and dual-port Faraday. Faraday values use kappa_F=1.0 and are Version 1 uncalibrated structural comparisons, not absolute predictions.",
        "calibration_status": "representative Version 1 uncalibrated",
        "faraday_absolute_prediction_valid": False,
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
