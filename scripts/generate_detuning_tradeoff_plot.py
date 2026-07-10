"""Generate a dissertation detuning trade-off physics plot."""

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
    dimensionless_detuning,
    faraday_rotation_angle,
    reabsorption_fraction,
    residual_optical_depth,
    scalar_phase_shift,
    scattered_photons_per_atom,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK_CONFIG = REPO_ROOT / "configs" / "notebook_v1_defaults.json"


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _git_commit() -> str:
    commands = [
        ["git", "rev-parse", "HEAD"],
        [r"C:\Program Files\Git\cmd\git.exe", "rev-parse", "HEAD"],
    ]
    for command in commands:
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


def _repo_relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def _detuning_values_hz(config: dict[str, Any]) -> np.ndarray:
    detuning = config["detuning"]
    if detuning["spacing"] != "log":
        raise ValueError("Only logarithmic detuning spacing is currently supported")
    values_ghz = np.geomspace(
        float(detuning["min_ghz"]),
        float(detuning["max_ghz"]),
        int(detuning["num_points"]),
    )
    if values_ghz.size < 4:
        raise ValueError("detuning.num_points must be at least 4")
    return values_ghz * 1e9


def _notebook_physical_parameters(notebook_config: dict[str, Any], plot_config: dict[str, Any]) -> dict[str, Any]:
    constants = notebook_config["constants"]
    atom = notebook_config["atom"]
    condensate = notebook_config["condensate"]
    geometry = notebook_config["imaging_geometry"]
    multishot = notebook_config["multishot_recovery"]
    faraday = notebook_config["faraday_recovery"]

    atomic_mass = atom["mass_number"] * constants["atomic_mass_unit"]
    scattering_length = condensate["scattering_length_bohr"] * constants["bohr_radius_m"]
    tf_state = build_thomas_fermi_state(
        condensate["atom_number"],
        scattering_length,
        condensate["trap_frequencies_hz"],
        atomic_mass,
        constants["hbar"],
        constants["boltzmann_constant"],
    )

    imaging_axis = int(multishot["imaging_axis"])
    if imaging_axis < 0 or imaging_axis >= len(tf_state.column_density):
        raise ValueError("multishot_recovery.imaging_axis is outside the column-density vector")

    h_planck = 2 * np.pi * constants["hbar"]
    saturation_intensity = (
        np.pi
        * h_planck
        * constants["speed_of_light"]
        * atom["natural_linewidth_rad_s"]
        / (3 * atom["transition_wavelength_m"] ** 3)
    )

    use_peak_intensity = bool(plot_config["fixed_operating_point"].get("use_peak_intensity", True))
    return {
        "species": atom["species"],
        "atom_number": float(condensate["atom_number"]),
        "imaging_axis": imaging_axis,
        "imaging_axis_label": geometry["axis_labels"][imaging_axis],
        "column_density_peak_m2": float(tf_state.column_density[imaging_axis]),
        "column_densities_m2": [float(value) for value in tf_state.column_density],
        "resonant_cross_section_m2": float(atom["resonant_cross_section_m2"]),
        "gamma_rad_per_s": float(atom["natural_linewidth_rad_s"]),
        "transition_wavelength_m": float(atom["transition_wavelength_m"]),
        "probe_diameter_m": float(geometry["probe_diameter_m"]),
        "probe_power_mw": float(multishot["probe_power_mw"]),
        "pulse_duration_s": float(multishot["pulse_duration_us"]) * 1e-6,
        "pulse_duration_us": float(multishot["pulse_duration_us"]),
        "kappa_F": float(faraday["kappa_F"]),
        "saturation_intensity_w_m2": float(saturation_intensity),
        "use_peak_intensity": use_peak_intensity,
    }


def _normalise(values: np.ndarray) -> np.ndarray:
    finite = np.asarray(values, dtype=float)
    scale = np.nanmax(np.abs(finite[np.isfinite(finite)]))
    if not np.isfinite(scale) or scale == 0:
        return np.zeros_like(finite)
    return finite / scale


def build_tradeoff_data(
    notebook_config: dict[str, Any],
    plot_config: dict[str, Any],
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    params = _notebook_physical_parameters(notebook_config, plot_config)
    detuning_hz = _detuning_values_hz(plot_config)

    delta = np.asarray(
        [dimensionless_detuning(float(value), params["gamma_rad_per_s"]) for value in detuning_hz],
        dtype=float,
    )
    scalar_phase = np.asarray(
        [
            scalar_phase_shift(
                float(value),
                params["column_density_peak_m2"],
                params["resonant_cross_section_m2"],
                params["gamma_rad_per_s"],
            )
            for value in detuning_hz
        ],
        dtype=float,
    )
    faraday_rotation = np.asarray(
        [
            faraday_rotation_angle(
                float(value),
                params["column_density_peak_m2"],
                params["resonant_cross_section_m2"],
                params["gamma_rad_per_s"],
                params["kappa_F"],
            )
            for value in detuning_hz
        ],
        dtype=float,
    )
    residual_od = np.asarray(
        [
            residual_optical_depth(
                float(value),
                params["column_density_peak_m2"],
                params["resonant_cross_section_m2"],
                params["gamma_rad_per_s"],
            )
            for value in detuning_hz
        ],
        dtype=float,
    )
    scattered = np.asarray(
        [
            scattered_photons_per_atom(
                float(value),
                params["probe_power_mw"],
                params["pulse_duration_s"],
                params["saturation_intensity_w_m2"],
                params["gamma_rad_per_s"],
                params["probe_diameter_m"],
                use_peak_intensity=params["use_peak_intensity"],
            )
            for value in detuning_hz
        ],
        dtype=float,
    )
    reabsorption = np.asarray(
        [
            reabsorption_fraction(
                float(value),
                np.asarray(params["column_densities_m2"], dtype=float),
                params["resonant_cross_section_m2"],
                params["gamma_rad_per_s"],
            )
            for value in detuning_hz
        ],
        dtype=float,
    )

    destructiveness = scattered * (1 + reabsorption)
    abs_phase = np.abs(scalar_phase)
    signal_per_scattered = np.divide(abs_phase, scattered, out=np.full_like(abs_phase, np.nan), where=scattered > 0)
    signal_to_destruction = np.divide(
        abs_phase,
        destructiveness,
        out=np.full_like(abs_phase, np.nan),
        where=destructiveness > 0,
    )

    dispersive_response = np.abs(delta / (1 + delta**2))
    absorptive_response = 1 / (1 + delta**2)

    data = {
        "detuning_hz": detuning_hz,
        "detuning_ghz": detuning_hz / 1e9,
        "dimensionless_delta": delta,
        "scalar_phase_rad": scalar_phase,
        "abs_scalar_phase_rad": abs_phase,
        "faraday_rotation_rad": faraday_rotation,
        "residual_optical_depth": residual_od,
        "scattered_photons_per_atom": scattered,
        "reabsorption_fraction": reabsorption,
        "destructiveness_metric": destructiveness,
        "signal_per_scattered_photon": signal_per_scattered,
        "signal_to_destruction": signal_to_destruction,
        "two_level_dispersive_response_abs": dispersive_response,
        "two_level_absorptive_response": absorptive_response,
    }
    data.update(
        {
            "normalised_abs_scalar_phase": _normalise(abs_phase),
            "normalised_scattered_photons_per_atom": _normalise(scattered),
            "normalised_residual_optical_depth": _normalise(residual_od),
            "normalised_signal_per_scattered_photon": _normalise(signal_per_scattered),
            "normalised_signal_to_destruction": _normalise(signal_to_destruction),
        }
    )
    return data, params


def _log_slope(x: np.ndarray, y: np.ndarray, fraction: float = 0.35) -> float:
    count = max(4, int(np.ceil(x.size * fraction)))
    tail_x = x[-count:]
    tail_y = y[-count:]
    mask = (tail_x > 0) & (tail_y > 0) & np.isfinite(tail_x) & np.isfinite(tail_y)
    if np.count_nonzero(mask) < 2:
        return float("nan")
    return float(np.polyfit(np.log(tail_x[mask]), np.log(tail_y[mask]), 1)[0])


def write_data_csv(path: Path, data: dict[str, np.ndarray]) -> None:
    fieldnames = list(data)
    row_count = len(data["detuning_hz"])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for index in range(row_count):
            writer.writerow({key: f"{values[index]:.12g}" for key, values in data.items()})


def write_figure(path: Path, plot_config: dict[str, Any], data: dict[str, np.ndarray]) -> None:
    figure_config = plot_config["figure"]
    fig, axis = plt.subplots(
        1,
        1,
        figsize=tuple(figure_config["figsize_inches"]),
        constrained_layout=True,
    )

    x = data["detuning_ghz"]
    axis.loglog(x, data["normalised_abs_scalar_phase"], lw=2.2, label=r"$|\phi|$")
    axis.loglog(
        x,
        data["normalised_scattered_photons_per_atom"],
        lw=2.2,
        label=r"$N_\gamma$",
    )
    axis.loglog(
        x,
        data["normalised_signal_per_scattered_photon"],
        lw=2.2,
        color="C2",
        label=r"$|\phi|/N_\gamma$",
    )
    axis.set_xlabel(r"$|\Delta|/2\pi$ (GHz)")
    axis.set_ylabel("Normalised value")
    axis.grid(alpha=0.25, which="both")
    axis.legend(fontsize=8.5)

    reference_ghz = 1.5
    if x.min() <= reference_ghz <= x.max():
        axis.axvline(reference_ghz, color="0.35", lw=1.0, ls=":", alpha=0.75)
    axis.text(
        reference_ghz * 1.03,
        0.72,
        "1.5 GHz",
        fontsize=8,
        color="0.3",
    )

    fig.suptitle(figure_config["title"], fontsize=12.5)
    fig.savefig(path, format="svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def build_summary(data: dict[str, np.ndarray], params: dict[str, Any], plot_config: dict[str, Any]) -> dict[str, Any]:
    phase_slope = _log_slope(data["detuning_hz"], data["abs_scalar_phase_rad"])
    scattering_slope = _log_slope(data["detuning_hz"], data["scattered_photons_per_atom"])
    residual_od_slope = _log_slope(data["detuning_hz"], data["residual_optical_depth"])
    ratio_slope = _log_slope(data["detuning_hz"], data["signal_per_scattered_photon"])
    best_index = int(np.nanargmax(data["signal_per_scattered_photon"]))
    return {
        "label": plot_config["label"],
        "status": "Version 1 representative / uncalibrated detuning trade-off analysis",
        "detuning_min_ghz": float(data["detuning_ghz"][0]),
        "detuning_max_ghz": float(data["detuning_ghz"][-1]),
        "num_points": int(data["detuning_ghz"].size),
        "fixed_probe_power_mw": params["probe_power_mw"],
        "fixed_pulse_duration_us": params["pulse_duration_us"],
        "imaging_axis": params["imaging_axis"],
        "imaging_axis_label": params["imaging_axis_label"],
        "column_density_peak_m2": params["column_density_peak_m2"],
        "kappa_F": params["kappa_F"],
        "far_detuned_log_slope_abs_phase": phase_slope,
        "far_detuned_log_slope_scattered_photons": scattering_slope,
        "far_detuned_log_slope_residual_od": residual_od_slope,
        "far_detuned_log_slope_signal_per_scattered_photon": ratio_slope,
        "best_signal_per_scattered_photon_index": best_index,
        "best_signal_per_scattered_photon_detuning_ghz": float(data["detuning_ghz"][best_index]),
        "best_signal_per_scattered_photon": float(data["signal_per_scattered_photon"][best_index]),
        "normalisation": plot_config["normalisation"],
        "plotted_quantities": plot_config["plotted_quantities"],
        "figure_layout": "single-panel log-log plot with one shared detuning axis and one shared normalised-value axis",
        "data_note": "CSV retains additional absolute quantities such as residual OD and signal_to_destruction for provenance; the SVG intentionally plots only |phi|, N_gamma, and |phi|/N_gamma.",
        "interpretation_note": (
            "In the far-detuned limit the phase-like signal scales approximately as 1/Delta, "
            "while scattering and residual OD scale approximately as 1/Delta^2."
        ),
    }


def build_metadata(
    config_path: Path,
    notebook_config_path: Path,
    plot_config: dict[str, Any],
    params: dict[str, Any],
    data: dict[str, np.ndarray],
    outputs: dict[str, Path],
) -> dict[str, Any]:
    return {
        "label": plot_config["label"],
        "status": "representative Version 1 analysis, not notebook-aligned figure recovery",
        "git_commit_hash": _git_commit(),
        "config_files_used": {
            "plot_config": _repo_relative(config_path),
            "notebook_defaults": _repo_relative(notebook_config_path),
        },
        "source_helper_functions": [
            "build_thomas_fermi_state",
            "dimensionless_detuning",
            "scalar_phase_shift",
            "residual_optical_depth",
            "scattered_photons_per_atom",
            "reabsorption_fraction",
            "faraday_rotation_angle",
        ],
        "notebook_expressions_used": {
            "dimensionless_detuning": "delta = 2 * Delta_Hz * 2*pi / Gamma_rad_s",
            "scalar_phase": "phi = sigma0 * n_col_peak * delta / (2 * (1 + delta**2))",
            "residual_od": "OD = sigma0 * n_col_peak / (1 + delta**2)",
            "scattering": "N_gamma = (Gamma/2) * s/(1 + s + delta**2) * tau",
            "faraday": "theta_F = kappa_F * phi",
        },
        "detuning_convention": plot_config["detuning"]["convention"],
        "detuning_range": {
            "min_ghz": float(data["detuning_ghz"][0]),
            "max_ghz": float(data["detuning_ghz"][-1]),
            "num_points": int(data["detuning_ghz"].size),
            "spacing": plot_config["detuning"]["spacing"],
        },
        "fixed_physical_parameters": {
            "species": params["species"],
            "atom_number": params["atom_number"],
            "column_density_peak_m2": params["column_density_peak_m2"],
            "column_densities_m2": params["column_densities_m2"],
            "resonant_cross_section_m2": params["resonant_cross_section_m2"],
            "gamma_rad_per_s": params["gamma_rad_per_s"],
            "kappa_F": params["kappa_F"],
        },
        "fixed_optical_parameters": {
            "probe_power_mw": params["probe_power_mw"],
            "pulse_duration_s": params["pulse_duration_s"],
            "probe_diameter_m": params["probe_diameter_m"],
            "saturation_intensity_w_m2": params["saturation_intensity_w_m2"],
            "use_peak_intensity": params["use_peak_intensity"],
        },
        "stage_specific_parameters": {
            "source": plot_config["fixed_operating_point"]["source"],
            "reason": plot_config["fixed_operating_point"]["reason"],
            "probe_power_mw": params["probe_power_mw"],
            "pulse_duration_us": params["pulse_duration_us"],
            "imaging_axis": params["imaging_axis"],
        },
        "units": {
            "detuning_hz": "Hz",
            "detuning_ghz": "GHz",
            "dimensionless_delta": "dimensionless",
            "scalar_phase_rad": "rad",
            "faraday_rotation_rad": "rad",
            "residual_optical_depth": "dimensionless",
            "scattered_photons_per_atom": "photons per atom per shot",
            "signal_per_scattered_photon": "rad per scattered photon",
        },
        "normalisation": plot_config["normalisation"],
        "plotted_quantities": plot_config["plotted_quantities"],
        "calibration_status": plot_config["caveats"]["calibration_status"],
        "prediction_status": plot_config["caveats"]["prediction_status"],
        "model_status": plot_config["caveats"]["model_status"],
        "outputs": {key: _repo_relative(path) for key, path in outputs.items()},
    }


def generate(config_path: Path) -> dict[str, Path]:
    plot_config = _load_json(config_path)
    notebook_config_path = REPO_ROOT / plot_config.get("notebook_defaults_config", str(NOTEBOOK_CONFIG))
    notebook_config = _load_json(notebook_config_path)
    data, params = build_tradeoff_data(notebook_config, plot_config)

    output_dir = REPO_ROOT / plot_config["output_directory"]
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "figure": output_dir / "detuning_tradeoff.svg",
        "data": output_dir / "detuning_tradeoff_data.csv",
        "summary": output_dir / "detuning_tradeoff_summary.json",
        "metadata": output_dir / "metadata.json",
    }

    write_data_csv(outputs["data"], data)
    write_figure(outputs["figure"], plot_config, data)
    _write_json(outputs["summary"], build_summary(data, params, plot_config))
    _write_json(
        outputs["metadata"],
        build_metadata(config_path, notebook_config_path, plot_config, params, data, outputs),
    )
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        required=True,
        type=Path,
        help="Path to the dissertation plot config JSON.",
    )
    args = parser.parse_args()
    outputs = generate(args.config)
    print("Generated detuning trade-off outputs:")
    for key, path in outputs.items():
        print(f"- {key}: {path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
