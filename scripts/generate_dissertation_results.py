"""Generate configurable dissertation-ready representative result outputs."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from html import escape
from pathlib import Path
from typing import Any

import numpy as np

from non_destructive_image import (
    evaluate_faraday_operating_point,
    summarise_faraday_sweep,
    sweep_faraday_detuning,
    sweep_faraday_exposure_time,
    sweep_faraday_intensity,
)


HELPERS_USED = [
    "evaluate_faraday_operating_point",
    "sweep_faraday_detuning",
    "sweep_faraday_intensity",
    "sweep_faraday_exposure_time",
    "summarise_faraday_sweep",
]


def _load_config(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _git_commit() -> str:
    commands = [
        ["git", "rev-parse", "HEAD"],
        [r"C:\Program Files\Git\cmd\git.exe", "rev-parse", "HEAD"],
    ]
    for command in commands:
        try:
            return subprocess.check_output(
                command,
                text=True,
                stderr=subprocess.DEVNULL,
            ).strip()
        except (OSError, subprocess.CalledProcessError):
            continue
    return "unknown"


def _base_parameters(config: dict[str, Any]) -> dict[str, Any]:
    physical = config["physical_parameters"]
    return {
        "column_density_peak": physical["column_density_peak"],
        "resonant_cross_section": physical["resonant_cross_section"],
        "gamma_rad_per_s": physical["gamma_rad_per_s"],
        "saturation_intensity": physical["saturation_intensity"],
        "probe_diameter_m": physical["probe_diameter_m"],
        "kappa_f": physical["kappa_F"],
        "column_densities_for_reabsorption": physical.get("column_densities_for_reabsorption"),
        "use_peak_intensity": physical.get("use_peak_intensity", True),
        "photons_per_camera_pixel": physical.get("photons_per_camera_pixel"),
    }


def _json_ready(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: _json_ready(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    return value


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(_json_ready(payload), handle, indent=2, sort_keys=True)
        handle.write("\n")


def _write_table(path: Path, result: dict[str, Any], parameter_key: str) -> None:
    metric_keys = [
        "faraday_signal_rad",
        "faraday_signal_scale",
        "scattered_photons_per_atom",
        "reabsorption_fraction",
        "destructiveness_metric",
        "estimated_per_frame_snr",
        "signal_per_scattered_photon",
        "information_per_scattered_photon",
        "signal_to_destruction",
    ]
    fieldnames = [parameter_key, *metric_keys]
    parameter_values = np.asarray(result[parameter_key], dtype=float)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for index, parameter_value in enumerate(parameter_values):
            row = {parameter_key: f"{parameter_value:.12g}"}
            for key in metric_keys:
                values = np.asarray(result[key], dtype=float)
                row[key] = f"{values[index]:.12g}"
            writer.writerow(row)


def _normalise(values: np.ndarray) -> np.ndarray:
    finite_values = values[np.isfinite(values)]
    if finite_values.size == 0:
        return np.zeros_like(values, dtype=float)
    minimum = float(np.min(finite_values))
    maximum = float(np.max(finite_values))
    if maximum == minimum:
        return np.ones_like(values, dtype=float)
    return (values - minimum) / (maximum - minimum)


def _polyline(points: list[tuple[float, float]], colour: str) -> str:
    serialised = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
    return f'<polyline points="{serialised}" fill="none" stroke="{colour}" stroke-width="2.5" />'


def _write_tradeoff_svg(
    path: Path,
    *,
    title: str,
    x_label: str,
    x_values: np.ndarray,
    metric_values: np.ndarray,
    destructiveness_values: np.ndarray,
    metric_label: str,
    label: str,
) -> None:
    width = 820
    height = 520
    left = 82
    right = 32
    top = 64
    bottom = 76
    plot_width = width - left - right
    plot_height = height - top - bottom

    x_min = float(np.min(x_values))
    x_max = float(np.max(x_values))
    x_span = x_max - x_min if x_max != x_min else 1.0
    metric_norm = _normalise(metric_values)
    destructiveness_norm = _normalise(destructiveness_values)

    def point(x_value: float, y_value: float) -> tuple[float, float]:
        x_pos = left + ((x_value - x_min) / x_span) * plot_width
        y_pos = top + (1.0 - y_value) * plot_height
        return x_pos, y_pos

    metric_points = [point(float(x), float(y)) for x, y in zip(x_values, metric_norm)]
    destructiveness_points = [point(float(x), float(y)) for x, y in zip(x_values, destructiveness_norm)]

    escaped_title = escape(title)
    escaped_label = escape(label)
    escaped_x_label = escape(x_label)
    escaped_metric_label = escape(metric_label)

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <title>{escaped_title}</title>
  <desc>Representative Version 1 uncalibrated dissertation placeholder generated from config label {escaped_label}.</desc>
  <rect width="100%" height="100%" fill="white" />
  <text x="{width / 2:.0f}" y="30" text-anchor="middle" font-family="Arial, sans-serif" font-size="18">{escaped_title}</text>
  <text x="{width / 2:.0f}" y="52" text-anchor="middle" font-family="Arial, sans-serif" font-size="12">Representative V1 uncalibrated result; regenerate after closed-loop calibration</text>
  <line x1="{left}" y1="{top + plot_height}" x2="{left + plot_width}" y2="{top + plot_height}" stroke="black" />
  <line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}" stroke="black" />
  <text x="{left}" y="{top + plot_height + 22}" text-anchor="middle" font-family="Arial, sans-serif" font-size="11">{x_min:.3g}</text>
  <text x="{left + plot_width}" y="{top + plot_height + 22}" text-anchor="middle" font-family="Arial, sans-serif" font-size="11">{x_max:.3g}</text>
  <text x="{left + plot_width / 2:.0f}" y="{height - 22}" text-anchor="middle" font-family="Arial, sans-serif" font-size="14">{escaped_x_label}</text>
  <text x="22" y="{top + plot_height / 2:.0f}" text-anchor="middle" font-family="Arial, sans-serif" font-size="14" transform="rotate(-90 22 {top + plot_height / 2:.0f})">Normalised value (0 to 1)</text>
  <text x="{left - 14}" y="{top + plot_height}" text-anchor="end" font-family="Arial, sans-serif" font-size="11">0</text>
  <text x="{left - 14}" y="{top + 4}" text-anchor="end" font-family="Arial, sans-serif" font-size="11">1</text>
  {_polyline(metric_points, "#1f77b4")}
  {_polyline(destructiveness_points, "#d62728")}
  <rect x="{left + 22}" y="{top + 16}" width="300" height="54" fill="white" stroke="#cccccc" />
  <line x1="{left + 36}" y1="{top + 34}" x2="{left + 74}" y2="{top + 34}" stroke="#1f77b4" stroke-width="2.5" />
  <text x="{left + 84}" y="{top + 38}" font-family="Arial, sans-serif" font-size="12">{escaped_metric_label}</text>
  <line x1="{left + 36}" y1="{top + 56}" x2="{left + 74}" y2="{top + 56}" stroke="#d62728" stroke-width="2.5" />
  <text x="{left + 84}" y="{top + 60}" font-family="Arial, sans-serif" font-size="12">destructiveness_metric</text>
</svg>
"""
    path.write_text(svg, encoding="utf-8")


def _generate(config: dict[str, Any], config_path: Path) -> dict[str, Path]:
    output_dir = Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    base = _base_parameters(config)
    fixed = config["fixed_operating_point"]
    sweeps = config["sweeps"]
    metric_key = config["metric_key"]
    label = config["label"]

    reference = evaluate_faraday_operating_point(
        fixed["detuning_hz"],
        base["column_density_peak"],
        base["resonant_cross_section"],
        base["gamma_rad_per_s"],
        fixed["probe_power_mw"],
        fixed["pulse_duration_s"],
        base["saturation_intensity"],
        base["probe_diameter_m"],
        kappa_f=base["kappa_f"],
        column_densities_for_reabsorption=base["column_densities_for_reabsorption"],
        use_peak_intensity=base["use_peak_intensity"],
        photons_per_camera_pixel=base["photons_per_camera_pixel"],
    )

    detuning = sweep_faraday_detuning(
        sweeps["detuning_values_hz"],
        base["column_density_peak"],
        base["resonant_cross_section"],
        base["gamma_rad_per_s"],
        fixed["probe_power_mw"],
        fixed["pulse_duration_s"],
        base["saturation_intensity"],
        base["probe_diameter_m"],
        kappa_f=base["kappa_f"],
        column_densities_for_reabsorption=base["column_densities_for_reabsorption"],
        use_peak_intensity=base["use_peak_intensity"],
        photons_per_camera_pixel=base["photons_per_camera_pixel"],
        objective_key=metric_key,
    )
    intensity = sweep_faraday_intensity(
        sweeps["probe_power_mw_values"],
        fixed["detuning_hz"],
        base["column_density_peak"],
        base["resonant_cross_section"],
        base["gamma_rad_per_s"],
        fixed["pulse_duration_s"],
        base["saturation_intensity"],
        base["probe_diameter_m"],
        kappa_f=base["kappa_f"],
        column_densities_for_reabsorption=base["column_densities_for_reabsorption"],
        use_peak_intensity=base["use_peak_intensity"],
        photons_per_camera_pixel=base["photons_per_camera_pixel"],
        objective_key=metric_key,
    )
    exposure = sweep_faraday_exposure_time(
        sweeps["pulse_duration_s_values"],
        fixed["detuning_hz"],
        base["column_density_peak"],
        base["resonant_cross_section"],
        base["gamma_rad_per_s"],
        fixed["probe_power_mw"],
        base["saturation_intensity"],
        base["probe_diameter_m"],
        kappa_f=base["kappa_f"],
        column_densities_for_reabsorption=base["column_densities_for_reabsorption"],
        use_peak_intensity=base["use_peak_intensity"],
        photons_per_camera_pixel=base["photons_per_camera_pixel"],
        objective_key=metric_key,
    )

    detuning_summary = summarise_faraday_sweep(detuning, metric_key)
    intensity_summary = summarise_faraday_sweep(intensity, metric_key)
    exposure_summary = summarise_faraday_sweep(exposure, metric_key)

    paths = {
        "detuning_table": output_dir / "detuning_sweep_table.csv",
        "intensity_table": output_dir / "intensity_sweep_table.csv",
        "exposure_table": output_dir / "exposure_time_sweep_table.csv",
        "summary": output_dir / "summary.json",
        "metadata": output_dir / "metadata.json",
        "detuning_figure": output_dir / "detuning_tradeoff.svg",
        "intensity_figure": output_dir / "intensity_tradeoff.svg",
        "exposure_figure": output_dir / "exposure_time_tradeoff.svg",
    }

    _write_table(paths["detuning_table"], detuning, "detuning_hz")
    _write_table(paths["intensity_table"], intensity, "probe_power_mw")
    _write_table(paths["exposure_table"], exposure, "pulse_duration_s")

    _write_json(
        paths["summary"],
        {
            "label": label,
            "metric_key": metric_key,
            "reference_operating_point": reference,
            "detuning_sweep": detuning_summary,
            "intensity_sweep": intensity_summary,
            "exposure_time_sweep": exposure_summary,
        },
    )

    metadata = {
        "label": label,
        "status": "Version 1 representative / uncalibrated results",
        "kappa_F": base["kappa_f"],
        "kappa_F_note": "kappa_F is the Version 1 placeholder unless changed in the config.",
        "calibration_status": "No experimental RAI / absorption calibration has yet been applied.",
        "regeneration_note": "Outputs are intended to be regenerated after closed-loop calibration updates parameters.",
        "git_commit_hash": _git_commit(),
        "config_file_used": str(config_path),
        "helper_functions_used": HELPERS_USED,
        "outputs": {key: str(path) for key, path in paths.items()},
    }
    _write_json(paths["metadata"], metadata)

    _write_tradeoff_svg(
        paths["detuning_figure"],
        title="Faraday Detuning Trade-off",
        x_label="detuning_hz",
        x_values=np.asarray(detuning["detuning_hz"], dtype=float),
        metric_values=np.asarray(detuning[metric_key], dtype=float),
        destructiveness_values=np.asarray(detuning["destructiveness_metric"], dtype=float),
        metric_label=metric_key,
        label=label,
    )
    _write_tradeoff_svg(
        paths["intensity_figure"],
        title="Faraday Probe-Power Trade-off",
        x_label="probe_power_mw",
        x_values=np.asarray(intensity["probe_power_mw"], dtype=float),
        metric_values=np.asarray(intensity[metric_key], dtype=float),
        destructiveness_values=np.asarray(intensity["destructiveness_metric"], dtype=float),
        metric_label=metric_key,
        label=label,
    )
    _write_tradeoff_svg(
        paths["exposure_figure"],
        title="Faraday Exposure-Time Trade-off",
        x_label="pulse_duration_s",
        x_values=np.asarray(exposure["pulse_duration_s"], dtype=float),
        metric_values=np.asarray(exposure[metric_key], dtype=float),
        destructiveness_values=np.asarray(exposure["destructiveness_metric"], dtype=float),
        metric_label=metric_key,
        label=label,
    )

    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        required=True,
        type=Path,
        help="Path to a dissertation result generation JSON config.",
    )
    args = parser.parse_args()

    config_path = args.config
    config = _load_config(config_path)
    paths = _generate(config, config_path)
    print("Generated dissertation result outputs:")
    for key, path in paths.items():
        print(f"- {key}: {path}")


if __name__ == "__main__":
    main()
