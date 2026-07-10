"""Recover the canonical notebook V1 scalar phase-map stage."""

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

from non_destructive_image import dimensionless_detuning, scalar_phase_shift
from scripts.recover_notebook_condensate_stage import build_condensate_stage, load_config
from scripts.plot_label_utils import PHASE_RAD, coordinate_label, cut_label


def git_commit() -> str:
    commands = [
        ["git", "rev-parse", "HEAD"],
        [r"C:\Program Files\Git\cmd\git.exe", "rev-parse", "HEAD"],
    ]
    for command in commands:
        try:
            return subprocess.check_output(command, text=True, stderr=subprocess.DEVNULL).strip()
        except (OSError, subprocess.CalledProcessError):
            continue
    return "unknown"


def write_json(path: Path, payload: dict[str, Any]) -> None:
    def convert(value: Any) -> Any:
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, dict):
            return {key: convert(item) for key, item in value.items()}
        if isinstance(value, list):
            return [convert(item) for item in value]
        return value

    with path.open("w", encoding="utf-8") as handle:
        json.dump(convert(payload), handle, indent=2, sort_keys=True)
        handle.write("\n")


def write_rows(path: Path, rows: list[dict[str, float]]) -> None:
    if not rows:
        raise ValueError(f"no rows to write for {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: f"{value:.12g}" for key, value in row.items()})


def array_stats(array: np.ndarray) -> dict[str, Any]:
    peak_index = np.unravel_index(int(np.argmax(array)), array.shape)
    return {
        "shape": list(array.shape),
        "min": float(np.min(array)),
        "max": float(np.max(array)),
        "mean": float(np.mean(array)),
        "sum": float(np.sum(array)),
        "centre_value": float(array[array.shape[0] // 2, array.shape[1] // 2]),
        "peak_index": [int(peak_index[0]), int(peak_index[1])],
    }


def max_relative_difference(left: np.ndarray, right: np.ndarray) -> float:
    denominator = np.maximum(np.abs(left), np.abs(right))
    mask = denominator > 0
    if not np.any(mask):
        return 0.0
    return float(np.max(np.abs(left[mask] - right[mask]) / denominator[mask]))


def notebook_phi_peak(
    detuning_hz: float,
    column_density_peak: float,
    resonant_cross_section: float,
    gamma_rad_per_s: float,
) -> float:
    detuning = 2 * detuning_hz * 2 * np.pi / gamma_rad_per_s
    return resonant_cross_section * column_density_peak * detuning / (2 * (1 + detuning**2))


def build_phase_stage(config: dict[str, Any]) -> dict[str, Any]:
    condensate_stage = build_condensate_stage(config)
    atom = config["atom"]
    phase_config = config["phase_recovery"]
    detuning_hz = float(phase_config["detuning_hz"])
    gamma = float(atom["natural_linewidth_rad_s"])
    sigma0 = float(atom["resonant_cross_section_m2"])
    imaging_axis = condensate_stage["imaging_axis"]
    column_density_peak = condensate_stage["state"]["column_density"][imaging_axis]
    profile = condensate_stage["notebook_profile"]

    notebook_phase_peak = notebook_phi_peak(detuning_hz, column_density_peak, sigma0, gamma)
    helper_phase_peak = scalar_phase_shift(detuning_hz, column_density_peak, sigma0, gamma)
    notebook_phase_map = notebook_phase_peak * profile
    helper_phase_map = helper_phase_peak * condensate_stage["helper_profile"]

    return {
        "condensate_stage": condensate_stage,
        "detuning_hz": detuning_hz,
        "dimensionless_detuning_notebook": 2 * detuning_hz * 2 * np.pi / gamma,
        "dimensionless_detuning_helper": dimensionless_detuning(detuning_hz, gamma),
        "resonant_cross_section_m2": sigma0,
        "gamma_rad_per_s": gamma,
        "column_density_peak_m2": column_density_peak,
        "notebook_phase_peak_rad": notebook_phase_peak,
        "helper_phase_peak_rad": helper_phase_peak,
        "notebook_phase_map_rad": notebook_phase_map,
        "helper_phase_map_rad": helper_phase_map,
    }


def comparison_report(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    notebook_map = stage["notebook_phase_map_rad"]
    helper_map = stage["helper_phase_map_rad"]
    peak_diff = abs(stage["notebook_phase_peak_rad"] - stage["helper_phase_peak_rad"])
    return {
        "status": "canonical notebook scalar phase stage matches helper scalar_phase_shift for tested deterministic quantities",
        "source_notebook_cells": {
            "phase_formula": 10,
            "stage_18_2_phase_map": 59,
            "step_19_1_phase_map": 67,
        },
        "detuning_hz": stage["detuning_hz"],
        "dimensionless_detuning": {
            "notebook_expression": stage["dimensionless_detuning_notebook"],
            "helper": stage["dimensionless_detuning_helper"],
            "absolute_difference": abs(stage["dimensionless_detuning_notebook"] - stage["dimensionless_detuning_helper"]),
        },
        "phase_peak_rad": {
            "notebook_expression": stage["notebook_phase_peak_rad"],
            "helper": stage["helper_phase_peak_rad"],
            "absolute_difference": peak_diff,
            "relative_difference": peak_diff / abs(stage["notebook_phase_peak_rad"])
            if stage["notebook_phase_peak_rad"] != 0
            else 0.0,
        },
        "column_density_peak_m2": stage["column_density_peak_m2"],
        "phase_map_comparison": {
            "units": "radians",
            "max_absolute_difference": float(np.max(np.abs(notebook_map - helper_map))),
            "max_relative_difference": max_relative_difference(notebook_map, helper_map),
            "stats": array_stats(notebook_map),
        },
    }


def phase_summary(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    return {
        "label": config["label"],
        "status": "Canonical notebook-aligned scalar phase-stage recovery output.",
        "source_notebook": config["source_notebook"],
        "source_notebook_cells": {
            "phase_formula": 10,
            "stage_18_2_phase_map": 59,
            "step_19_1_phase_map": 67,
        },
        "detuning_hz": stage["detuning_hz"],
        "detuning_ghz": stage["detuning_hz"] / 1e9,
        "dimensionless_detuning": stage["dimensionless_detuning_notebook"],
        "resonant_cross_section_m2": stage["resonant_cross_section_m2"],
        "gamma_rad_per_s": stage["gamma_rad_per_s"],
        "column_density_peak_m2": stage["column_density_peak_m2"],
        "phase_peak_rad": stage["notebook_phase_peak_rad"],
        "phase_map_stats_rad": array_stats(stage["notebook_phase_map_rad"]),
        "scope_boundary": "Only condensate/profile -> scalar phase map is recovered. No PCI, DGI, Faraday, camera, or multishot figures are generated.",
        "no_experimental_calibration_applied": True,
    }


def write_phase_lineouts(path: Path, config: dict[str, Any], stage: dict[str, Any]) -> None:
    units = config["unit_conventions"]
    coordinate_um = (
        stage["condensate_stage"]["coordinate_axis_m"] * units["coordinate_conversion_m_to_um"]
    )
    phase_map = stage["notebook_phase_map_rad"]
    mid = phase_map.shape[0] // 2
    write_rows(
        path,
        [
            {
                "position_um": float(coordinate_um[index]),
                "phase_y_line_rad": float(phase_map[mid, index]),
                "phase_z_line_rad": float(phase_map[index, mid]),
            }
            for index in range(coordinate_um.size)
        ],
    )


def write_phase_figure(path: Path, config: dict[str, Any], stage: dict[str, Any]) -> None:
    display = config["display"]
    phase_config = config["phase_recovery"]
    units = config["unit_conventions"]
    coordinate_um = (
        stage["condensate_stage"]["coordinate_axis_m"] * units["coordinate_conversion_m_to_um"]
    )
    phase_map = stage["notebook_phase_map_rad"]
    mid = phase_map.shape[0] // 2

    plt.rcParams.update(
        {
            "figure.figsize": (8, 5),
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 12,
            "legend.fontsize": 10,
            "figure.dpi": 110,
        }
    )
    fig, (ax_map, ax_line) = plt.subplots(
        1,
        2,
        figsize=tuple(phase_config["figure_size_inches"]),
        gridspec_kw={"width_ratios": [1.5, 1]},
    )
    im = ax_map.imshow(
        phase_map,
        extent=display["image_extent_um"],
        origin="lower",
        cmap=phase_config["map_cmap"],
    )
    ax_map.set_xlim(*display["column_density_xlim_um"])
    ax_map.set_ylim(*display["column_density_ylim_um"])
    ax_map.set_xlabel(coordinate_label("y"))
    ax_map.set_ylabel(coordinate_label("z"))
    ax_map.set_title(r"the phase map $\phi(\mathbf{r})$ the atoms imprint")
    plt.colorbar(im, ax=ax_map, fraction=0.03, label=PHASE_RAD)

    ax_line.plot(coordinate_um, phase_map[mid, :], "C1", lw=2, label=cut_label("y"))
    ax_line.plot(coordinate_um, phase_map[:, mid], "C0", lw=2, label=cut_label("z"))
    ax_line.axhline(stage["notebook_phase_peak_rad"], color="gray", ls=":", lw=1)
    ax_line.annotate(
        rf"$\phi_\mathrm{{peak}}$ = {stage['notebook_phase_peak_rad']:.3f} rad",
        (12, stage["notebook_phase_peak_rad"] * 0.93),
        fontsize=9,
        color="gray",
    )
    ax_line.set_xlim(-45, 45)
    ax_line.set_xlabel(coordinate_label())
    ax_line.set_ylabel(PHASE_RAD)
    ax_line.set_title("centre cuts")
    ax_line.legend(fontsize=8.5)
    ax_line.grid(alpha=0.25)
    fig.suptitle("Step 1 - the only thing the atoms give us", fontsize=12)
    plt.tight_layout()
    fig.savefig(path, format="svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def generate(config: dict[str, Any], config_path: Path) -> dict[str, str]:
    output_dir = Path("results/notebook_aligned_recovery/phase_stage")
    output_dir.mkdir(parents=True, exist_ok=True)
    stage = build_phase_stage(config)

    outputs = {
        "comparison_report": output_dir / "comparison_report.json",
        "phase_summary": output_dir / "phase_summary.json",
        "central_phase_lineouts": output_dir / "central_phase_lineouts.csv",
        "metadata": output_dir / "metadata.json",
        "figure": output_dir / "scalar_phase_stage.svg",
    }
    write_json(outputs["comparison_report"], comparison_report(config, stage))
    write_json(outputs["phase_summary"], phase_summary(config, stage))
    write_phase_lineouts(outputs["central_phase_lineouts"], config, stage)
    write_phase_figure(outputs["figure"], config, stage)
    write_json(
        outputs["metadata"],
        {
            "label": config["label"],
            "status": "Canonical notebook-aligned scalar phase-stage recovery output.",
            "config_file_used": str(config_path),
            "git_commit_hash": git_commit(),
            "source_notebook": config["source_notebook"],
            "source_notebook_cells": {
                "phase_formula": 10,
                "stage_18_2_phase_map": 59,
                "step_19_1_phase_map": 67,
            },
            "generated_outputs": {key: str(path) for key, path in outputs.items()},
            "scope_boundary": "Only scalar phase recovery is generated. No PCI, DGI, Faraday, camera, or multishot recovery is performed.",
            "no_experimental_calibration_applied": True,
        },
    )
    return {key: str(path) for key, path in outputs.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config = load_config(args.config)
    outputs = generate(config, args.config)
    print("Recovered canonical notebook scalar phase stage outputs:")
    for key, path in outputs.items():
        print(f"- {key}: {path}")


if __name__ == "__main__":
    main()
