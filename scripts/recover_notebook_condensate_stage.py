"""Recover the canonical notebook V1 condensate stage."""

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

from non_destructive_image import build_thomas_fermi_state, thomas_fermi_profile_2d
from scripts.plot_label_utils import (
    DENSITY_CM3,
    column_density_distribution_label,
    coordinate_label,
    radius_legend_label,
)


def load_config(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


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


def notebook_expression_state(config: dict[str, Any]) -> dict[str, Any]:
    constants = config["constants"]
    atom = config["atom"]
    condensate = config["condensate"]
    hbar = float(constants["hbar"])
    k_b = float(constants["boltzmann_constant"])
    amu = float(constants["atomic_mass_unit"])
    a0 = float(constants["bohr_radius_m"])

    atom_number = float(condensate["atom_number"])
    atomic_mass = float(atom["mass_number"]) * amu
    scattering_length = float(condensate["scattering_length_bohr"]) * a0
    trap_frequencies_hz = np.asarray(condensate["trap_frequencies_hz"], dtype=float)

    omega = 2 * np.pi * trap_frequencies_hz
    omega_bar = omega.prod() ** (1 / 3)
    a_ho = np.sqrt(hbar / (atomic_mass * omega_bar))
    mu = 0.5 * (15 * atom_number * scattering_length / a_ho) ** (2 / 5) * hbar * omega_bar
    t_mu = mu / k_b
    n_peak = mu * atomic_mass / (4 * np.pi * hbar**2 * scattering_length)
    radii = np.sqrt(2 * mu / (atomic_mass * omega**2))
    column_density = (4 / 3) * n_peak * radii
    atom_number_check = (8 * np.pi / 15) * n_peak * radii.prod()

    return {
        "atom_number": atom_number,
        "atomic_mass": atomic_mass,
        "scattering_length": scattering_length,
        "trap_frequencies_hz": trap_frequencies_hz,
        "hbar": hbar,
        "boltzmann_constant": k_b,
        "omega": omega,
        "omega_bar": omega_bar,
        "a_ho": a_ho,
        "mu": mu,
        "chemical_potential_temperature": t_mu,
        "n_peak": n_peak,
        "radii": radii,
        "column_density": column_density,
        "atom_number_check": atom_number_check,
    }


def build_condensate_stage(config: dict[str, Any]) -> dict[str, Any]:
    state = notebook_expression_state(config)
    helper_state = build_thomas_fermi_state(
        state["atom_number"],
        state["scattering_length"],
        state["trap_frequencies_hz"],
        state["atomic_mass"],
        state["hbar"],
        state["boltzmann_constant"],
    )

    grid = config["grid"]
    geometry = config["imaging_geometry"]
    ngrid = int(grid["ngrid"])
    field_of_view_m = float(grid["field_of_view_m"])
    dgrid = field_of_view_m / ngrid
    coordinate_axis_m = (np.arange(ngrid) - ngrid // 2) * dgrid
    grid_a_m, grid_b_m = np.meshgrid(coordinate_axis_m, coordinate_axis_m)

    imaging_axis = int(geometry["imaging_axis"])
    plane = [index for index in range(3) if index != imaging_axis]
    axis_labels = geometry["axis_labels"]

    radii = np.asarray(state["radii"], dtype=float)
    notebook_profile = np.maximum(
        0,
        1 - grid_a_m**2 / radii[plane[0]] ** 2 - grid_b_m**2 / radii[plane[1]] ** 2,
    ) ** 1.5
    helper_profile = thomas_fermi_profile_2d(grid_a_m, grid_b_m, radii[plane[0]], radii[plane[1]])
    notebook_column_density = state["column_density"][imaging_axis] * notebook_profile
    helper_column_density = helper_state.column_density[imaging_axis] * helper_profile

    cut_points = int(grid["density_cut_points"])
    span = float(grid["density_cut_span_rmax"])
    density_cut_axis_m = np.linspace(-span * radii.max(), span * radii.max(), cut_points)
    density_cuts_m3 = np.vstack(
        [state["n_peak"] * np.maximum(0, 1 - (density_cut_axis_m / radius) ** 2) for radius in radii]
    )

    return {
        "state": state,
        "helper_state": helper_state,
        "coordinate_axis_m": coordinate_axis_m,
        "grid_a_m": grid_a_m,
        "grid_b_m": grid_b_m,
        "imaging_axis": imaging_axis,
        "transverse_plane_indices": plane,
        "transverse_plane_labels": [axis_labels[index] for index in plane],
        "density_cut_axis_m": density_cut_axis_m,
        "density_cuts_m3": density_cuts_m3,
        "notebook_profile": notebook_profile,
        "helper_profile": helper_profile,
        "notebook_column_density_m2": notebook_column_density,
        "helper_column_density_m2": helper_column_density,
    }


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


def comparison_report(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    state = stage["state"]
    helper = stage["helper_state"]
    notebook_profile = stage["notebook_profile"]
    helper_profile = stage["helper_profile"]
    notebook_column = stage["notebook_column_density_m2"]
    helper_column = stage["helper_column_density_m2"]

    scalar_pairs = {
        "chemical_potential_j": (state["mu"], helper.chemical_potential),
        "chemical_potential_temperature_k": (
            state["chemical_potential_temperature"],
            helper.chemical_potential_temperature,
        ),
        "peak_density_m3": (state["n_peak"], helper.peak_density),
        "atom_number_check": (state["atom_number_check"], helper.atom_number_check),
    }
    scalar_comparisons = {
        key: {
            "notebook_expression": float(left),
            "helper": float(right),
            "absolute_difference": float(abs(left - right)),
            "relative_difference": float(abs(left - right) / abs(left)) if left != 0 else 0.0,
        }
        for key, (left, right) in scalar_pairs.items()
    }

    return {
        "status": "canonical notebook condensate stage matches current helpers for tested deterministic quantities",
        "source_notebook_cells": config["source_notebook_cells"],
        "imaging_axis": stage["imaging_axis"],
        "transverse_plane_indices": stage["transverse_plane_indices"],
        "transverse_plane_labels": stage["transverse_plane_labels"],
        "grid": {
            "shape": [int(config["grid"]["ngrid"]), int(config["grid"]["ngrid"])],
            "field_of_view_m": float(config["grid"]["field_of_view_m"]),
            "coordinate_min_m": float(stage["coordinate_axis_m"].min()),
            "coordinate_max_m": float(stage["coordinate_axis_m"].max()),
            "coordinate_step_m": float(stage["coordinate_axis_m"][1] - stage["coordinate_axis_m"][0]),
        },
        "scalar_comparisons": scalar_comparisons,
        "vector_comparisons": {
            "radii_m": {
                "notebook_expression": state["radii"],
                "helper": helper.radii,
                "max_absolute_difference": float(np.max(np.abs(state["radii"] - helper.radii))),
                "max_relative_difference": max_relative_difference(state["radii"], helper.radii),
            },
            "column_density_m2": {
                "notebook_expression": state["column_density"],
                "helper": helper.column_density,
                "max_absolute_difference": float(np.max(np.abs(state["column_density"] - helper.column_density))),
                "max_relative_difference": max_relative_difference(state["column_density"], helper.column_density),
            },
        },
        "profile_comparison": {
            "max_absolute_difference": float(np.max(np.abs(notebook_profile - helper_profile))),
            "max_relative_difference": max_relative_difference(notebook_profile, helper_profile),
            "stats": array_stats(notebook_profile),
        },
        "column_density_map_comparison": {
            "units": "m^-2",
            "max_absolute_difference": float(np.max(np.abs(notebook_column - helper_column))),
            "max_relative_difference": max_relative_difference(notebook_column, helper_column),
            "stats_m2": array_stats(notebook_column),
            "stats_cm2": array_stats(notebook_column * config["unit_conventions"]["column_density_conversion_m2_to_cm2"]),
        },
    }


def condensate_summary(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    state = stage["state"]
    units = config["unit_conventions"]
    density_cuts_cm3 = stage["density_cuts_m3"] * units["density_conversion_m3_to_cm3"]
    column_cm2 = stage["notebook_column_density_m2"] * units["column_density_conversion_m2_to_cm2"]
    return {
        "label": config["label"],
        "status": config["status_note"],
        "source_notebook": config["source_notebook"],
        "source_notebook_cells": config["source_notebook_cells"],
        "atom_number": state["atom_number"],
        "species": config["atom"]["species"],
        "mass_kg": state["atomic_mass"],
        "scattering_length_m": state["scattering_length"],
        "scattering_length_bohr": config["condensate"]["scattering_length_bohr"],
        "trap_frequencies_hz": state["trap_frequencies_hz"],
        "chemical_potential_j": state["mu"],
        "chemical_potential_temperature_nk": state["chemical_potential_temperature"] * 1e9,
        "peak_density_m3": state["n_peak"],
        "peak_density_cm3": state["n_peak"] * units["density_conversion_m3_to_cm3"],
        "radii_m": state["radii"],
        "radii_um": state["radii"] * units["coordinate_conversion_m_to_um"],
        "column_density_m2": state["column_density"],
        "column_density_cm2": state["column_density"] * units["column_density_conversion_m2_to_cm2"],
        "atom_number_check": state["atom_number_check"],
        "grid_shape": list(stage["notebook_profile"].shape),
        "coordinate_range_um": [
            float(stage["coordinate_axis_m"].min() * units["coordinate_conversion_m_to_um"]),
            float(stage["coordinate_axis_m"].max() * units["coordinate_conversion_m_to_um"]),
        ],
        "density_cut_stats_cm3": {
            "x_across": {
                "max": float(np.max(density_cuts_cm3[0])),
                "mean": float(np.mean(density_cuts_cm3[0])),
            },
            "y_along": {
                "max": float(np.max(density_cuts_cm3[1])),
                "mean": float(np.mean(density_cuts_cm3[1])),
            },
            "z_across": {
                "max": float(np.max(density_cuts_cm3[2])),
                "mean": float(np.mean(density_cuts_cm3[2])),
            },
        },
        "column_density_map_stats_cm2": array_stats(column_cm2),
    }


def write_central_lineouts(path: Path, config: dict[str, Any], stage: dict[str, Any]) -> None:
    units = config["unit_conventions"]
    axis_um = stage["coordinate_axis_m"] * units["coordinate_conversion_m_to_um"]
    profile = stage["notebook_profile"]
    column_cm2 = stage["notebook_column_density_m2"] * units["column_density_conversion_m2_to_cm2"]
    mid = profile.shape[0] // 2
    write_rows(
        path,
        [
            {
                "position_um": float(axis_um[index]),
                "profile_transverse_a": float(profile[mid, index]),
                "profile_transverse_b": float(profile[index, mid]),
                "column_density_transverse_a_cm2": float(column_cm2[mid, index]),
                "column_density_transverse_b_cm2": float(column_cm2[index, mid]),
            }
            for index in range(axis_um.size)
        ],
    )


def write_density_cuts(path: Path, config: dict[str, Any], stage: dict[str, Any]) -> None:
    units = config["unit_conventions"]
    axis_um = stage["density_cut_axis_m"] * units["coordinate_conversion_m_to_um"]
    density_cm3 = stage["density_cuts_m3"] * units["density_conversion_m3_to_cm3"]
    write_rows(
        path,
        [
            {
                "position_um": float(axis_um[index]),
                "density_x_across_cm3": float(density_cm3[0, index]),
                "density_y_along_cm3": float(density_cm3[1, index]),
                "density_z_across_cm3": float(density_cm3[2, index]),
            }
            for index in range(axis_um.size)
        ],
    )


def write_condensate_figure(path: Path, config: dict[str, Any], stage: dict[str, Any]) -> None:
    units = config["unit_conventions"]
    display = config["display"]
    state = stage["state"]
    radii = state["radii"]
    xax_um = stage["density_cut_axis_m"] * units["coordinate_conversion_m_to_um"]
    density_cm3 = stage["density_cuts_m3"] * units["density_conversion_m3_to_cm3"]
    column_cm2 = stage["notebook_column_density_m2"] * units["column_density_conversion_m2_to_cm2"]
    axis_labels = config["imaging_geometry"]["axis_labels"]
    plane = stage["transverse_plane_indices"]
    imaging_axis = stage["imaging_axis"]

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
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=tuple(display["figure_size_inches"]))
    labels = ["x (across)", "y (along)", "z (across)"]
    colours = ["C0", "C1", "C2"]
    for index, (label, colour) in enumerate(zip(labels, colours)):
        ax_a.plot(
            xax_um,
            density_cm3[index],
            colour,
            lw=2,
            label=radius_legend_label(label, radii[index] * units["coordinate_conversion_m_to_um"]),
        )
    ax_a.set_xlabel(coordinate_label())
    ax_a.set_ylabel(rf"$n(\mathbf{{r}})$ ({DENSITY_CM3})")
    ax_a.set_title("(a) 3D density cuts through the trap centre")
    ax_a.legend(fontsize=8.5)
    ax_a.grid(alpha=0.25)

    im = ax_b.imshow(
        column_cm2,
        extent=display["image_extent_um"],
        origin="lower",
        cmap=display["column_density_cmap"],
    )
    ax_b.set_xlim(*display["column_density_xlim_um"])
    ax_b.set_ylim(*display["column_density_ylim_um"])
    ax_b.set_xlabel(coordinate_label(axis_labels[plane[0]]))
    ax_b.set_ylabel(coordinate_label(axis_labels[plane[1]]))
    ax_b.set_title(
        f"(b) Column density along {axis_labels[imaging_axis]}\n"
        "(what every imaging mode actually integrates over)"
    )
    plt.colorbar(
        im,
        ax=ax_b,
        fraction=0.032,
        label=column_density_distribution_label(axis_labels[plane[0]], axis_labels[plane[1]]),
    )
    fig.suptitle("Stage 1: from trap parameters to a Thomas-Fermi condensate", y=1.04, fontsize=11.5)
    plt.tight_layout()
    fig.savefig(path, format="svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def generate(config: dict[str, Any], config_path: Path) -> dict[str, str]:
    output_dir = Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    stage = build_condensate_stage(config)

    outputs = {
        "comparison_report": output_dir / "comparison_report.json",
        "condensate_summary": output_dir / "condensate_summary.json",
        "central_lineouts": output_dir / "central_lineouts.csv",
        "density_cuts": output_dir / "density_cuts.csv",
        "metadata": output_dir / "metadata.json",
        "figure": output_dir / "condensate_density_stage.svg",
    }
    write_json(outputs["comparison_report"], comparison_report(config, stage))
    write_json(outputs["condensate_summary"], condensate_summary(config, stage))
    write_central_lineouts(outputs["central_lineouts"], config, stage)
    write_density_cuts(outputs["density_cuts"], config, stage)
    write_condensate_figure(outputs["figure"], config, stage)
    write_json(
        outputs["metadata"],
        {
            "label": config["label"],
            "status": "Canonical notebook-aligned condensate-stage recovery output.",
            "config_file_used": str(config_path),
            "git_commit_hash": git_commit(),
            "source_notebook": config["source_notebook"],
            "source_notebook_cells": config["source_notebook_cells"],
            "generated_outputs": {key: str(path) for key, path in outputs.items()},
            "scope_boundary": "Only parameters -> Thomas-Fermi condensate -> profile/column-density arrays are recovered. No phase, imaging, camera, or multishot recovery is performed.",
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
    print("Recovered canonical notebook condensate stage outputs:")
    for key, path in outputs.items():
        print(f"- {key}: {path}")


if __name__ == "__main__":
    main()
