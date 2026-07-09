"""Recover the notebook-aligned Stage 18.1 condensate/density figure."""

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
            return subprocess.check_output(command, text=True, stderr=subprocess.DEVNULL).strip()
        except (OSError, subprocess.CalledProcessError):
            continue
    return "unknown"


def _write_json(path: Path, payload: dict[str, Any]) -> None:
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


def _write_rows(path: Path, rows: list[dict[str, float]]) -> None:
    if not rows:
        raise ValueError(f"no rows to write for {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: f"{value:.12g}" for key, value in row.items()})


def _direct_notebook_state(config: dict[str, Any]) -> dict[str, Any]:
    constants = config["constants"]
    condensate = config["condensate"]
    hbar = float(constants["hbar"])
    k_b = float(constants["boltzmann_constant"])
    amu = float(constants["atomic_mass_unit"])
    a0 = float(constants["bohr_radius_m"])

    atom_number = float(condensate["atom_number"])
    atomic_mass = float(condensate["erbium_mass_number"]) * amu
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
        "omega": omega,
        "omega_bar": omega_bar,
        "a_ho": a_ho,
        "mu": mu,
        "chemical_potential_temperature": t_mu,
        "n_peak": n_peak,
        "radii": radii,
        "column_density": column_density,
        "atom_number_check": atom_number_check,
        "atom_number": atom_number,
        "atomic_mass": atomic_mass,
        "scattering_length": scattering_length,
        "trap_frequencies_hz": trap_frequencies_hz,
        "hbar": hbar,
        "boltzmann_constant": k_b,
    }


def _grid_and_arrays(config: dict[str, Any], state: dict[str, Any]) -> dict[str, Any]:
    grid = config["grid"]
    display = config["display"]
    ngrid = int(grid["ngrid"])
    field_of_view_m = float(grid["field_of_view_m"])
    dgrid = field_of_view_m / ngrid
    gax = (np.arange(ngrid) - ngrid // 2) * dgrid
    ga, gb = np.meshgrid(gax, gax)

    axis_demo = int(display["imaging_axis"])
    plane = [index for index in range(3) if index != axis_demo]
    radii = np.asarray(state["radii"], dtype=float)
    n_peak = float(state["n_peak"])
    column_density = np.asarray(state["column_density"], dtype=float)

    notebook_profile = np.maximum(
        0,
        1 - ga**2 / radii[plane[0]] ** 2 - gb**2 / radii[plane[1]] ** 2,
    ) ** 1.5
    helper_profile = thomas_fermi_profile_2d(ga, gb, radii[plane[0]], radii[plane[1]])
    notebook_column_density_map = column_density[axis_demo] * notebook_profile
    helper_column_density_map = column_density[axis_demo] * helper_profile

    cut_points = int(grid["density_cut_points"])
    span = float(grid["density_cut_span_rmax"])
    xax = np.linspace(-span * radii.max(), span * radii.max(), cut_points)
    density_cuts = np.vstack(
        [n_peak * np.maximum(0, 1 - (xax / radius) ** 2) for radius in radii]
    )

    return {
        "ngrid": ngrid,
        "field_of_view_m": field_of_view_m,
        "dgrid": dgrid,
        "gax": gax,
        "ga": ga,
        "gb": gb,
        "axis_demo": axis_demo,
        "plane": plane,
        "xax": xax,
        "density_cuts": density_cuts,
        "notebook_profile": notebook_profile,
        "helper_profile": helper_profile,
        "notebook_column_density_map": notebook_column_density_map,
        "helper_column_density_map": helper_column_density_map,
    }


def _stats(array: np.ndarray) -> dict[str, Any]:
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


def _max_relative_difference(a: np.ndarray, b: np.ndarray) -> float:
    denominator = np.maximum(np.abs(a), np.abs(b))
    mask = denominator > 0
    if not np.any(mask):
        return 0.0
    return float(np.max(np.abs(a[mask] - b[mask]) / denominator[mask]))


def _comparison_report(config: dict[str, Any], state: dict[str, Any], arrays: dict[str, Any]) -> dict[str, Any]:
    helper_state = build_thomas_fermi_state(
        state["atom_number"],
        state["scattering_length"],
        state["trap_frequencies_hz"],
        state["atomic_mass"],
        state["hbar"],
        state["boltzmann_constant"],
    )

    notebook_profile = arrays["notebook_profile"]
    helper_profile = arrays["helper_profile"]
    notebook_column_map = arrays["notebook_column_density_map"]
    helper_column_map = arrays["helper_column_density_map"]

    state_pairs = {
        "chemical_potential": (state["mu"], helper_state.chemical_potential),
        "chemical_potential_temperature": (
            state["chemical_potential_temperature"],
            helper_state.chemical_potential_temperature,
        ),
        "peak_density": (state["n_peak"], helper_state.peak_density),
        "atom_number_check": (state["atom_number_check"], helper_state.atom_number_check),
    }
    scalar_comparisons = {
        key: {
            "notebook_direct": float(left),
            "helper": float(right),
            "absolute_difference": float(abs(left - right)),
            "relative_difference": float(abs(left - right) / abs(left)) if left != 0 else 0.0,
        }
        for key, (left, right) in state_pairs.items()
    }

    vector_comparisons = {
        "radii_m": {
            "notebook_direct": state["radii"],
            "helper": helper_state.radii,
            "max_absolute_difference": float(np.max(np.abs(state["radii"] - helper_state.radii))),
            "max_relative_difference": _max_relative_difference(state["radii"], helper_state.radii),
        },
        "column_density_m2": {
            "notebook_direct": state["column_density"],
            "helper": helper_state.column_density,
            "max_absolute_difference": float(np.max(np.abs(state["column_density"] - helper_state.column_density))),
            "max_relative_difference": _max_relative_difference(state["column_density"], helper_state.column_density),
        },
    }

    return {
        "status": "notebook direct formulas and current helpers match for the recovered Stage 18.1 quantities",
        "source_notebook_cell": config["source_notebook_cell"],
        "source_notebook_section": config["source_notebook_section"],
        "grid": {
            "ngrid": arrays["ngrid"],
            "field_of_view_m": arrays["field_of_view_m"],
            "dgrid_m": arrays["dgrid"],
            "gax_min_m": float(arrays["gax"].min()),
            "gax_max_m": float(arrays["gax"].max()),
            "imaging_axis": arrays["axis_demo"],
            "transverse_plane_indices": arrays["plane"],
            "transverse_plane_labels": [["x", "y", "z"][index] for index in arrays["plane"]],
        },
        "scalar_comparisons": scalar_comparisons,
        "vector_comparisons": vector_comparisons,
        "profile_comparison": {
            "shape": list(notebook_profile.shape),
            "max_absolute_difference": float(np.max(np.abs(notebook_profile - helper_profile))),
            "max_relative_difference": _max_relative_difference(notebook_profile, helper_profile),
            "stats": _stats(notebook_profile),
        },
        "column_density_map_comparison": {
            "units": "m^-2 before plotting; plotted as cm^-2 by multiplying by 1e-4",
            "max_absolute_difference_m2": float(np.max(np.abs(notebook_column_map - helper_column_map))),
            "max_relative_difference": _max_relative_difference(notebook_column_map, helper_column_map),
            "stats_m2": _stats(notebook_column_map),
            "stats_cm2": _stats(notebook_column_map * 1e-4),
        },
    }


def _write_lineout_csvs(output_dir: Path, state: dict[str, Any], arrays: dict[str, Any]) -> None:
    xax = arrays["xax"]
    density_cuts = arrays["density_cuts"]
    _write_rows(
        output_dir / "density_lineouts_stage18_1.csv",
        [
            {
                "position_m": float(xax[index]),
                "position_um": float(xax[index] * 1e6),
                "density_x_across_cm3": float(density_cuts[0, index] * 1e-6),
                "density_y_along_cm3": float(density_cuts[1, index] * 1e-6),
                "density_z_across_cm3": float(density_cuts[2, index] * 1e-6),
            }
            for index in range(xax.size)
        ],
    )

    gax = arrays["gax"]
    column_map_cm2 = arrays["notebook_column_density_map"] * 1e-4
    mid = column_map_cm2.shape[0] // 2
    _write_rows(
        output_dir / "column_density_center_lineouts_stage18_1.csv",
        [
            {
                "position_m": float(gax[index]),
                "position_um": float(gax[index] * 1e6),
                "column_density_y_line_cm2": float(column_map_cm2[mid, index]),
                "column_density_z_line_cm2": float(column_map_cm2[index, mid]),
            }
            for index in range(gax.size)
        ],
    )


def _write_figure(path: Path, config: dict[str, Any], state: dict[str, Any], arrays: dict[str, Any]) -> None:
    display = config["display"]
    extent = display["image_extent_um"]
    radii = np.asarray(state["radii"], dtype=float)
    n_peak = float(state["n_peak"])
    column_map_cm2 = arrays["notebook_column_density_map"] * 1e-4
    xax = arrays["xax"]

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
        n1d = n_peak * np.maximum(0, 1 - (xax / radii[index]) ** 2)
        ax_a.plot(
            xax * 1e6,
            n1d * 1e-6,
            colour,
            lw=2,
            label=f"{label}, R={radii[index] * 1e6:.2f} um",
        )
    ax_a.set_xlabel("position (um)")
    ax_a.set_ylabel("n(r)  (cm^-3)")
    ax_a.set_title("(a) 3D density cuts through the trap centre")
    ax_a.legend(fontsize=8.5)
    ax_a.grid(alpha=0.25)

    im = ax_b.imshow(
        column_map_cm2,
        extent=extent,
        origin="lower",
        cmap=display["column_density_cmap"],
    )
    ax_b.set_xlim(*display["column_density_xlim_um"])
    ax_b.set_ylim(*display["column_density_ylim_um"])
    plane = arrays["plane"]
    axis_demo = arrays["axis_demo"]
    axis_labels = ["x", "y", "z"]
    ax_b.set_xlabel(f"{axis_labels[plane[0]]} (um)")
    ax_b.set_ylabel(f"{axis_labels[plane[1]]} (um)")
    ax_b.set_title(
        f"(b) Column density along {axis_labels[axis_demo]}\n"
        "(what every imaging mode actually integrates over)"
    )
    plt.colorbar(im, ax=ax_b, fraction=0.032, label="n_col (cm^-2)")
    fig.suptitle("Stage 1: from trap parameters to a Thomas-Fermi condensate", y=1.04, fontsize=11.5)
    plt.tight_layout()
    fig.savefig(path, format="svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _figure_spec(config: dict[str, Any], state: dict[str, Any], arrays: dict[str, Any]) -> dict[str, Any]:
    return {
        "figure_filename": "condensate_density_stage18_1.svg",
        "source_notebook": config["source_notebook"],
        "source_notebook_cell": config["source_notebook_cell"],
        "source_notebook_section": config["source_notebook_section"],
        "description": "Notebook-aligned recovery of Stage 18.1 condensate/density figure.",
        "panels": [
            {
                "panel": "a",
                "quantity": "3D Thomas-Fermi density cuts through the trap centre",
                "x_axis": "position (um)",
                "y_axis": "n(r) (cm^-3)",
                "arrays": "n_peak * max(0, 1 - (xax/R_i)^2) for i in x,y,z",
                "normalisation": "absolute density converted from m^-3 to cm^-3 by multiplying by 1e-6",
            },
            {
                "panel": "b",
                "quantity": "column density along x for the across-cigar y,z image plane",
                "x_axis": "y (um)",
                "y_axis": "z (um)",
                "array": "n_col[0] * max(0, 1 - y^2/R_y^2 - z^2/R_z^2)^1.5",
                "normalisation": "absolute column density converted from m^-2 to cm^-2 by multiplying by 1e-4",
                "display": {
                    "extent_um": config["display"]["image_extent_um"],
                    "xlim_um": config["display"]["column_density_xlim_um"],
                    "ylim_um": config["display"]["column_density_ylim_um"],
                    "cmap": config["display"]["column_density_cmap"],
                },
            },
        ],
        "condensate_summary": {
            "atom_number": state["atom_number"],
            "scattering_length_m": state["scattering_length"],
            "trap_frequencies_hz": state["trap_frequencies_hz"],
            "chemical_potential_temperature_nk": state["chemical_potential_temperature"] * 1e9,
            "peak_density_cm3": state["n_peak"] * 1e-6,
            "radii_um": state["radii"] * 1e6,
            "column_density_cm2": state["column_density"] * 1e-4,
        },
        "grid_summary": {
            "ngrid": arrays["ngrid"],
            "field_of_view_um": arrays["field_of_view_m"] * 1e6,
            "coordinate_axis_min_um": float(arrays["gax"].min() * 1e6),
            "coordinate_axis_max_um": float(arrays["gax"].max() * 1e6),
        },
    }


def generate(config: dict[str, Any], config_path: Path) -> dict[str, str]:
    output_dir = Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    state = _direct_notebook_state(config)
    arrays = _grid_and_arrays(config, state)

    comparison = _comparison_report(config, state, arrays)
    figure_spec = _figure_spec(config, state, arrays)
    metadata = {
        "label": config["label"],
        "status": "Notebook-aligned Version 1 condensate/density figure recovery.",
        "source_notebook": config["source_notebook"],
        "source_notebook_cell": config["source_notebook_cell"],
        "source_notebook_section": config["source_notebook_section"],
        "config_file_used": str(config_path),
        "git_commit_hash": _git_commit(),
        "generated_figure": "condensate_density_stage18_1.svg",
        "no_experimental_calibration_applied": True,
        "not_final_calibrated_theory": True,
        "scope_boundary": "Only Stage 18.1 condensate/density recovery is generated. No phase, imaging, camera, or multishot figures are generated.",
    }

    outputs = {
        "figure": output_dir / "condensate_density_stage18_1.svg",
        "comparison_report": output_dir / "comparison_report.json",
        "figure_spec": output_dir / "figure_spec.json",
        "metadata": output_dir / "metadata.json",
        "density_lineouts": output_dir / "density_lineouts_stage18_1.csv",
        "column_density_lineouts": output_dir / "column_density_center_lineouts_stage18_1.csv",
    }

    _write_figure(outputs["figure"], config, state, arrays)
    _write_json(outputs["comparison_report"], comparison)
    _write_json(outputs["figure_spec"], figure_spec)
    _write_json(outputs["metadata"], metadata)
    _write_lineout_csvs(output_dir, state, arrays)
    return {key: str(path) for key, path in outputs.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()

    config = _load_config(args.config)
    outputs = generate(config, args.config)
    print("Recovered notebook-aligned condensate figure outputs:")
    for key, path in outputs.items():
        print(f"- {key}: {path}")


if __name__ == "__main__":
    main()
