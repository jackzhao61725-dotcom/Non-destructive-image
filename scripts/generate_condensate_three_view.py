"""Generate a three-view condensate column-density projection figure.

This is a notebook-aligned condensate-model extension built from the recovered
canonical condensate stage. It is not an exact recovered notebook figure unless
the notebook later proves to contain the same three-view layout.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from non_destructive_image import thomas_fermi_profile_2d
from scripts.plot_label_utils import (
    column_density_distribution_label,
    coordinate_label,
    peak_column_density_symbol,
)
from scripts.recover_notebook_condensate_stage import (
    array_stats,
    build_condensate_stage,
    git_commit,
    load_config,
    write_json,
    write_rows,
)


OUTPUT_DIR = Path("results/notebook_aligned_recovery/condensate_three_view")


def _projection_view(
    config: dict[str, Any],
    condensate_stage: dict[str, Any],
    integrated_axis: int,
) -> dict[str, Any]:
    axis_labels = config["imaging_geometry"]["axis_labels"]
    plane = [index for index in range(3) if index != integrated_axis]
    coordinate_axis_m = condensate_stage["coordinate_axis_m"]
    grid_a_m, grid_b_m = np.meshgrid(coordinate_axis_m, coordinate_axis_m)
    state = condensate_stage["state"]
    radii = np.asarray(state["radii"], dtype=float)

    profile = thomas_fermi_profile_2d(grid_a_m, grid_b_m, radii[plane[0]], radii[plane[1]])
    column_density_m2 = state["column_density"][integrated_axis] * profile

    return {
        "integrated_axis_index": integrated_axis,
        "integrated_axis_label": axis_labels[integrated_axis],
        "display_plane_indices": plane,
        "display_plane_labels": [axis_labels[index] for index in plane],
        "profile": profile,
        "column_density_m2": column_density_m2,
        "peak_column_density_m2": float(state["column_density"][integrated_axis]),
        "absolute_column_density": True,
        "normalised": False,
    }


def build_three_view_stage(config: dict[str, Any]) -> dict[str, Any]:
    """Return all three TF column-density projections from the canonical stage."""

    condensate_stage = build_condensate_stage(config)
    views = [_projection_view(config, condensate_stage, axis) for axis in range(3)]
    return {
        "condensate_stage": condensate_stage,
        "views": views,
        "axis_labels": config["imaging_geometry"]["axis_labels"],
        "coordinate_axis_m": condensate_stage["coordinate_axis_m"],
    }


def three_view_summary(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    units = config["unit_conventions"]
    return {
        "label": config["label"],
        "status": "Notebook-aligned condensate-model three-view projection extension.",
        "source_stage": "scripts/recover_notebook_condensate_stage.py",
        "source_notebook": config["source_notebook"],
        "not_exact_notebook_figure": True,
        "no_experimental_calibration_applied": True,
        "not_final_calibrated_prediction": True,
        "quantity": "absolute Thomas-Fermi column density",
        "internal_unit": "m^-2",
        "normalised": False,
        "grid_shape": [int(config["grid"]["ngrid"]), int(config["grid"]["ngrid"])],
        "field_of_view_m": float(config["grid"]["field_of_view_m"]),
        "display_extent_um": config["display"]["image_extent_um"],
        "coordinate_range_um": [
            float(stage["coordinate_axis_m"].min() * units["coordinate_conversion_m_to_um"]),
            float(stage["coordinate_axis_m"].max() * units["coordinate_conversion_m_to_um"]),
        ],
        "views": [
            {
                "integrated_axis": view["integrated_axis_label"],
                "display_plane": view["display_plane_labels"],
                "axis_convention": (
                    f"integrate along {view['integrated_axis_label']} -> "
                    f"display {view['display_plane_labels'][0]}-{view['display_plane_labels'][1]} plane"
                ),
                "distribution_label": (
                    f"n_col({view['display_plane_labels'][0]},{view['display_plane_labels'][1]})"
                ),
                "peak_scalar_symbol": f"n_tilde_{view['integrated_axis_label']}",
                "peak_scalar_mathtext": peak_column_density_symbol(view["integrated_axis_label"]),
                "notation_note": (
                    "The plotted 2D distribution is n_col over the displayed plane. "
                    "The peak scalar value corresponds to the thesis-table tilde-n symbol "
                    f"for integration along {view['integrated_axis_label']}."
                ),
                "absolute_column_density": True,
                "normalised": False,
                "peak_column_density_m2": view["peak_column_density_m2"],
                "map_stats_m2": array_stats(view["column_density_m2"]),
            }
            for view in stage["views"]
        ],
    }


def write_central_lineouts(path: Path, config: dict[str, Any], stage: dict[str, Any]) -> None:
    units = config["unit_conventions"]
    position_um = stage["coordinate_axis_m"] * units["coordinate_conversion_m_to_um"]
    mid = position_um.size // 2
    rows = []
    for index, position in enumerate(position_um):
        row: dict[str, float] = {"position_um": float(position)}
        for view in stage["views"]:
            axis = view["integrated_axis_label"]
            plane_a, plane_b = view["display_plane_labels"]
            image = view["column_density_m2"]
            row[f"{axis}_integrated_{plane_a}_line_m2"] = float(image[mid, index])
            row[f"{axis}_integrated_{plane_b}_line_m2"] = float(image[index, mid])
        rows.append(row)
    write_rows(path, rows)


def write_three_view_figure(path: Path, config: dict[str, Any], stage: dict[str, Any]) -> None:
    display = config["display"]
    figure_size = (14.5, 4.3)
    plt.rcParams.update(
        {
            "figure.figsize": figure_size,
            "font.size": 10.5,
            "axes.labelsize": 11,
            "axes.titlesize": 11,
            "figure.dpi": 110,
        }
    )

    fig, axes = plt.subplots(1, 3, figsize=figure_size, constrained_layout=True)
    for axis, view in zip(axes, stage["views"]):
        plane_a, plane_b = view["display_plane_labels"]
        image = view["column_density_m2"]
        im = axis.imshow(
            image,
            extent=display["image_extent_um"],
            origin="lower",
            cmap=display["column_density_cmap"],
        )
        axis.set_xlabel(coordinate_label(plane_a))
        axis.set_ylabel(coordinate_label(plane_b))
        axis.set_xlim(*display["image_extent_um"][:2])
        axis.set_ylim(*display["image_extent_um"][2:])
        axis.set_title(
            rf"integrate along ${view['integrated_axis_label']}$"
            "\n"
            rf"display ${plane_a}$-${plane_b}$ plane"
        )
        plt.colorbar(
            im,
            ax=axis,
            fraction=0.046,
            pad=0.035,
            label=column_density_distribution_label(plane_a, plane_b, display_unit="m^-2"),
        )

    fig.suptitle(
        "Notebook-aligned condensate-model extension: three absolute column-density projections",
        fontsize=12,
    )
    fig.savefig(path, format="svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def generate(config: dict[str, Any], config_path: Path) -> dict[str, str]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    stage = build_three_view_stage(config)
    outputs = {
        "figure": OUTPUT_DIR / "condensate_three_view.svg",
        "summary": OUTPUT_DIR / "condensate_three_view_summary.json",
        "central_lineouts": OUTPUT_DIR / "central_lineouts.csv",
        "metadata": OUTPUT_DIR / "metadata.json",
    }

    write_three_view_figure(outputs["figure"], config, stage)
    write_json(outputs["summary"], three_view_summary(config, stage))
    write_central_lineouts(outputs["central_lineouts"], config, stage)
    write_json(
        outputs["metadata"],
        {
            "label": config["label"],
            "status": "Notebook-aligned condensate-model extension.",
            "config_file_used": str(config_path),
            "git_commit_hash": git_commit(),
            "source_stage": "scripts/recover_notebook_condensate_stage.py",
            "not_exact_notebook_figure": True,
            "no_experimental_calibration_applied": True,
            "not_final_calibrated_prediction": True,
            "generated_outputs": {key: str(value) for key, value in outputs.items()},
            "axis_conventions": [
                {
                    "integrated_axis": view["integrated_axis_label"],
                    "display_plane": view["display_plane_labels"],
                    "distribution_label": (
                        f"n_col({view['display_plane_labels'][0]},{view['display_plane_labels'][1]})"
                    ),
                    "peak_scalar_symbol": f"n_tilde_{view['integrated_axis_label']}",
                    "absolute_column_density": True,
                    "normalised": False,
                }
                for view in stage["views"]
            ],
        },
    )
    return {key: str(value) for key, value in outputs.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config = load_config(args.config)
    outputs = generate(config, args.config)
    print("Generated notebook-aligned condensate three-view extension outputs:")
    for key, value in outputs.items():
        print(f"- {key}: {value}")


if __name__ == "__main__":
    main()
