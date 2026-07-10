"""Recover the canonical notebook V1 DGI image stage."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from non_destructive_image import simulate_dgi_image
from scripts.recover_notebook_condensate_stage import load_config
from scripts.recover_notebook_pci_stage import (
    build_pupil,
    complex_array_stats,
    field_comparison,
    git_commit,
    real_array_stats,
    write_json,
    write_rows,
)
from scripts.recover_notebook_phase_stage import build_phase_stage
from scripts.plot_label_utils import NORMALISED_INTENSITY, coordinate_label, cut_label


def build_dgi_stage(config: dict[str, Any]) -> dict[str, Any]:
    phase_stage = build_phase_stage(config)
    dgi_config = config["dgi_recovery"]
    pupil_stage = build_pupil(config)
    phase_map = phase_stage["notebook_phase_map_rad"]
    stop_optical_depth = float(dgi_config["stop_optical_depth"])

    notebook_object_field = np.exp(1j * phase_map)
    notebook_scattered_field = notebook_object_field - 1
    notebook_propagated_scattered_field = np.fft.ifft2(
        np.fft.fft2(notebook_scattered_field) * pupil_stage["pupil"]
    )
    notebook_reference_field = 10 ** (-stop_optical_depth / 2)
    notebook_dgi_image_intensity = np.abs(notebook_reference_field + notebook_propagated_scattered_field) ** 2

    helper_result = simulate_dgi_image(
        phase_map,
        pupil_stage["pupil"],
        stop_optical_depth=stop_optical_depth,
        return_intermediates=True,
    )

    return {
        "phase_stage": phase_stage,
        "pupil_stage": pupil_stage,
        "stop_optical_depth": stop_optical_depth,
        "dgi_reference_intensity": notebook_reference_field**2,
        "notebook_object_field": notebook_object_field,
        "notebook_scattered_field": notebook_scattered_field,
        "notebook_pupil": pupil_stage["pupil"],
        "notebook_reference_field": notebook_reference_field,
        "notebook_propagated_scattered_field": notebook_propagated_scattered_field,
        "notebook_dgi_image_intensity": notebook_dgi_image_intensity,
        "helper_object_field": helper_result["object_field"],
        "helper_scattered_field": helper_result["scattered_field"],
        "helper_reference_field": helper_result["dgi_reference_field"],
        "helper_propagated_scattered_field": helper_result["propagated_scattered_field"],
        "helper_dgi_image_intensity": helper_result["dgi_image_intensity"],
    }


def comparison_report(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    intensity = stage["notebook_dgi_image_intensity"]
    return {
        "status": "canonical notebook DGI stage matches simulate_dgi_image for tested deterministic quantities",
        "source_notebook_cells": {
            "dgi_transfer_curve": 16,
            "fourier_optics_sim_image": 18,
        },
        "notebook_formula": {
            "object_field": "np.exp(1j * phase_map)",
            "scattered_field": "object_field - 1",
            "propagation": "np.fft.ifft2(np.fft.fft2(scattered_field) * pupil)",
            "reference_field": "10 ** (-OD / 2)",
            "dgi_image_intensity": "np.abs(reference_field + propagated_scattered_field) ** 2",
        },
        "pupil": {
            "shape": list(stage["notebook_pupil"].shape),
            "nonzero_pixels": int(np.count_nonzero(stage["notebook_pupil"])),
            "numerical_aperture": stage["pupil_stage"]["numerical_aperture"],
            "dgrid_m": stage["pupil_stage"]["dgrid_m"],
        },
        "dgi_stop": {
            "stop_optical_depth": stage["stop_optical_depth"],
            "reference_field_notebook": stage["notebook_reference_field"],
            "reference_field_helper": stage["helper_reference_field"],
            "reference_field_absolute_difference": abs(
                stage["notebook_reference_field"] - stage["helper_reference_field"]
            ),
            "reference_intensity": stage["dgi_reference_intensity"],
        },
        "object_field": {
            "stats": complex_array_stats(stage["notebook_object_field"]),
            **field_comparison(stage["notebook_object_field"], stage["helper_object_field"]),
        },
        "propagated_scattered_field": {
            "stats": complex_array_stats(stage["notebook_propagated_scattered_field"]),
            **field_comparison(
                stage["notebook_propagated_scattered_field"],
                stage["helper_propagated_scattered_field"],
            ),
        },
        "dgi_image_intensity": {
            "units": "incident-I0-normalised intensity",
            "stats": real_array_stats(intensity),
            **field_comparison(intensity, stage["helper_dgi_image_intensity"]),
        },
        "scope_boundary": "Only scalar phase map -> DGI image is recovered. Faraday, camera, and multishot stages are not generated.",
    }


def dgi_summary(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    return {
        "label": config["label"],
        "status": "Canonical notebook-aligned DGI-stage recovery output.",
        "source_notebook": config["source_notebook"],
        "source_notebook_cells": {
            "dgi_transfer_curve": 16,
            "fourier_optics_sim_image": 18,
        },
        "detuning_hz": stage["phase_stage"]["detuning_hz"],
        "phase_peak_rad": stage["phase_stage"]["notebook_phase_peak_rad"],
        "stop_optical_depth": stage["stop_optical_depth"],
        "dgi_reference_field": stage["notebook_reference_field"],
        "dgi_reference_intensity": stage["dgi_reference_intensity"],
        "numerical_aperture": stage["pupil_stage"]["numerical_aperture"],
        "pupil_nonzero_pixels": int(np.count_nonzero(stage["notebook_pupil"])),
        "dgi_image_stats": real_array_stats(stage["notebook_dgi_image_intensity"]),
        "scope_boundary": "Only phase_map -> DGI image is recovered. No Faraday, camera, or multishot figures are generated.",
        "no_experimental_calibration_applied": True,
    }


def write_dgi_lineouts(path: Path, config: dict[str, Any], stage: dict[str, Any]) -> None:
    units = config["unit_conventions"]
    coordinate_um = (
        stage["phase_stage"]["condensate_stage"]["coordinate_axis_m"]
        * units["coordinate_conversion_m_to_um"]
    )
    intensity = stage["notebook_dgi_image_intensity"]
    mid = intensity.shape[0] // 2
    write_rows(
        path,
        [
            {
                "position_um": float(coordinate_um[index]),
                "dgi_intensity_y_line_incident_i0": float(intensity[mid, index]),
                "dgi_intensity_z_line_incident_i0": float(intensity[index, mid]),
                "dgi_signal_y_line_minus_reference_intensity": float(
                    intensity[mid, index] - stage["dgi_reference_intensity"]
                ),
                "dgi_signal_z_line_minus_reference_intensity": float(
                    intensity[index, mid] - stage["dgi_reference_intensity"]
                ),
            }
            for index in range(coordinate_um.size)
        ],
    )


def write_dgi_figure(path: Path, config: dict[str, Any], stage: dict[str, Any]) -> None:
    display = config["display"]
    dgi_config = config["dgi_recovery"]
    units = config["unit_conventions"]
    coordinate_um = (
        stage["phase_stage"]["condensate_stage"]["coordinate_axis_m"]
        * units["coordinate_conversion_m_to_um"]
    )
    intensity = stage["notebook_dgi_image_intensity"]
    mid = intensity.shape[0] // 2

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
        figsize=tuple(dgi_config["figure_size_inches"]),
        gridspec_kw={"width_ratios": [1.5, 1]},
    )
    im = ax_map.imshow(
        intensity,
        extent=display["image_extent_um"],
        origin="lower",
        cmap=dgi_config["map_cmap"],
    )
    ax_map.set_xlim(*display["column_density_xlim_um"])
    ax_map.set_ylim(*display["column_density_ylim_um"])
    ax_map.set_xlabel(coordinate_label("y"))
    ax_map.set_ylabel(coordinate_label("z"))
    ax_map.set_title("DGI image intensity, incident-$I_0$ convention")
    plt.colorbar(im, ax=ax_map, fraction=0.03, label=NORMALISED_INTENSITY)

    ax_line.plot(coordinate_um, intensity[mid, :], "C1", lw=2, label=cut_label("y"))
    ax_line.plot(coordinate_um, intensity[:, mid], "C0", lw=2, label=cut_label("z"))
    ax_line.axhline(stage["dgi_reference_intensity"], color="gray", ls=":", lw=1)
    ax_line.annotate(
        f"$10^{{-OD}}$ = {stage['dgi_reference_intensity']:.4g}",
        (12, max(stage["dgi_reference_intensity"], float(np.max(intensity)) * 0.08)),
        fontsize=9,
        color="gray",
    )
    ax_line.set_xlim(-45, 45)
    ax_line.set_xlabel(coordinate_label())
    ax_line.set_ylabel(NORMALISED_INTENSITY)
    ax_line.set_title("centre cuts")
    ax_line.legend(fontsize=8.5)
    ax_line.grid(alpha=0.25)
    fig.suptitle("Notebook-aligned DGI recovery: scalar phase map -> image", fontsize=12)
    plt.tight_layout()
    fig.savefig(path, format="svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def generate(config: dict[str, Any], config_path: Path) -> dict[str, str]:
    output_dir = Path("results/notebook_aligned_recovery/dgi_stage")
    output_dir.mkdir(parents=True, exist_ok=True)
    stage = build_dgi_stage(config)

    outputs = {
        "comparison_report": output_dir / "comparison_report.json",
        "dgi_summary": output_dir / "dgi_summary.json",
        "central_dgi_lineouts": output_dir / "central_dgi_lineouts.csv",
        "metadata": output_dir / "metadata.json",
        "figure": output_dir / "dgi_image_stage.svg",
    }
    write_json(outputs["comparison_report"], comparison_report(config, stage))
    write_json(outputs["dgi_summary"], dgi_summary(config, stage))
    write_dgi_lineouts(outputs["central_dgi_lineouts"], config, stage)
    write_dgi_figure(outputs["figure"], config, stage)
    write_json(
        outputs["metadata"],
        {
            "label": config["label"],
            "status": "Canonical notebook-aligned DGI-stage recovery output.",
            "config_file_used": str(config_path),
            "git_commit_hash": git_commit(),
            "source_notebook": config["source_notebook"],
            "source_notebook_cells": {
                "dgi_transfer_curve": 16,
                "fourier_optics_sim_image": 18,
            },
            "generated_outputs": {key: str(path) for key, path in outputs.items()},
            "scope_boundary": "Only DGI recovery is generated. Faraday, camera, and multishot recovery are not performed.",
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
    print("Recovered canonical notebook DGI stage outputs:")
    for key, path in outputs.items():
        print(f"- {key}: {path}")


if __name__ == "__main__":
    main()
