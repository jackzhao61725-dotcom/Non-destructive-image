"""Recover the canonical notebook V1 Faraday ideal-image stage."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from non_destructive_image import simulate_faraday_image
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
from scripts.plot_label_utils import DARK_FIELD_INTENSITY, DUAL_PORT_SIGNAL, coordinate_label, cut_label


BASELINE_PATH = Path("regression/baseline/imaging/faraday_imaging_baseline_v1.npz")


def build_faraday_stage(config: dict[str, Any]) -> dict[str, Any]:
    phase_stage = build_phase_stage(config)
    faraday_config = config["faraday_recovery"]
    pupil_stage = build_pupil(config)
    kappa_f = float(faraday_config["kappa_F"])
    theta_f_peak = kappa_f * phase_stage["notebook_phase_peak_rad"]
    theta_f_map = kappa_f * phase_stage["notebook_phase_map_rad"]

    notebook_sigma_plus_object_field = np.exp(1j * theta_f_map)
    notebook_sigma_minus_object_field = np.exp(-1j * theta_f_map)
    notebook_sigma_plus_scattered_field = notebook_sigma_plus_object_field - 1
    notebook_sigma_minus_scattered_field = notebook_sigma_minus_object_field - 1
    notebook_sigma_plus_propagated_scattered_field = np.fft.ifft2(
        np.fft.fft2(notebook_sigma_plus_scattered_field) * pupil_stage["pupil"]
    )
    notebook_sigma_minus_propagated_scattered_field = np.fft.ifft2(
        np.fft.fft2(notebook_sigma_minus_scattered_field) * pupil_stage["pupil"]
    )
    notebook_sigma_plus_field = 1 + notebook_sigma_plus_propagated_scattered_field
    notebook_sigma_minus_field = 1 + notebook_sigma_minus_propagated_scattered_field
    notebook_output_ex_field = (notebook_sigma_plus_field + notebook_sigma_minus_field) / 2
    notebook_output_ey_field = 1j * (notebook_sigma_plus_field - notebook_sigma_minus_field) / 2
    notebook_dark_field_intensity = np.abs(notebook_output_ey_field) ** 2
    notebook_dual_port_u_intensity = np.abs(notebook_output_ex_field + notebook_output_ey_field) ** 2 / 2
    notebook_dual_port_v_intensity = np.abs(notebook_output_ex_field - notebook_output_ey_field) ** 2 / 2
    notebook_dual_port_signal = (notebook_dual_port_v_intensity - notebook_dual_port_u_intensity) / (
        notebook_dual_port_v_intensity + notebook_dual_port_u_intensity
    )

    helper_result = simulate_faraday_image(
        theta_f_map,
        pupil_stage["pupil"],
        return_intermediates=True,
    )

    return {
        "phase_stage": phase_stage,
        "pupil_stage": pupil_stage,
        "kappa_F": kappa_f,
        "theta_f_peak_rad": theta_f_peak,
        "notebook_theta_f_map_rad": theta_f_map,
        "notebook_pupil": pupil_stage["pupil"],
        "notebook_sigma_plus_object_field": notebook_sigma_plus_object_field,
        "notebook_sigma_minus_object_field": notebook_sigma_minus_object_field,
        "notebook_sigma_plus_scattered_field": notebook_sigma_plus_scattered_field,
        "notebook_sigma_minus_scattered_field": notebook_sigma_minus_scattered_field,
        "notebook_sigma_plus_propagated_scattered_field": notebook_sigma_plus_propagated_scattered_field,
        "notebook_sigma_minus_propagated_scattered_field": notebook_sigma_minus_propagated_scattered_field,
        "notebook_sigma_plus_field": notebook_sigma_plus_field,
        "notebook_sigma_minus_field": notebook_sigma_minus_field,
        "notebook_output_ex_field": notebook_output_ex_field,
        "notebook_output_ey_field": notebook_output_ey_field,
        "notebook_dark_field_intensity": notebook_dark_field_intensity,
        "notebook_dual_port_u_intensity": notebook_dual_port_u_intensity,
        "notebook_dual_port_v_intensity": notebook_dual_port_v_intensity,
        "notebook_dual_port_signal": notebook_dual_port_signal,
        "helper_result": helper_result,
    }


def _baseline_comparison(stage: dict[str, Any]) -> dict[str, Any]:
    if not BASELINE_PATH.exists():
        return {"available": False, "path": str(BASELINE_PATH)}

    comparisons: dict[str, Any] = {"available": True, "path": str(BASELINE_PATH), "arrays": {}}
    mapping = {
        "theta_f_map_rad": "notebook_theta_f_map_rad",
        "sigma_plus_object_field": "notebook_sigma_plus_object_field",
        "sigma_minus_object_field": "notebook_sigma_minus_object_field",
        "sigma_plus_field": "notebook_sigma_plus_field",
        "sigma_minus_field": "notebook_sigma_minus_field",
        "output_ex_field": "notebook_output_ex_field",
        "output_ey_field": "notebook_output_ey_field",
        "dark_field_intensity": "notebook_dark_field_intensity",
        "dual_port_u_intensity": "notebook_dual_port_u_intensity",
        "dual_port_v_intensity": "notebook_dual_port_v_intensity",
        "dual_port_signal": "notebook_dual_port_signal",
    }
    with np.load(BASELINE_PATH) as baseline:
        metadata = json.loads(str(baseline["metadata_json"]))
        comparisons["metadata"] = {
            "baseline_name": metadata["baseline_name"],
            "kappa_F": metadata["kappa_F"],
            "microscopic_faraday_model": metadata["microscopic_faraday_model"],
            "notebook_reference": metadata["notebook_reference"],
        }
        for baseline_key, stage_key in mapping.items():
            comparisons["arrays"][baseline_key] = field_comparison(stage[stage_key], baseline[baseline_key])
    return comparisons


def comparison_report(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    helper = stage["helper_result"]
    return {
        "status": "canonical notebook Faraday stage matches simulate_faraday_image for tested deterministic quantities",
        "source_notebook_cells": {
            "kappa_F_placeholder": 39,
            "theta_F_peak": 43,
            "sim_faraday_fields_and_faraday_maps": 51,
            "step_10_dark_field_rotation": 85,
            "step_12_dual_port_signal": 89,
        },
        "notebook_formula": {
            "theta_f_map": "kappa_F * scalar_phase_map",
            "sigma_plus_object_field": "np.exp(+1j * theta_f_map)",
            "sigma_minus_object_field": "np.exp(-1j * theta_f_map)",
            "propagation": "1 + np.fft.ifft2(np.fft.fft2(object_field - 1) * pupil)",
            "output_ex_field": "(sigma_plus_field + sigma_minus_field) / 2",
            "output_ey_field": "1j * (sigma_plus_field - sigma_minus_field) / 2",
            "dark_field_intensity": "np.abs(output_ey_field) ** 2",
            "dual_port_signal": "(I_v - I_u) / (I_v + I_u)",
        },
        "faraday_model": {
            "kappa_F": stage["kappa_F"],
            "kappa_F_status": "Version 1 phenomenological placeholder; no experimental calibration applied",
            "microscopic_faraday_model": False,
            "theta_f_peak_rad": stage["theta_f_peak_rad"],
        },
        "pupil": {
            "shape": list(stage["notebook_pupil"].shape),
            "nonzero_pixels": int(np.count_nonzero(stage["notebook_pupil"])),
            "numerical_aperture": stage["pupil_stage"]["numerical_aperture"],
            "dgrid_m": stage["pupil_stage"]["dgrid_m"],
        },
        "theta_f_map": {
            "units": "radians",
            "stats": real_array_stats(stage["notebook_theta_f_map_rad"]),
            **field_comparison(stage["notebook_theta_f_map_rad"], helper["theta_f_map_rad"]),
        },
        "sigma_plus_object_field": {
            "stats": complex_array_stats(stage["notebook_sigma_plus_object_field"]),
            **field_comparison(stage["notebook_sigma_plus_object_field"], helper["sigma_plus_object_field"]),
        },
        "sigma_minus_object_field": {
            "stats": complex_array_stats(stage["notebook_sigma_minus_object_field"]),
            **field_comparison(stage["notebook_sigma_minus_object_field"], helper["sigma_minus_object_field"]),
        },
        "sigma_plus_field": {
            "stats": complex_array_stats(stage["notebook_sigma_plus_field"]),
            **field_comparison(stage["notebook_sigma_plus_field"], helper["sigma_plus_field"]),
        },
        "sigma_minus_field": {
            "stats": complex_array_stats(stage["notebook_sigma_minus_field"]),
            **field_comparison(stage["notebook_sigma_minus_field"], helper["sigma_minus_field"]),
        },
        "output_ex_field": {
            "stats": complex_array_stats(stage["notebook_output_ex_field"]),
            **field_comparison(stage["notebook_output_ex_field"], helper["output_ex_field"]),
        },
        "output_ey_field": {
            "stats": complex_array_stats(stage["notebook_output_ey_field"]),
            **field_comparison(stage["notebook_output_ey_field"], helper["output_ey_field"]),
        },
        "dark_field_intensity": {
            "units": "incident-I0-normalised intensity",
            "stats": real_array_stats(stage["notebook_dark_field_intensity"]),
            **field_comparison(stage["notebook_dark_field_intensity"], helper["dark_field_intensity"]),
        },
        "dual_port_u_intensity": {
            "units": "incident-I0-normalised intensity",
            "stats": real_array_stats(stage["notebook_dual_port_u_intensity"]),
            **field_comparison(stage["notebook_dual_port_u_intensity"], helper["dual_port_u_intensity"]),
        },
        "dual_port_v_intensity": {
            "units": "incident-I0-normalised intensity",
            "stats": real_array_stats(stage["notebook_dual_port_v_intensity"]),
            **field_comparison(stage["notebook_dual_port_v_intensity"], helper["dual_port_v_intensity"]),
        },
        "dual_port_signal": {
            "units": "normalised difference",
            "stats": real_array_stats(stage["notebook_dual_port_signal"]),
            **field_comparison(stage["notebook_dual_port_signal"], helper["dual_port_signal"]),
        },
        "baseline_comparison": _baseline_comparison(stage),
        "scope_boundary": "Only scalar phase map -> ideal Faraday outputs are recovered. Camera, noise, and multishot stages are not generated.",
    }


def faraday_summary(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    return {
        "label": config["label"],
        "status": "Canonical notebook-aligned Faraday-stage recovery output.",
        "source_notebook": config["source_notebook"],
        "source_notebook_cells": {
            "kappa_F_placeholder": 39,
            "theta_F_peak": 43,
            "sim_faraday_fields_and_faraday_maps": 51,
            "step_10_dark_field_rotation": 85,
            "step_12_dual_port_signal": 89,
        },
        "detuning_hz": stage["phase_stage"]["detuning_hz"],
        "scalar_phase_peak_rad": stage["phase_stage"]["notebook_phase_peak_rad"],
        "kappa_F": stage["kappa_F"],
        "theta_f_peak_rad": stage["theta_f_peak_rad"],
        "numerical_aperture": stage["pupil_stage"]["numerical_aperture"],
        "pupil_nonzero_pixels": int(np.count_nonzero(stage["notebook_pupil"])),
        "theta_f_map_stats": real_array_stats(stage["notebook_theta_f_map_rad"]),
        "dark_field_intensity_stats": real_array_stats(stage["notebook_dark_field_intensity"]),
        "dual_port_u_intensity_stats": real_array_stats(stage["notebook_dual_port_u_intensity"]),
        "dual_port_v_intensity_stats": real_array_stats(stage["notebook_dual_port_v_intensity"]),
        "dual_port_signal_stats": real_array_stats(stage["notebook_dual_port_signal"]),
        "scope_boundary": "Only ideal Faraday image outputs are recovered. No camera, noisy-frame, or multishot figures are generated.",
        "no_experimental_calibration_applied": True,
        "microscopic_faraday_model": False,
    }


def write_faraday_lineouts(path: Path, config: dict[str, Any], stage: dict[str, Any]) -> None:
    units = config["unit_conventions"]
    coordinate_um = (
        stage["phase_stage"]["condensate_stage"]["coordinate_axis_m"]
        * units["coordinate_conversion_m_to_um"]
    )
    mid = stage["notebook_theta_f_map_rad"].shape[0] // 2
    write_rows(
        path,
        [
            {
                "position_um": float(coordinate_um[index]),
                "theta_f_y_line_rad": float(stage["notebook_theta_f_map_rad"][mid, index]),
                "theta_f_z_line_rad": float(stage["notebook_theta_f_map_rad"][index, mid]),
                "dark_field_y_line_i0": float(stage["notebook_dark_field_intensity"][mid, index]),
                "dark_field_z_line_i0": float(stage["notebook_dark_field_intensity"][index, mid]),
                "dual_port_signal_y_line": float(stage["notebook_dual_port_signal"][mid, index]),
                "dual_port_signal_z_line": float(stage["notebook_dual_port_signal"][index, mid]),
            }
            for index in range(coordinate_um.size)
        ],
    )


def _write_map_and_line_figure(
    path: Path,
    *,
    config: dict[str, Any],
    stage: dict[str, Any],
    image: np.ndarray,
    y_line: np.ndarray,
    z_line: np.ndarray,
    cmap: str,
    colorbar_label: str,
    ylabel: str,
    map_title: str,
    line_title: str,
    suptitle: str,
    vmin: float | None = None,
    vmax: float | None = None,
) -> None:
    display = config["display"]
    faraday_config = config["faraday_recovery"]
    units = config["unit_conventions"]
    coordinate_um = (
        stage["phase_stage"]["condensate_stage"]["coordinate_axis_m"]
        * units["coordinate_conversion_m_to_um"]
    )

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
        figsize=tuple(faraday_config["figure_size_inches"]),
        gridspec_kw={"width_ratios": [1.5, 1]},
    )
    im = ax_map.imshow(
        image,
        extent=display["image_extent_um"],
        origin="lower",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )
    ax_map.set_xlim(*display["column_density_xlim_um"])
    ax_map.set_ylim(*display["column_density_ylim_um"])
    ax_map.set_xlabel(coordinate_label("y"))
    ax_map.set_ylabel(coordinate_label("z"))
    ax_map.set_title(map_title)
    plt.colorbar(im, ax=ax_map, fraction=0.03, label=colorbar_label)

    ax_line.plot(coordinate_um, y_line, "C1", lw=2, label=cut_label("y"))
    ax_line.plot(coordinate_um, z_line, "C0", lw=2, label=cut_label("z"))
    ax_line.axhline(0, color="gray", ls=":", lw=1)
    ax_line.set_xlim(-45, 45)
    ax_line.set_xlabel(coordinate_label())
    ax_line.set_ylabel(ylabel)
    ax_line.set_title(line_title)
    ax_line.legend(fontsize=8.5)
    ax_line.grid(alpha=0.25)
    fig.suptitle(suptitle, fontsize=12)
    plt.tight_layout()
    fig.savefig(path, format="svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_faraday_figures(output_dir: Path, config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Path]:
    faraday_config = config["faraday_recovery"]
    mid = stage["notebook_theta_f_map_rad"].shape[0] // 2
    dark_path = output_dir / "faraday_dark_field_stage.svg"
    signal_path = output_dir / "faraday_dual_port_signal_stage.svg"

    _write_map_and_line_figure(
        dark_path,
        config=config,
        stage=stage,
        image=stage["notebook_dark_field_intensity"],
        y_line=stage["notebook_dark_field_intensity"][mid, :],
        z_line=stage["notebook_dark_field_intensity"][:, mid],
        cmap=faraday_config["dark_field_cmap"],
        colorbar_label=DARK_FIELD_INTENSITY,
        ylabel=DARK_FIELD_INTENSITY,
        map_title="dark-field Faraday intensity",
        line_title="centre cuts",
        suptitle="Notebook-aligned Faraday recovery: dark-field ideal output",
        vmin=-0.005,
        vmax=0.03,
    )
    _write_map_and_line_figure(
        signal_path,
        config=config,
        stage=stage,
        image=stage["notebook_dual_port_signal"],
        y_line=stage["notebook_dual_port_signal"][mid, :],
        z_line=stage["notebook_dual_port_signal"][:, mid],
        cmap=faraday_config["dual_port_cmap"],
        colorbar_label=DUAL_PORT_SIGNAL,
        ylabel="$S$",
        map_title="dual-port Faraday signal",
        line_title="centre cuts",
        suptitle="Notebook-aligned Faraday recovery: dual-port ideal output",
        vmin=-0.45,
        vmax=0.45,
    )
    return {"dark_field_figure": dark_path, "dual_port_signal_figure": signal_path}


def generate(config: dict[str, Any], config_path: Path) -> dict[str, str]:
    output_dir = Path("results/notebook_aligned_recovery/faraday_stage")
    output_dir.mkdir(parents=True, exist_ok=True)
    stage = build_faraday_stage(config)

    outputs = {
        "comparison_report": output_dir / "comparison_report.json",
        "faraday_summary": output_dir / "faraday_summary.json",
        "central_faraday_lineouts": output_dir / "central_faraday_lineouts.csv",
        "metadata": output_dir / "metadata.json",
    }
    figures = write_faraday_figures(output_dir, config, stage)
    outputs.update(figures)

    write_json(outputs["comparison_report"], comparison_report(config, stage))
    write_json(outputs["faraday_summary"], faraday_summary(config, stage))
    write_faraday_lineouts(outputs["central_faraday_lineouts"], config, stage)
    write_json(
        outputs["metadata"],
        {
            "label": config["label"],
            "status": "Canonical notebook-aligned Faraday-stage recovery output.",
            "config_file_used": str(config_path),
            "git_commit_hash": git_commit(),
            "source_notebook": config["source_notebook"],
            "source_notebook_cells": {
                "kappa_F_placeholder": 39,
                "theta_F_peak": 43,
                "sim_faraday_fields_and_faraday_maps": 51,
                "step_10_dark_field_rotation": 85,
                "step_12_dual_port_signal": 89,
            },
            "generated_outputs": {key: str(path) for key, path in outputs.items()},
            "scope_boundary": "Only ideal Faraday recovery is generated. Camera, noisy-frame, and multishot recovery are not performed.",
            "kappa_F": stage["kappa_F"],
            "kappa_F_status": "Version 1 phenomenological placeholder; not experimentally calibrated.",
            "microscopic_faraday_model": False,
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
    print("Recovered canonical notebook Faraday stage outputs:")
    for key, path in outputs.items():
        print(f"- {key}: {path}")


if __name__ == "__main__":
    main()
