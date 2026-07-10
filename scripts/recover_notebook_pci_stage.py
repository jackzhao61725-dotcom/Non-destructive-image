"""Recover the canonical notebook V1 PCI image stage."""

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

from non_destructive_image import simulate_pci_image
from scripts.recover_notebook_condensate_stage import load_config
from scripts.recover_notebook_phase_stage import build_phase_stage
from scripts.plot_label_utils import NORMALISED_INTENSITY, coordinate_label, cut_label


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
        if isinstance(value, (complex, np.complexfloating)):
            return {"real": float(np.real(value)), "imag": float(np.imag(value))}
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


def real_array_stats(array: np.ndarray) -> dict[str, Any]:
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


def complex_array_stats(array: np.ndarray) -> dict[str, Any]:
    magnitude = np.abs(array)
    phase = np.angle(array)
    return {
        "shape": list(array.shape),
        "magnitude_min": float(np.min(magnitude)),
        "magnitude_max": float(np.max(magnitude)),
        "magnitude_mean": float(np.mean(magnitude)),
        "real_mean": float(np.mean(np.real(array))),
        "imag_mean": float(np.mean(np.imag(array))),
        "phase_min": float(np.min(phase)),
        "phase_max": float(np.max(phase)),
    }


def max_relative_difference(left: np.ndarray, right: np.ndarray) -> float:
    denominator = np.maximum(np.abs(left), np.abs(right))
    mask = denominator > 0
    if not np.any(mask):
        return 0.0
    return float(np.max(np.abs(left[mask] - right[mask]) / denominator[mask]))


def build_pupil(config: dict[str, Any]) -> dict[str, np.ndarray | float]:
    atom = config["atom"]
    grid = config["grid"]
    geometry = config["imaging_geometry"]
    ngrid = int(grid["ngrid"])
    fov_m = float(grid["field_of_view_m"])
    dgrid_m = fov_m / ngrid
    spatial_frequency_axis_m_inv = np.fft.fftfreq(ngrid, dgrid_m)
    spatial_frequency_x_m_inv, spatial_frequency_y_m_inv = np.meshgrid(
        spatial_frequency_axis_m_inv,
        spatial_frequency_axis_m_inv,
    )
    numerical_aperture = float(geometry["probe_diameter_m"]) / (2 * float(geometry["focal_length_1_m"]))
    wavelength_m = float(atom["transition_wavelength_m"])
    pupil = (
        np.sqrt(spatial_frequency_x_m_inv**2 + spatial_frequency_y_m_inv**2)
        <= numerical_aperture / wavelength_m
    ).astype(float)
    return {
        "dgrid_m": dgrid_m,
        "numerical_aperture": numerical_aperture,
        "spatial_frequency_axis_m_inv": spatial_frequency_axis_m_inv,
        "spatial_frequency_x_m_inv": spatial_frequency_x_m_inv,
        "spatial_frequency_y_m_inv": spatial_frequency_y_m_inv,
        "pupil": pupil,
    }


def build_pci_stage(config: dict[str, Any]) -> dict[str, Any]:
    phase_stage = build_phase_stage(config)
    pci_config = config["pci_recovery"]
    pupil_stage = build_pupil(config)
    phase_map = phase_stage["notebook_phase_map_rad"]
    t_p = float(pci_config["phase_plate_transmittance"])
    theta = float(pci_config["phase_plate_phase_rad"])

    notebook_object_field = np.exp(1j * phase_map)
    notebook_scattered_field = notebook_object_field - 1
    notebook_propagated_scattered_field = np.fft.ifft2(
        np.fft.fft2(notebook_scattered_field) * pupil_stage["pupil"]
    )
    notebook_reference_field = t_p * np.exp(1j * theta)
    notebook_pci_image_intensity = np.abs(notebook_reference_field + notebook_propagated_scattered_field) ** 2

    helper_result = simulate_pci_image(
        phase_map,
        pupil_stage["pupil"],
        phase_plate_transmittance=t_p,
        phase_plate_phase=theta,
        return_intermediates=True,
    )

    return {
        "phase_stage": phase_stage,
        "pupil_stage": pupil_stage,
        "phase_plate_transmittance": t_p,
        "phase_plate_phase_rad": theta,
        "plate_background_intensity": t_p**2,
        "notebook_object_field": notebook_object_field,
        "notebook_scattered_field": notebook_scattered_field,
        "notebook_pupil": pupil_stage["pupil"],
        "notebook_reference_field": notebook_reference_field,
        "notebook_propagated_scattered_field": notebook_propagated_scattered_field,
        "notebook_pci_image_intensity": notebook_pci_image_intensity,
        "helper_object_field": helper_result["object_field"],
        "helper_scattered_field": helper_result["scattered_field"],
        "helper_reference_field": helper_result["pci_reference_field"],
        "helper_propagated_scattered_field": helper_result["propagated_scattered_field"],
        "helper_pci_image_intensity": helper_result["pci_image_intensity"],
    }


def field_comparison(left: np.ndarray | complex, right: np.ndarray | complex) -> dict[str, float]:
    left_array = np.asarray(left)
    right_array = np.asarray(right)
    return {
        "max_absolute_difference": float(np.max(np.abs(left_array - right_array))),
        "max_relative_difference": max_relative_difference(left_array, right_array),
    }


def comparison_report(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    intensity = stage["notebook_pci_image_intensity"]
    return {
        "status": "canonical notebook PCI stage matches simulate_pci_image for tested deterministic quantities",
        "source_notebook_cells": {
            "pci_transfer_curve": 16,
            "fourier_optics_sim_image": 18,
        },
        "notebook_formula": {
            "object_field": "np.exp(1j * phase_map)",
            "scattered_field": "object_field - 1",
            "propagation": "np.fft.ifft2(np.fft.fft2(scattered_field) * pupil)",
            "reference_field": "t_p * np.exp(1j * theta)",
            "pci_image_intensity": "np.abs(reference_field + propagated_scattered_field) ** 2",
        },
        "pupil": {
            "shape": list(stage["notebook_pupil"].shape),
            "nonzero_pixels": int(np.count_nonzero(stage["notebook_pupil"])),
            "numerical_aperture": stage["pupil_stage"]["numerical_aperture"],
            "dgrid_m": stage["pupil_stage"]["dgrid_m"],
        },
        "phase_plate": {
            "phase_plate_transmittance": stage["phase_plate_transmittance"],
            "phase_plate_phase_rad": stage["phase_plate_phase_rad"],
            "plate_background_intensity": stage["plate_background_intensity"],
            "reference_field_notebook": stage["notebook_reference_field"],
            "reference_field_helper": stage["helper_reference_field"],
            "reference_field_absolute_difference": abs(
                stage["notebook_reference_field"] - stage["helper_reference_field"]
            ),
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
        "pci_image_intensity": {
            "units": "incident-I0-normalised intensity",
            "stats": real_array_stats(intensity),
            **field_comparison(intensity, stage["helper_pci_image_intensity"]),
        },
        "scope_boundary": "Only scalar phase map -> PCI image is recovered. DGI, Faraday, camera, and multishot stages are not generated.",
    }


def pci_summary(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    return {
        "label": config["label"],
        "status": "Canonical notebook-aligned PCI-stage recovery output.",
        "source_notebook": config["source_notebook"],
        "source_notebook_cells": {
            "pci_transfer_curve": 16,
            "fourier_optics_sim_image": 18,
        },
        "detuning_hz": stage["phase_stage"]["detuning_hz"],
        "phase_peak_rad": stage["phase_stage"]["notebook_phase_peak_rad"],
        "phase_plate_transmittance": stage["phase_plate_transmittance"],
        "phase_plate_phase_rad": stage["phase_plate_phase_rad"],
        "plate_background_intensity": stage["plate_background_intensity"],
        "numerical_aperture": stage["pupil_stage"]["numerical_aperture"],
        "pupil_nonzero_pixels": int(np.count_nonzero(stage["notebook_pupil"])),
        "pci_image_stats": real_array_stats(stage["notebook_pci_image_intensity"]),
        "scope_boundary": "Only phase_map -> PCI image is recovered. No DGI, Faraday, camera, or multishot figures are generated.",
        "no_experimental_calibration_applied": True,
    }


def write_pci_lineouts(path: Path, config: dict[str, Any], stage: dict[str, Any]) -> None:
    units = config["unit_conventions"]
    coordinate_um = (
        stage["phase_stage"]["condensate_stage"]["coordinate_axis_m"]
        * units["coordinate_conversion_m_to_um"]
    )
    intensity = stage["notebook_pci_image_intensity"]
    mid = intensity.shape[0] // 2
    write_rows(
        path,
        [
            {
                "position_um": float(coordinate_um[index]),
                "pci_intensity_y_line_incident_i0": float(intensity[mid, index]),
                "pci_intensity_z_line_incident_i0": float(intensity[index, mid]),
                "pci_signal_y_line_minus_plate_background": float(
                    intensity[mid, index] - stage["plate_background_intensity"]
                ),
                "pci_signal_z_line_minus_plate_background": float(
                    intensity[index, mid] - stage["plate_background_intensity"]
                ),
            }
            for index in range(coordinate_um.size)
        ],
    )


def write_pci_figure(path: Path, config: dict[str, Any], stage: dict[str, Any]) -> None:
    display = config["display"]
    pci_config = config["pci_recovery"]
    units = config["unit_conventions"]
    coordinate_um = (
        stage["phase_stage"]["condensate_stage"]["coordinate_axis_m"]
        * units["coordinate_conversion_m_to_um"]
    )
    intensity = stage["notebook_pci_image_intensity"]
    signal = intensity - stage["plate_background_intensity"]
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
        figsize=tuple(pci_config["figure_size_inches"]),
        gridspec_kw={"width_ratios": [1.5, 1]},
    )
    im = ax_map.imshow(
        intensity,
        extent=display["image_extent_um"],
        origin="lower",
        cmap=pci_config["map_cmap"],
    )
    ax_map.set_xlim(*display["column_density_xlim_um"])
    ax_map.set_ylim(*display["column_density_ylim_um"])
    ax_map.set_xlabel(coordinate_label("y"))
    ax_map.set_ylabel(coordinate_label("z"))
    ax_map.set_title("PCI image intensity, incident-$I_0$ convention")
    plt.colorbar(im, ax=ax_map, fraction=0.03, label=NORMALISED_INTENSITY)

    ax_line.plot(coordinate_um, intensity[mid, :], "C1", lw=2, label=cut_label("y"))
    ax_line.plot(coordinate_um, intensity[:, mid], "C0", lw=2, label=cut_label("z"))
    ax_line.axhline(stage["plate_background_intensity"], color="gray", ls=":", lw=1)
    ax_line.annotate(
        f"$t_p^2$ = {stage['plate_background_intensity']:.3f}",
        (12, stage["plate_background_intensity"] * 0.99),
        fontsize=9,
        color="gray",
    )
    ax_line.set_xlim(-45, 45)
    ax_line.set_xlabel(coordinate_label())
    ax_line.set_ylabel(NORMALISED_INTENSITY)
    ax_line.set_title("centre cuts")
    ax_line.legend(fontsize=8.5)
    ax_line.grid(alpha=0.25)
    fig.suptitle("Notebook-aligned PCI recovery: scalar phase map -> image", fontsize=12)
    plt.tight_layout()
    fig.savefig(path, format="svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def generate(config: dict[str, Any], config_path: Path) -> dict[str, str]:
    output_dir = Path("results/notebook_aligned_recovery/pci_stage")
    output_dir.mkdir(parents=True, exist_ok=True)
    stage = build_pci_stage(config)

    outputs = {
        "comparison_report": output_dir / "comparison_report.json",
        "pci_summary": output_dir / "pci_summary.json",
        "central_pci_lineouts": output_dir / "central_pci_lineouts.csv",
        "metadata": output_dir / "metadata.json",
        "figure": output_dir / "pci_image_stage.svg",
    }
    write_json(outputs["comparison_report"], comparison_report(config, stage))
    write_json(outputs["pci_summary"], pci_summary(config, stage))
    write_pci_lineouts(outputs["central_pci_lineouts"], config, stage)
    write_pci_figure(outputs["figure"], config, stage)
    write_json(
        outputs["metadata"],
        {
            "label": config["label"],
            "status": "Canonical notebook-aligned PCI-stage recovery output.",
            "config_file_used": str(config_path),
            "git_commit_hash": git_commit(),
            "source_notebook": config["source_notebook"],
            "source_notebook_cells": {
                "pci_transfer_curve": 16,
                "fourier_optics_sim_image": 18,
            },
            "generated_outputs": {key: str(path) for key, path in outputs.items()},
            "scope_boundary": "Only PCI recovery is generated. DGI, Faraday, camera, and multishot recovery are not performed.",
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
    print("Recovered canonical notebook PCI stage outputs:")
    for key, path in outputs.items():
        print(f"- {key}: {path}")


if __name__ == "__main__":
    main()
