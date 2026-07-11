"""Recover the canonical notebook V1 deterministic camera stage."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from non_destructive_image import intensity_at_atoms, simulate_camera_image
from scripts.recover_notebook_condensate_stage import load_config
from scripts.recover_notebook_pci_stage import (
    field_comparison,
    git_commit,
    real_array_stats,
    write_json,
    write_rows,
)
from scripts.recover_notebook_pci_stage import build_pci_stage
from scripts.plot_label_utils import NORMALISED_INTENSITY, coordinate_label, cut_label


def photon_energy_j(config: dict[str, Any]) -> float:
    constants = config["constants"]
    atom = config["atom"]
    h_planck = 2 * np.pi * float(constants["hbar"])
    return h_planck * float(constants["speed_of_light"]) / float(atom["transition_wavelength_m"])


def notebook_photons_per_pixel(config: dict[str, Any]) -> dict[str, float]:
    geometry = config["imaging_geometry"]
    camera_config = config["camera_recovery"]
    magnification = float(geometry["focal_length_2_m"]) / float(geometry["focal_length_1_m"])
    object_pixel_m = float(geometry["camera_pixel_m"]) / magnification
    intensity_w_m2 = intensity_at_atoms(
        float(camera_config["probe_power_mw"]),
        float(geometry["probe_diameter_m"]),
        use_peak_intensity=True,
    )
    photons_per_pixel = (
        intensity_w_m2
        * object_pixel_m**2
        * float(camera_config["default_exposure_s"])
        * float(camera_config["quantum_efficiency"])
        / photon_energy_j(config)
    )
    return {
        "magnification": magnification,
        "object_pixel_m": object_pixel_m,
        "intensity_at_atoms_w_m2": float(intensity_w_m2),
        "photon_energy_j": photon_energy_j(config),
        "photons_per_pixel": float(photons_per_pixel),
    }


def build_camera_stage(config: dict[str, Any]) -> dict[str, Any]:
    camera_config = config["camera_recovery"]
    pci_stage = build_pci_stage(config)
    input_ideal_image = pci_stage["notebook_pci_image_intensity"]
    bin_size = int(camera_config["bin_size"])
    rows = (input_ideal_image.shape[0] // bin_size) * bin_size
    cols = (input_ideal_image.shape[1] // bin_size) * bin_size
    notebook_binned_image = (
        input_ideal_image[:rows, :cols]
        .reshape(rows // bin_size, bin_size, cols // bin_size, bin_size)
        .mean(axis=(1, 3))
    )
    photon_scale = notebook_photons_per_pixel(config)
    photons_per_pixel = photon_scale["photons_per_pixel"]
    notebook_deterministic_counts = notebook_binned_image * photons_per_pixel
    notebook_camera_image = notebook_deterministic_counts / photons_per_pixel

    helper_result = simulate_camera_image(
        input_ideal_image,
        bin_size=bin_size,
        photons_per_pixel=photons_per_pixel,
        return_intermediates=True,
    )

    return {
        "pci_stage": pci_stage,
        "input_stage": camera_config["input_stage"],
        "input_ideal_image": input_ideal_image,
        "probe_power_mw": float(camera_config["probe_power_mw"]),
        "bin_size": bin_size,
        "trimmed_shape": [rows, cols],
        "camera_shape": list(notebook_binned_image.shape),
        "default_exposure_s": float(camera_config["default_exposure_s"]),
        "quantum_efficiency": float(camera_config["quantum_efficiency"]),
        "read_noise_electrons": float(camera_config["read_noise_electrons"]),
        "photon_scale": photon_scale,
        "notebook_binned_image": notebook_binned_image,
        "notebook_deterministic_counts": notebook_deterministic_counts,
        "notebook_camera_image": notebook_camera_image,
        "helper_binned_image": helper_result["binned_image"],
        "helper_deterministic_counts": helper_result["deterministic_counts"],
        "helper_camera_image": helper_result["camera_image"],
    }


def comparison_report(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    return {
        "status": "canonical notebook deterministic camera stage matches simulate_camera_image for tested quantities",
        "source_notebook_cells": {
            "camera_helpers": 20,
        },
        "notebook_formula": {
            "magnification": "Mag = f2 / f1",
            "object_pixel": "pix_obj = pix_cam / Mag",
            "photons_per_pixel": "intensity_at_atoms(P_mW) * pix_obj**2 * tau_s * QE / E_phot",
            "binning": "Iratio[:nb, :nb].reshape(nb//15, 15, nb//15, 15).mean(axis=(1, 3))",
            "deterministic_counts": "binned * N_phot_pix",
            "camera_image": "deterministic_counts / N_phot_pix",
        },
        "scope_boundary": "Only the deterministic camera pipeline is recovered. Poisson photon noise, Gaussian read noise, noisy frames, and multishot sequences are not generated.",
        "input_ideal_image": {
            "source": "recovered PCI image intensity",
            "units": "incident-I0-normalised intensity",
            "stats": real_array_stats(stage["input_ideal_image"]),
        },
        "camera_parameters": {
            "probe_power_mw": stage["probe_power_mw"],
            "bin_size": stage["bin_size"],
            "trimmed_shape": stage["trimmed_shape"],
            "camera_shape": stage["camera_shape"],
            "default_exposure_s": stage["default_exposure_s"],
            "quantum_efficiency": stage["quantum_efficiency"],
            "read_noise_electrons": stage["read_noise_electrons"],
            "read_noise_status": "recorded from notebook but not applied in deterministic recovery",
            **stage["photon_scale"],
        },
        "binned_image": {
            "units": "incident-I0-normalised intensity",
            "stats": real_array_stats(stage["notebook_binned_image"]),
            **field_comparison(stage["notebook_binned_image"], stage["helper_binned_image"]),
        },
        "deterministic_counts": {
            "units": "detected photoelectrons before stochastic noise",
            "stats": real_array_stats(stage["notebook_deterministic_counts"]),
            **field_comparison(stage["notebook_deterministic_counts"], stage["helper_deterministic_counts"]),
        },
        "normalised_camera_image": {
            "units": "incident-I0-normalised intensity",
            "stats": real_array_stats(stage["notebook_camera_image"]),
            **field_comparison(stage["notebook_camera_image"], stage["helper_camera_image"]),
        },
        "normalisation_identity": {
            "description": "With no stochastic noise, counts / N_phot_pix returns the binned ideal image.",
            **field_comparison(stage["notebook_camera_image"], stage["notebook_binned_image"]),
        },
    }


def camera_summary(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    return {
        "label": config["label"],
        "status": "Canonical notebook-aligned deterministic camera-stage recovery output.",
        "source_notebook": config["source_notebook"],
        "source_notebook_cells": {"camera_helpers": 20},
        "input_stage": stage["input_stage"],
        "camera_parameters": {
            "probe_power_mw": stage["probe_power_mw"],
            "bin_size": stage["bin_size"],
            "default_exposure_s": stage["default_exposure_s"],
            "quantum_efficiency": stage["quantum_efficiency"],
            "read_noise_electrons": stage["read_noise_electrons"],
            "read_noise_applied": False,
            **stage["photon_scale"],
        },
        "input_ideal_image_stats": real_array_stats(stage["input_ideal_image"]),
        "binned_image_stats": real_array_stats(stage["notebook_binned_image"]),
        "deterministic_counts_stats": real_array_stats(stage["notebook_deterministic_counts"]),
        "normalised_camera_image_stats": real_array_stats(stage["notebook_camera_image"]),
        "scope_boundary": "Only ideal image -> deterministic camera image is recovered. No stochastic noisy frame or multishot sequence is generated.",
        "no_experimental_calibration_applied": True,
    }


def write_camera_lineouts(path: Path, config: dict[str, Any], stage: dict[str, Any]) -> None:
    units = config["unit_conventions"]
    coordinate_um = (
        stage["pci_stage"]["phase_stage"]["condensate_stage"]["coordinate_axis_m"]
        * units["coordinate_conversion_m_to_um"]
    )
    bin_size = stage["bin_size"]
    rows, cols = stage["trimmed_shape"]
    binned_axis_y_um = coordinate_um[:cols].reshape(cols // bin_size, bin_size).mean(axis=1)
    binned_axis_z_um = coordinate_um[:rows].reshape(rows // bin_size, bin_size).mean(axis=1)
    image = stage["notebook_camera_image"]
    counts = stage["notebook_deterministic_counts"]
    mid_row = image.shape[0] // 2
    mid_col = image.shape[1] // 2
    write_rows(
        path,
        [
            {
                "position_um": float(binned_axis_y_um[index]),
                "camera_y_line_normalised_i0": float(image[mid_row, index]),
                "camera_z_line_normalised_i0": float(image[index, mid_col]),
                "binned_y_line_normalised_i0": float(stage["notebook_binned_image"][mid_row, index]),
                "binned_z_line_normalised_i0": float(stage["notebook_binned_image"][index, mid_col]),
                "deterministic_counts_y_line_electrons": float(counts[mid_row, index]),
                "deterministic_counts_z_line_electrons": float(counts[index, mid_col]),
            }
            for index in range(image.shape[0])
        ],
    )


def write_camera_figure(path: Path, config: dict[str, Any], stage: dict[str, Any]) -> None:
    display = config["display"]
    camera_config = config["camera_recovery"]
    units = config["unit_conventions"]
    coordinate_um = (
        stage["pci_stage"]["phase_stage"]["condensate_stage"]["coordinate_axis_m"]
        * units["coordinate_conversion_m_to_um"]
    )
    rows, cols = stage["trimmed_shape"]
    bin_size = stage["bin_size"]
    binned_axis_y_um = coordinate_um[:cols].reshape(cols // bin_size, bin_size).mean(axis=1)
    binned_axis_z_um = coordinate_um[:rows].reshape(rows // bin_size, bin_size).mean(axis=1)
    image = stage["notebook_camera_image"]
    mid_row = image.shape[0] // 2
    mid_col = image.shape[1] // 2

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
        figsize=tuple(camera_config["figure_size_inches"]),
        gridspec_kw={"width_ratios": [1.5, 1]},
    )
    im = ax_map.imshow(
        image,
        extent=display["image_extent_um"],
        origin="lower",
        cmap=camera_config["map_cmap"],
    )
    ax_map.set_xlim(*display["column_density_xlim_um"])
    ax_map.set_ylim(*display["column_density_ylim_um"])
    ax_map.set_xlabel(coordinate_label("y"))
    ax_map.set_ylabel(coordinate_label("z"))
    ax_map.set_title("Camera image")
    plt.colorbar(im, ax=ax_map, fraction=0.03, label=NORMALISED_INTENSITY)

    ax_line.plot(binned_axis_y_um, image[mid_row, :], "C1", lw=2, label=f"camera {cut_label('y')}")
    ax_line.plot(binned_axis_z_um, image[:, mid_col], "C0", lw=2, label=f"camera {cut_label('z')}")
    ax_line.set_xlim(-45, 45)
    ax_line.set_xlabel(coordinate_label())
    ax_line.set_ylabel(NORMALISED_INTENSITY)
    ax_line.set_title("Central lineouts")
    ax_line.legend(fontsize=8.5)
    ax_line.grid(alpha=0.25)
    fig.suptitle("Deterministic camera image", fontsize=12)
    plt.tight_layout()
    fig.savefig(path, format="svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def generate(config: dict[str, Any], config_path: Path) -> dict[str, str]:
    output_dir = Path("results/notebook_aligned_recovery/camera_stage")
    output_dir.mkdir(parents=True, exist_ok=True)
    stage = build_camera_stage(config)

    outputs = {
        "comparison_report": output_dir / "comparison_report.json",
        "camera_summary": output_dir / "camera_summary.json",
        "central_camera_lineouts": output_dir / "central_camera_lineouts.csv",
        "metadata": output_dir / "metadata.json",
        "figure": output_dir / "camera_deterministic_stage.svg",
    }
    write_json(outputs["comparison_report"], comparison_report(config, stage))
    write_json(outputs["camera_summary"], camera_summary(config, stage))
    write_camera_lineouts(outputs["central_camera_lineouts"], config, stage)
    write_camera_figure(outputs["figure"], config, stage)
    write_json(
        outputs["metadata"],
        {
            "label": config["label"],
            "status": "Canonical notebook-aligned deterministic camera-stage recovery output.",
            "config_file_used": str(config_path),
            "git_commit_hash": git_commit(),
            "source_notebook": config["source_notebook"],
            "source_notebook_cells": {"camera_helpers": 20},
            "generated_outputs": {key: str(path) for key, path in outputs.items()},
            "scope_boundary": "Only deterministic camera recovery is generated. Stochastic noisy frames and multishot recovery are not performed.",
            "input_stage": stage["input_stage"],
            "probe_power_mw": stage["probe_power_mw"],
            "bin_size": stage["bin_size"],
            "photons_per_pixel": stage["photon_scale"]["photons_per_pixel"],
            "read_noise_applied": False,
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
    print("Recovered canonical notebook deterministic camera stage outputs:")
    for key, path in outputs.items():
        print(f"- {key}: {path}")


if __name__ == "__main__":
    main()
