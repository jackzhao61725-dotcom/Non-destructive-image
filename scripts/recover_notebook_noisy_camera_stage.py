"""Recover the canonical notebook V1 stochastic camera/noise stage."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from non_destructive_image import simulate_noisy_camera_image
from scripts.recover_notebook_camera_stage import build_camera_stage
from scripts.recover_notebook_condensate_stage import load_config
from scripts.recover_notebook_pci_stage import field_comparison, git_commit, real_array_stats, write_json, write_rows
from scripts.plot_label_utils import NORMALISED_INTENSITY, coordinate_label


def build_noisy_camera_stage(config: dict[str, Any]) -> dict[str, Any]:
    camera_stage = build_camera_stage(config)
    noise_config = config["noisy_camera_recovery"]
    seed = int(noise_config["rng_seed"])
    photons_per_pixel = camera_stage["photon_scale"]["photons_per_pixel"]
    read_noise = camera_stage["read_noise_electrons"]
    binned_image = camera_stage["notebook_binned_image"]

    notebook_rng = np.random.default_rng(seed)
    notebook_noisy_counts = notebook_rng.poisson(np.clip(binned_image, 0, None) * photons_per_pixel) + notebook_rng.normal(
        0,
        read_noise,
        binned_image.shape,
    )
    notebook_noisy_image = notebook_noisy_counts / photons_per_pixel

    helper_result = simulate_noisy_camera_image(
        binned_image,
        photons_per_pixel=photons_per_pixel,
        rng=np.random.default_rng(seed),
        read_noise_electrons=read_noise,
        input_is_binned=True,
        normalize=True,
        return_intermediates=True,
    )
    different_seed = simulate_noisy_camera_image(
        binned_image,
        photons_per_pixel=photons_per_pixel,
        rng=np.random.default_rng(seed + 1),
        read_noise_electrons=read_noise,
        input_is_binned=True,
        normalize=True,
        return_intermediates=True,
    )

    expected_counts = np.clip(binned_image, 0, None) * photons_per_pixel
    expected_count_variance = expected_counts + read_noise**2
    expected_image_std = np.sqrt(expected_count_variance) / photons_per_pixel
    residual = notebook_noisy_image - camera_stage["notebook_camera_image"]

    return {
        "camera_stage": camera_stage,
        "input_stage": noise_config["input_stage"],
        "rng_seed": seed,
        "different_rng_seed": seed + 1,
        "read_noise_electrons": read_noise,
        "photons_per_pixel": photons_per_pixel,
        "notebook_binned_image": binned_image,
        "notebook_deterministic_camera_image": camera_stage["notebook_camera_image"],
        "notebook_expected_counts": expected_counts,
        "notebook_expected_count_variance": expected_count_variance,
        "notebook_expected_image_std": expected_image_std,
        "notebook_noisy_counts": notebook_noisy_counts,
        "notebook_noisy_image": notebook_noisy_image,
        "notebook_residual": residual,
        "helper_binned_image": helper_result["binned_image"],
        "helper_noisy_counts": helper_result["noisy_counts"],
        "helper_noisy_image": helper_result["noisy_image"],
        "different_seed_noisy_image": different_seed["noisy_image"],
    }


def scalar_stats(values: np.ndarray) -> dict[str, float]:
    return {
        "min": float(np.min(values)),
        "max": float(np.max(values)),
        "mean": float(np.mean(values)),
        "std": float(np.std(values)),
    }


def comparison_report(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    residual = stage["notebook_residual"]
    expected_std = stage["notebook_expected_image_std"]
    return {
        "status": "canonical notebook stochastic camera recipe matches simulate_noisy_camera_image for explicit-seed deterministic replay",
        "source_notebook_cells": {
            "rng_initialisation": 0,
            "camera_helpers": 20,
        },
        "notebook_formula": {
            "counts": "rng.poisson(np.clip(binned, 0, None) * Nd) + rng.normal(0, read_e, binned.shape)",
            "noisy_image": "counts / Nd",
        },
        "rng_policy": {
            "notebook_rng": "global np.random.default_rng(7)",
            "recovery_rng": "explicit np.random.default_rng(seed)",
            "rng_seed": stage["rng_seed"],
            "exact_notebook_frame_reproduction": "Exact for this isolated first PCI camera call under clean notebook execution order; ambiguous for arbitrary interactive notebook RNG state.",
            "hidden_rng_introduced": False,
        },
        "camera_parameters": {
            "input_stage": stage["input_stage"],
            "shape": list(stage["notebook_noisy_image"].shape),
            "photons_per_pixel": stage["photons_per_pixel"],
            "read_noise_electrons": stage["read_noise_electrons"],
        },
        "binned_input": {
            "stats": real_array_stats(stage["notebook_binned_image"]),
            **field_comparison(stage["notebook_binned_image"], stage["helper_binned_image"]),
        },
        "noisy_counts": {
            "units": "detected photoelectrons after Poisson photon noise plus Gaussian read noise",
            "stats": real_array_stats(stage["notebook_noisy_counts"]),
            **field_comparison(stage["notebook_noisy_counts"], stage["helper_noisy_counts"]),
        },
        "noisy_camera_image": {
            "units": "incident-I0-normalised intensity",
            "stats": real_array_stats(stage["notebook_noisy_image"]),
            **field_comparison(stage["notebook_noisy_image"], stage["helper_noisy_image"]),
        },
        "stochastic_statistics": {
            "deterministic_mean": float(np.mean(stage["notebook_deterministic_camera_image"])),
            "noisy_mean": float(np.mean(stage["notebook_noisy_image"])),
            "mean_difference": float(np.mean(stage["notebook_noisy_image"]) - np.mean(stage["notebook_deterministic_camera_image"])),
            "residual_stats": scalar_stats(residual),
            "expected_per_pixel_image_std_stats": scalar_stats(expected_std),
            "sample_residual_std_over_mean_expected_std": float(np.std(residual) / np.mean(expected_std)),
            "different_seed_max_absolute_difference": float(
                np.max(np.abs(stage["notebook_noisy_image"] - stage["different_seed_noisy_image"]))
            ),
        },
        "scope_boundary": "Only deterministic camera image -> one stochastic noisy camera frame is recovered. Multishot frame sequences are not generated.",
    }


def noisy_camera_summary(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    report = comparison_report(config, stage)
    return {
        "label": config["label"],
        "status": "Canonical notebook-aligned stochastic camera-stage recovery output.",
        "source_notebook": config["source_notebook"],
        "source_notebook_cells": report["source_notebook_cells"],
        "rng_policy": report["rng_policy"],
        "camera_parameters": report["camera_parameters"],
        "noisy_camera_image_stats": real_array_stats(stage["notebook_noisy_image"]),
        "noisy_counts_stats": real_array_stats(stage["notebook_noisy_counts"]),
        "stochastic_statistics": report["stochastic_statistics"],
        "scope_boundary": report["scope_boundary"],
        "no_experimental_calibration_applied": True,
    }


def write_noisy_camera_statistics(path: Path, stage: dict[str, Any]) -> None:
    deterministic = stage["notebook_deterministic_camera_image"]
    noisy = stage["notebook_noisy_image"]
    counts = stage["notebook_noisy_counts"]
    residual = stage["notebook_residual"]
    expected_std = stage["notebook_expected_image_std"]
    write_rows(
        path,
        [
            {
                "frame_index": index,
                "deterministic_mean": float(np.mean(deterministic)),
                "noisy_mean": float(np.mean(noisy)),
                "noisy_counts_mean": float(np.mean(counts)),
                "residual_mean": float(np.mean(residual)),
                "residual_std": float(np.std(residual)),
                "expected_image_std_mean": float(np.mean(expected_std)),
                "expected_image_std_min": float(np.min(expected_std)),
                "expected_image_std_max": float(np.max(expected_std)),
            }
            for index in [0]
        ],
    )


def write_noisy_camera_figure(path: Path, config: dict[str, Any], stage: dict[str, Any]) -> None:
    display = config["display"]
    noise_config = config["noisy_camera_recovery"]
    units = config["unit_conventions"]
    coordinate_um = (
        stage["camera_stage"]["pci_stage"]["phase_stage"]["condensate_stage"]["coordinate_axis_m"]
        * units["coordinate_conversion_m_to_um"]
    )
    rows, cols = stage["camera_stage"]["trimmed_shape"]
    bin_size = stage["camera_stage"]["bin_size"]
    binned_axis_y_um = coordinate_um[:cols].reshape(cols // bin_size, bin_size).mean(axis=1)
    binned_axis_z_um = coordinate_um[:rows].reshape(rows // bin_size, bin_size).mean(axis=1)
    noisy = stage["notebook_noisy_image"]
    deterministic = stage["notebook_deterministic_camera_image"]
    mid_row = noisy.shape[0] // 2
    mid_col = noisy.shape[1] // 2

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
        figsize=tuple(noise_config["figure_size_inches"]),
        gridspec_kw={"width_ratios": [1.5, 1]},
    )
    im = ax_map.imshow(
        noisy,
        extent=display["image_extent_um"],
        origin="lower",
        cmap=noise_config["map_cmap"],
        vmin=noise_config["map_vmin"],
        vmax=noise_config["map_vmax"],
    )
    ax_map.set_xlim(*display["column_density_xlim_um"])
    ax_map.set_ylim(*display["column_density_ylim_um"])
    ax_map.set_xlabel(coordinate_label("y"))
    ax_map.set_ylabel(coordinate_label("z"))
    ax_map.set_title("Noisy PCI camera frame")
    plt.colorbar(im, ax=ax_map, fraction=0.03, label=NORMALISED_INTENSITY)

    ax_line.plot(binned_axis_y_um, deterministic[mid_row, :], "C0", lw=2, label="noise-free")
    ax_line.plot(binned_axis_y_um, noisy[mid_row, :], "C3.", ms=4, alpha=0.75, label="noisy frame")
    ax_line.set_xlim(-45, 45)
    ax_line.set_xlabel(coordinate_label())
    ax_line.set_ylabel(NORMALISED_INTENSITY)
    ax_line.set_title("Central lineout")
    ax_line.legend(fontsize=8.5)
    ax_line.grid(alpha=0.25)
    fig.suptitle("Noisy PCI camera frame", fontsize=12)
    plt.tight_layout()
    fig.savefig(path, format="svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def generate(config: dict[str, Any], config_path: Path) -> dict[str, str]:
    output_dir = Path("results/notebook_aligned_recovery/noisy_camera_stage")
    output_dir.mkdir(parents=True, exist_ok=True)
    stage = build_noisy_camera_stage(config)

    outputs = {
        "comparison_report": output_dir / "comparison_report.json",
        "noisy_camera_summary": output_dir / "noisy_camera_summary.json",
        "noisy_camera_statistics": output_dir / "noisy_camera_statistics.csv",
        "metadata": output_dir / "metadata.json",
        "figure": output_dir / "noisy_camera_stage.svg",
    }
    write_json(outputs["comparison_report"], comparison_report(config, stage))
    write_json(outputs["noisy_camera_summary"], noisy_camera_summary(config, stage))
    write_noisy_camera_statistics(outputs["noisy_camera_statistics"], stage)
    write_noisy_camera_figure(outputs["figure"], config, stage)
    write_json(
        outputs["metadata"],
        {
            "label": config["label"],
            "status": "Canonical notebook-aligned stochastic camera-stage recovery output.",
            "config_file_used": str(config_path),
            "git_commit_hash": git_commit(),
            "source_notebook": config["source_notebook"],
            "source_notebook_cells": {
                "rng_initialisation": 0,
                "camera_helpers": 20,
            },
            "generated_outputs": {key: str(path) for key, path in outputs.items()},
            "scope_boundary": "Only stochastic camera recovery is generated. Multishot recovery is not performed.",
            "rng_seed": stage["rng_seed"],
            "read_noise_electrons": stage["read_noise_electrons"],
            "photons_per_pixel": stage["photons_per_pixel"],
            "exact_notebook_frame_reproduction": "Exact for this isolated first PCI camera call under clean notebook execution order; ambiguous for arbitrary interactive notebook RNG state.",
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
    print("Recovered canonical notebook stochastic camera stage outputs:")
    for key, path in outputs.items():
        print(f"- {key}: {path}")


if __name__ == "__main__":
    main()
