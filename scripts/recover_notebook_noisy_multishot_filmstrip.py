"""Recover the notebook-aligned noisy PCI multishot filmstrip."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from non_destructive_image import simulate_noisy_camera_image
from scripts.recover_notebook_condensate_stage import load_config
from scripts.recover_notebook_multishot_stage import build_multishot_stage, photons_per_camera_pixel
from scripts.recover_notebook_pci_stage import build_pupil, git_commit, real_array_stats, write_json, write_rows


def _fresh_profile(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    grid = config["grid"]
    axis = int(config["multishot_recovery"]["imaging_axis"])
    constants = stage["constants"]
    ngrid = int(grid["ngrid"])
    fov = float(grid["field_of_view_m"])
    dgrid = fov / ngrid
    coordinate_axis = (np.arange(ngrid) - ngrid // 2) * dgrid
    grid_a, grid_b = np.meshgrid(coordinate_axis, coordinate_axis)
    plane = [index for index in range(3) if index != axis]
    radii = constants["radii"]
    profile = np.maximum(0, 1 - grid_a**2 / radii[plane[0]] ** 2 - grid_b**2 / radii[plane[1]] ** 2) ** 1.5
    return {"coordinate_axis_m": coordinate_axis, "profile": profile, "plane": plane}


def _ideal_pci_image(config: dict[str, Any], profile: np.ndarray, phase_value: float, pupil: np.ndarray) -> np.ndarray:
    pci = config["pci_recovery"]
    t_p = float(pci["phase_plate_transmittance"])
    theta = float(pci["phase_plate_phase_rad"])
    scattered = np.fft.ifft2(np.fft.fft2(np.exp(1j * phase_value * profile) - 1) * pupil)
    return np.abs(t_p * np.exp(1j * theta) + scattered) ** 2


def _notebook_noisy_frame(
    binned_image: np.ndarray,
    photons_per_pixel: float,
    read_noise: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    counts = rng.poisson(np.clip(binned_image, 0, None) * photons_per_pixel) + rng.normal(
        0,
        read_noise,
        binned_image.shape,
    )
    return counts / photons_per_pixel, counts


def _bin_image(image: np.ndarray, bin_size: int) -> np.ndarray:
    rows = (image.shape[0] // bin_size) * bin_size
    cols = (image.shape[1] // bin_size) * bin_size
    return image[:rows, :cols].reshape(rows // bin_size, bin_size, cols // bin_size, bin_size).mean(axis=(1, 3))


def build_noisy_multishot_filmstrip(config: dict[str, Any]) -> dict[str, Any]:
    multishot_stage = build_multishot_stage(config)
    filmstrip_config = config["noisy_multishot_filmstrip_recovery"]
    sequence = multishot_stage["notebook"]["heating"]
    selected = [int(frame) for frame in filmstrip_config["selected_frames"]]
    if any(frame >= len(sequence["shot"]) for frame in selected):
        raise ValueError("selected filmstrip frame is outside the recovered heating sequence")

    camera = config["camera_recovery"]
    multishot = config["multishot_recovery"]
    bin_size = int(camera["bin_size"])
    read_noise = float(camera["read_noise_electrons"])
    power_mw = float(multishot["probe_power_mw"])
    tau_s = float(multishot["pulse_duration_us"]) * 1e-6
    photons = photons_per_camera_pixel(config, power_mw, tau_s)
    pupil = build_pupil(config)["pupil"]
    profile_stage = _fresh_profile(config, multishot_stage)
    profile = profile_stage["profile"]
    rng_seed = int(filmstrip_config["rng_seed"])
    notebook_rng = np.random.default_rng(rng_seed)

    frames: list[dict[str, Any]] = []
    for frame_index in selected:
        phase = float(sequence["phi"][frame_index])
        ideal_image = _ideal_pci_image(config, profile, phase, pupil)
        binned = _bin_image(ideal_image, bin_size)
        noisy, counts = _notebook_noisy_frame(binned, photons, read_noise, notebook_rng)
        helper = simulate_noisy_camera_image(
            binned,
            photons_per_pixel=photons,
            rng=np.random.default_rng(rng_seed),
            read_noise_electrons=read_noise,
            input_is_binned=True,
            normalize=True,
            return_intermediates=True,
        )
        # Advance the helper RNG to the same point by replaying all frames up to this one.
        replay_rng = np.random.default_rng(rng_seed)
        replay_result = None
        for replay_frame in selected[: selected.index(frame_index) + 1]:
            replay_phase = float(sequence["phi"][replay_frame])
            replay_ideal = _ideal_pci_image(config, profile, replay_phase, pupil)
            replay_binned = _bin_image(replay_ideal, bin_size)
            replay_result = simulate_noisy_camera_image(
                replay_binned,
                photons_per_pixel=photons,
                rng=replay_rng,
                read_noise_electrons=read_noise,
                input_is_binned=True,
                normalize=True,
                return_intermediates=True,
            )
        assert replay_result is not None
        frames.append(
            {
                "frame_index": frame_index,
                "shot": float(sequence["shot"][frame_index]),
                "N0": float(sequence["N0"][frame_index]),
                "loss_fraction": float(sequence["frac"][frame_index]),
                "temperature_k": float(sequence["T"][frame_index]),
                "phi_rad": phase,
                "snr": float(sequence["snr"][frame_index]),
                "accumulated_snr": float(sequence["accumulated_snr"][frame_index]),
                "ideal_image": ideal_image,
                "binned_image": binned,
                "noisy_frame": noisy,
                "noisy_counts": counts,
                "helper_binned_image": replay_result["binned_image"],
                "helper_noisy_frame": replay_result["noisy_image"],
                "helper_noisy_counts": replay_result["noisy_counts"],
                "single_frame_helper_check_not_used_for_comparison": helper["noisy_image"],
            }
        )

    cell44_last = int(len(sequence["shot"]) - 1)
    cell44_frames = sorted(set([0, max(1, cell44_last // 2), cell44_last]))
    return {
        "multishot_stage": multishot_stage,
        "profile_stage": profile_stage,
        "selected_frames": selected,
        "cell44_related_frames": cell44_frames,
        "rng_seed": rng_seed,
        "photons_per_pixel": photons,
        "read_noise_electrons": read_noise,
        "bin_size": bin_size,
        "frames": frames,
    }


def _frame_comparison(frame: dict[str, Any]) -> dict[str, Any]:
    def compare(left: np.ndarray, right: np.ndarray) -> dict[str, float]:
        return {
            "max_absolute_difference": float(np.max(np.abs(left - right))),
            "max_relative_difference": float(
                np.max(np.abs(left - right) / np.maximum(np.maximum(np.abs(left), np.abs(right)), 1e-300))
            ),
        }

    return {
        "frame_index": frame["frame_index"],
        "binned_image": compare(frame["binned_image"], frame["helper_binned_image"]),
        "noisy_counts": compare(frame["noisy_counts"], frame["helper_noisy_counts"]),
        "noisy_frame": compare(frame["noisy_frame"], frame["helper_noisy_frame"]),
        "noisy_frame_stats": real_array_stats(frame["noisy_frame"]),
        "binned_image_stats": real_array_stats(frame["binned_image"]),
        "residual_stats": {
            "mean": float(np.mean(frame["noisy_frame"] - frame["binned_image"])),
            "std": float(np.std(frame["noisy_frame"] - frame["binned_image"])),
            "min": float(np.min(frame["noisy_frame"] - frame["binned_image"])),
            "max": float(np.max(frame["noisy_frame"] - frame["binned_image"])),
        },
    }


def comparison_report(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    return {
        "status": "notebook-aligned noisy PCI multishot filmstrip recovered for the tested explicit-seed recipe",
        "source_notebook_cells": {
            "cloud_fading_camera_variant": 44,
            "step_14_noisy_filmstrip": 93,
        },
        "primary_recovered_cell": 93,
        "related_cell44_frame_indices": stage["cell44_related_frames"],
        "selected_frame_indices": stage["selected_frames"],
        "notebook_formula": {
            "frame_selection": "st_show = [0, 5, 10, 14]",
            "phase_source": "phi_s = seq_h['phi'][s]",
            "ideal_image": "I_s, _ = sim_image(SEQ_axis, phi_s, 'PCI')",
            "binning": "I_s[:_nb2,:_nb2].reshape(_nb2//15,15,_nb2//15,15).mean(axis=(1,3))",
            "noise": "(rng.poisson(clip(b,0,None)*st_Nd_run) + rng.normal(0, read_e, b.shape)) / st_Nd_run",
            "display": "imshow(frame, extent=ext, origin='lower', cmap='inferno', vmin=0.90, vmax=1.18)",
        },
        "rng_policy": {
            "notebook_rng": "global np.random.default_rng(7)",
            "recovery_rng": "explicit np.random.default_rng(seed)",
            "rng_seed": stage["rng_seed"],
            "exact_full_notebook_rng_reproduction": False,
            "reason": "Cell 93 appears after earlier stochastic camera and detuning-strip draws, so exact global notebook state depends on full interactive execution history.",
            "tested_reproducibility": "The recovered filmstrip is exactly reproducible under the explicit-seed replay used by this script.",
        },
        "parameters": {
            "photons_per_pixel": stage["photons_per_pixel"],
            "read_noise_electrons": stage["read_noise_electrons"],
            "bin_size": stage["bin_size"],
            "frame_shape": list(stage["frames"][0]["noisy_frame"].shape),
        },
        "frame_comparisons": [_frame_comparison(frame) for frame in stage["frames"]],
        "scope_boundary": "Only the noisy PCI multishot filmstrip is recovered. Faraday camera panels, dual-port flicker, operating maps, and shot-noise maps are not generated.",
    }


def filmstrip_summary(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    return {
        "label": config["label"],
        "status": "Notebook-aligned noisy PCI multishot filmstrip recovery output.",
        "source_notebook": config["source_notebook"],
        "source_notebook_cells": {
            "cloud_fading_camera_variant": 44,
            "step_14_noisy_filmstrip": 93,
        },
        "selected_frame_indices": stage["selected_frames"],
        "cell44_related_frame_indices": stage["cell44_related_frames"],
        "rng_seed": stage["rng_seed"],
        "photons_per_pixel": stage["photons_per_pixel"],
        "read_noise_electrons": stage["read_noise_electrons"],
        "frame_summaries": [
            {
                "frame_index": frame["frame_index"],
                "N0": frame["N0"],
                "loss_fraction": frame["loss_fraction"],
                "temperature_k": frame["temperature_k"],
                "phi_rad": frame["phi_rad"],
                "snr": frame["snr"],
                "accumulated_snr": frame["accumulated_snr"],
                "noisy_frame_stats": real_array_stats(frame["noisy_frame"]),
            }
            for frame in stage["frames"]
        ],
        "scope_boundary": "No noisy Faraday, dual-port, operating-map, or shot-noise figures are generated.",
        "no_experimental_calibration_applied": True,
    }


def write_selected_frames(path: Path, stage: dict[str, Any]) -> None:
    write_rows(
        path,
        [
            {
                "frame_index": float(frame["frame_index"]),
                "shot": float(frame["shot"]),
                "N0": float(frame["N0"]),
                "loss_fraction": float(frame["loss_fraction"]),
                "temperature_k": float(frame["temperature_k"]),
                "phi_rad": float(frame["phi_rad"]),
                "snr": float(frame["snr"]),
                "accumulated_snr": float(frame["accumulated_snr"]),
            }
            for frame in stage["frames"]
        ],
    )


def write_frame_statistics(path: Path, stage: dict[str, Any]) -> None:
    rows = []
    for frame in stage["frames"]:
        noisy = frame["noisy_frame"]
        binned = frame["binned_image"]
        residual = noisy - binned
        rows.append(
            {
                "frame_index": float(frame["frame_index"]),
                "noisy_min": float(np.min(noisy)),
                "noisy_max": float(np.max(noisy)),
                "noisy_mean": float(np.mean(noisy)),
                "noisy_std": float(np.std(noisy)),
                "binned_mean": float(np.mean(binned)),
                "residual_mean": float(np.mean(residual)),
                "residual_std": float(np.std(residual)),
                "residual_min": float(np.min(residual)),
                "residual_max": float(np.max(residual)),
            }
        )
    write_rows(path, rows)


def write_filmstrip_figure(path: Path, config: dict[str, Any], stage: dict[str, Any]) -> None:
    filmstrip_config = config["noisy_multishot_filmstrip_recovery"]
    display = config["display"]
    frames = stage["frames"]
    fig, axes = plt.subplots(1, len(frames), figsize=tuple(filmstrip_config["figure_size_inches"]))
    if len(frames) == 1:
        axes = [axes]
    for axis, frame in zip(axes, frames):
        im = axis.imshow(
            frame["noisy_frame"],
            extent=display["image_extent_um"],
            origin="lower",
            cmap=filmstrip_config["map_cmap"],
            vmin=float(filmstrip_config["map_vmin"]),
            vmax=float(filmstrip_config["map_vmax"]),
        )
        axis.set_xlim(*display["column_density_xlim_um"])
        axis.set_ylim(*display["column_density_ylim_um"])
        axis.set_xlabel("y (um)")
        axis.set_title(
            f"shot {frame['frame_index']}\n"
            f"$N_0$ = {frame['N0'] / 1e3:.1f}k   T = {frame['temperature_k'] * 1e9:.0f} nK\n"
            f"$\\varphi$ = {frame['phi_rad']:.2f}",
            fontsize=9.5,
        )
    axes[0].set_ylabel("z (um)")
    plt.colorbar(im, ax=axes[-1], fraction=0.04, label="$I/I_0$")
    params = config["multishot_recovery"]
    fig.suptitle(
        "Step 14 - the same condensate across the run "
        f"(PCI, $\\Delta$ = {params['detuning_ghz']} GHz, {params['probe_power_mw']} mW, "
        f"{params['pulse_duration_us']:.0f} us; stops at {params['loss_fraction_limit']:.0%} loss)",
        fontsize=12,
    )
    plt.tight_layout()
    fig.savefig(path, format="svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def generate(config: dict[str, Any], config_path: Path) -> dict[str, str]:
    output_dir = Path("results/notebook_aligned_recovery/noisy_multishot_filmstrip")
    output_dir.mkdir(parents=True, exist_ok=True)
    stage = build_noisy_multishot_filmstrip(config)
    outputs = {
        "comparison_report": output_dir / "comparison_report.json",
        "filmstrip_summary": output_dir / "filmstrip_summary.json",
        "selected_frames": output_dir / "selected_frames.csv",
        "frame_statistics": output_dir / "frame_statistics.csv",
        "metadata": output_dir / "metadata.json",
        "figure": output_dir / "noisy_multishot_pci_filmstrip.svg",
    }
    write_json(outputs["comparison_report"], comparison_report(config, stage))
    write_json(outputs["filmstrip_summary"], filmstrip_summary(config, stage))
    write_selected_frames(outputs["selected_frames"], stage)
    write_frame_statistics(outputs["frame_statistics"], stage)
    write_filmstrip_figure(outputs["figure"], config, stage)
    write_json(
        outputs["metadata"],
        {
            "label": config["label"],
            "status": "Notebook-aligned noisy PCI multishot filmstrip recovery output.",
            "config_file_used": str(config_path),
            "git_commit_hash": git_commit(),
            "source_notebook": config["source_notebook"],
            "source_notebook_cells": {
                "cloud_fading_camera_variant": 44,
                "step_14_noisy_filmstrip": 93,
            },
            "generated_outputs": {key: str(path) for key, path in outputs.items()},
            "selected_frame_indices": stage["selected_frames"],
            "cell44_related_frame_indices": stage["cell44_related_frames"],
            "rng_seed": stage["rng_seed"],
            "exact_full_notebook_rng_reproduction": False,
            "scope_boundary": "Only noisy PCI filmstrip recovery is generated. Faraday camera panels, dual-port flicker, operating maps, and shot-noise maps are not generated.",
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
    print("Recovered notebook-aligned noisy PCI multishot filmstrip outputs:")
    for key, path in outputs.items():
        print(f"- {key}: {path}")


if __name__ == "__main__":
    main()
