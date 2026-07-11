"""Recover the notebook V1 Faraday camera-level reference panel."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from non_destructive_image import intensity_at_atoms, simulate_camera_image, simulate_noisy_camera_image
from scripts.plot_label_utils import DARK_FIELD_INTENSITY, DUAL_PORT_SIGNAL, coordinate_label
from scripts.recover_notebook_camera_stage import photon_energy_j
from scripts.recover_notebook_condensate_stage import load_config
from scripts.recover_notebook_faraday_stage import build_faraday_stage
from scripts.recover_notebook_pci_stage import (
    field_comparison,
    git_commit,
    real_array_stats,
    write_json,
    write_rows,
)


OUTPUT_DIR = Path("results/notebook_aligned_recovery/faraday_camera_panel")
CELL_51_PROBE_POWER_MW = 5.0
CELL_51_REFERENCE = "notebook_sections/06_faraday.py cell 51, section 17.2 / 17.3"


def photons_per_pixel_for_power(config: dict[str, Any], probe_power_mw: float) -> dict[str, float]:
    """Reproduce notebook N_phot_pix(P_mW) for the cell 51 probe power."""

    geometry = config["imaging_geometry"]
    camera = config["camera_recovery"]
    magnification = float(geometry["focal_length_2_m"]) / float(geometry["focal_length_1_m"])
    object_pixel_m = float(geometry["camera_pixel_m"]) / magnification
    intensity_w_m2 = intensity_at_atoms(
        probe_power_mw,
        float(geometry["probe_diameter_m"]),
        use_peak_intensity=True,
    )
    photon_energy = photon_energy_j(config)
    photons_per_pixel = (
        intensity_w_m2
        * object_pixel_m**2
        * float(camera["default_exposure_s"])
        * float(camera["quantum_efficiency"])
        / photon_energy
    )
    return {
        "probe_power_mw": float(probe_power_mw),
        "magnification": float(magnification),
        "object_pixel_m": float(object_pixel_m),
        "intensity_at_atoms_w_m2": float(intensity_w_m2),
        "photon_energy_j": float(photon_energy),
        "photons_per_pixel": float(photons_per_pixel),
    }


def camera_axis_um(config: dict[str, Any], high_res_axis_m: np.ndarray, camera_shape: tuple[int, int]) -> np.ndarray:
    bin_size = int(config["camera_recovery"]["bin_size"])
    cols = camera_shape[1] * bin_size
    conversion = float(config["unit_conventions"]["coordinate_conversion_m_to_um"])
    return (high_res_axis_m[:cols] * conversion).reshape(camera_shape[1], bin_size).mean(axis=1)


def _camera_pair(
    image: np.ndarray,
    *,
    bin_size: int,
    photons_per_pixel: float,
    read_noise_electrons: float,
    rng: np.random.Generator,
) -> dict[str, np.ndarray]:
    deterministic = simulate_camera_image(
        image,
        bin_size=bin_size,
        photons_per_pixel=photons_per_pixel,
        return_intermediates=True,
    )
    noisy = simulate_noisy_camera_image(
        deterministic["binned_image"],
        photons_per_pixel=photons_per_pixel,
        rng=rng,
        read_noise_electrons=read_noise_electrons,
        input_is_binned=True,
        normalize=True,
        return_intermediates=True,
    )
    return {
        "binned_image": deterministic["binned_image"],
        "deterministic_counts": deterministic["deterministic_counts"],
        "camera_image": deterministic["camera_image"],
        "noisy_counts": noisy["noisy_counts"],
        "noisy_image": noisy["noisy_image"],
    }


def build_faraday_camera_panel(config: dict[str, Any]) -> dict[str, Any]:
    faraday = build_faraday_stage(config)
    camera = config["camera_recovery"]
    noise = config["noisy_camera_recovery"]
    bin_size = int(camera["bin_size"])
    read_noise = float(camera["read_noise_electrons"])
    seed = int(noise["rng_seed"])
    photon_scale = photons_per_pixel_for_power(config, CELL_51_PROBE_POWER_MW)
    photons_per_pixel = photon_scale["photons_per_pixel"]

    rng = np.random.default_rng(seed)
    dark = _camera_pair(
        faraday["notebook_dark_field_intensity"],
        bin_size=bin_size,
        photons_per_pixel=photons_per_pixel,
        read_noise_electrons=read_noise,
        rng=rng,
    )
    port_u = _camera_pair(
        faraday["notebook_dual_port_u_intensity"],
        bin_size=bin_size,
        photons_per_pixel=photons_per_pixel,
        read_noise_electrons=read_noise,
        rng=rng,
    )
    port_v = _camera_pair(
        faraday["notebook_dual_port_v_intensity"],
        bin_size=bin_size,
        photons_per_pixel=photons_per_pixel,
        read_noise_electrons=read_noise,
        rng=rng,
    )
    noisy_signal = (port_v["noisy_image"] - port_u["noisy_image"]) / (
        port_v["noisy_image"] + port_u["noisy_image"]
    )
    ideal_signal = (port_v["camera_image"] - port_u["camera_image"]) / (
        port_v["camera_image"] + port_u["camera_image"]
    )

    replay_rng = np.random.default_rng(seed)
    replay_dark = simulate_noisy_camera_image(
        faraday["notebook_dark_field_intensity"],
        photons_per_pixel=photons_per_pixel,
        rng=replay_rng,
        read_noise_electrons=read_noise,
        bin_size=bin_size,
        normalize=True,
        return_intermediates=True,
    )
    replay_u = simulate_noisy_camera_image(
        faraday["notebook_dual_port_u_intensity"],
        photons_per_pixel=photons_per_pixel,
        rng=replay_rng,
        read_noise_electrons=read_noise,
        bin_size=bin_size,
        normalize=True,
        return_intermediates=True,
    )
    replay_v = simulate_noisy_camera_image(
        faraday["notebook_dual_port_v_intensity"],
        photons_per_pixel=photons_per_pixel,
        rng=replay_rng,
        read_noise_electrons=read_noise,
        bin_size=bin_size,
        normalize=True,
        return_intermediates=True,
    )
    replay_signal = (replay_v["noisy_image"] - replay_u["noisy_image"]) / (
        replay_v["noisy_image"] + replay_u["noisy_image"]
    )

    axis_um = camera_axis_um(
        config,
        faraday["phase_stage"]["condensate_stage"]["coordinate_axis_m"],
        dark["camera_image"].shape,
    )
    mid = dark["camera_image"].shape[0] // 2
    return {
        "faraday_stage": faraday,
        "photon_scale": photon_scale,
        "bin_size": bin_size,
        "read_noise_electrons": read_noise,
        "rng_seed": seed,
        "dark": dark,
        "port_u": port_u,
        "port_v": port_v,
        "dual_port_noisy_signal": noisy_signal,
        "dual_port_ideal_signal": ideal_signal,
        "replay": {
            "dark_noisy_image": replay_dark["noisy_image"],
            "port_u_noisy_image": replay_u["noisy_image"],
            "port_v_noisy_image": replay_v["noisy_image"],
            "dual_port_noisy_signal": replay_signal,
        },
        "axis_um": axis_um,
        "mid": mid,
    }


def _frame_row(name: str, array: np.ndarray, units: str, noisy: bool) -> dict[str, float | str | bool]:
    stats = real_array_stats(array)
    return {
        "quantity": name,
        "units": units,
        "is_noisy": noisy,
        "shape_rows": stats["shape"][0],
        "shape_cols": stats["shape"][1],
        "min": stats["min"],
        "max": stats["max"],
        "mean": stats["mean"],
        "sum": stats["sum"],
        "centre_value": stats["centre_value"],
        "peak_row": stats["peak_index"][0],
        "peak_col": stats["peak_index"][1],
    }


def comparison_report(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    dark = stage["dark"]
    port_u = stage["port_u"]
    port_v = stage["port_v"]
    return {
        "status": "notebook cell 51 Faraday camera panel recovered with explicit-seed replay",
        "notebook_counterpart": {
            "section": "17.2 / 17.3",
            "primary_cell": 51,
            "file": "notebook_sections/06_faraday.py",
            "reference": CELL_51_REFERENCE,
        },
        "notebook_formula": {
            "faraday_maps": "fm = faraday_maps(1.5e9, axis=0)",
            "dark_camera": "cam_dark, ideal_dark = to_camera(fm['I_dark'], 5.0)",
            "port_u_camera": "cam_u, ideal_u = to_camera(fm['I_u'], 5.0)",
            "port_v_camera": "cam_v, ideal_v = to_camera(fm['I_v'], 5.0)",
            "noisy_signal": "S_map = (cam_v - cam_u) / (cam_v + cam_u)",
            "ideal_signal": "S_ideal = (ideal_v - ideal_u) / (ideal_v + ideal_u)",
        },
        "camera_binning": {
            "bin_size": stage["bin_size"],
            "trimmed_high_resolution_shape": [1020, 1020],
            "camera_shape": list(dark["camera_image"].shape),
        },
        "photon_scale": stage["photon_scale"],
        "noise": {
            "read_noise_electrons": stage["read_noise_electrons"],
            "rng_seed": stage["rng_seed"],
            "notebook_rng_policy": "global np.random.default_rng(7)",
            "recovery_rng_policy": "explicit np.random.default_rng(seed), replaying the cell 51 dark/u/v call order",
            "exact_notebook_global_rng_reproduction": False,
            "exact_explicit_seed_replay": True,
        },
        "dark_field_noiseless_binned": {
            "units": "incident-I0-normalised intensity",
            "stats": real_array_stats(dark["camera_image"]),
        },
        "dark_field_noisy": {
            "units": "incident-I0-normalised intensity",
            "stats": real_array_stats(dark["noisy_image"]),
            **field_comparison(dark["noisy_image"], stage["replay"]["dark_noisy_image"]),
        },
        "dual_port_u_noiseless_binned": {
            "units": "incident-I0-normalised intensity",
            "stats": real_array_stats(port_u["camera_image"]),
        },
        "dual_port_v_noiseless_binned": {
            "units": "incident-I0-normalised intensity",
            "stats": real_array_stats(port_v["camera_image"]),
        },
        "dual_port_noiseless_signal": {
            "units": "normalised difference",
            "stats": real_array_stats(stage["dual_port_ideal_signal"]),
        },
        "dual_port_noisy_signal": {
            "units": "normalised difference",
            "stats": real_array_stats(stage["dual_port_noisy_signal"]),
            **field_comparison(stage["dual_port_noisy_signal"], stage["replay"]["dual_port_noisy_signal"]),
        },
        "central_lineouts": {
            "lineout_axis": "central camera row versus y coordinate",
            "row_index": stage["mid"],
            "position_units": "um",
            "dark_noisy_vs_noiseless_max_abs_difference": float(
                np.max(np.abs(dark["noisy_image"][stage["mid"], :] - dark["camera_image"][stage["mid"], :]))
            ),
            "dual_port_noisy_vs_noiseless_max_abs_difference": float(
                np.max(
                    np.abs(stage["dual_port_noisy_signal"][stage["mid"], :] - stage["dual_port_ideal_signal"][stage["mid"], :])
                )
            ),
        },
        "comparison_status": {
            "exact_helper_replay": "All tested noisy arrays match direct helper replay with the same explicit RNG seed and cell 51 call order.",
            "statistical_notebook_global_rng": "Arbitrary interactive notebook-global RNG state is not claimed.",
        },
    }


def summary(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    report = comparison_report(config, stage)
    faraday = stage["faraday_stage"]
    return {
        "label": config["label"],
        "status": "Notebook-aligned Faraday camera-level reference panel recovered.",
        "notebook_counterpart": report["notebook_counterpart"],
        "detuning_hz": faraday["phase_stage"]["detuning_hz"],
        "detuning_ghz": faraday["phase_stage"]["detuning_hz"] / 1e9,
        "kappa_F": faraday["kappa_F"],
        "theta_f_peak_rad": faraday["theta_f_peak_rad"],
        "probe_power_mw": stage["photon_scale"]["probe_power_mw"],
        "exposure_time_s": float(config["camera_recovery"]["default_exposure_s"]),
        "camera_shape": list(stage["dark"]["camera_image"].shape),
        "photons_per_pixel": stage["photon_scale"]["photons_per_pixel"],
        "read_noise_electrons": stage["read_noise_electrons"],
        "rng_seed": stage["rng_seed"],
        "dark_field_noisy_stats": real_array_stats(stage["dark"]["noisy_image"]),
        "dual_port_noisy_signal_stats": real_array_stats(stage["dual_port_noisy_signal"]),
        "calibration_status": "Version 1 notebook-aligned / uncalibrated; no experimental RAI or absorption calibration applied.",
        "model_status": "Finite-rotation propagation with kappa_F = 1.0 phenomenological placeholder; no microscopic Faraday model.",
    }


def write_lineouts(path: Path, stage: dict[str, Any]) -> None:
    mid = stage["mid"]
    write_rows(
        path,
        [
            {
                "position_um": float(stage["axis_um"][index]),
                "dark_field_noiseless_binned_i0": float(stage["dark"]["camera_image"][mid, index]),
                "dark_field_noisy_i0": float(stage["dark"]["noisy_image"][mid, index]),
                "dual_port_signal_noiseless": float(stage["dual_port_ideal_signal"][mid, index]),
                "dual_port_signal_noisy": float(stage["dual_port_noisy_signal"][mid, index]),
                "port_u_noiseless_binned_i0": float(stage["port_u"]["camera_image"][mid, index]),
                "port_v_noiseless_binned_i0": float(stage["port_v"]["camera_image"][mid, index]),
                "port_u_noisy_i0": float(stage["port_u"]["noisy_image"][mid, index]),
                "port_v_noisy_i0": float(stage["port_v"]["noisy_image"][mid, index]),
            }
            for index in range(stage["dark"]["camera_image"].shape[1])
        ],
    )


def write_frame_statistics(path: Path, stage: dict[str, Any]) -> None:
    rows = [
        _frame_row("dark_field_noiseless_binned", stage["dark"]["camera_image"], "I_dark/I0", False),
        _frame_row("dark_field_noisy", stage["dark"]["noisy_image"], "I_dark/I0", True),
        _frame_row("dual_port_u_noiseless_binned", stage["port_u"]["camera_image"], "I/I0", False),
        _frame_row("dual_port_u_noisy", stage["port_u"]["noisy_image"], "I/I0", True),
        _frame_row("dual_port_v_noiseless_binned", stage["port_v"]["camera_image"], "I/I0", False),
        _frame_row("dual_port_v_noisy", stage["port_v"]["noisy_image"], "I/I0", True),
        _frame_row("dual_port_signal_noiseless", stage["dual_port_ideal_signal"], "S", False),
        _frame_row("dual_port_signal_noisy", stage["dual_port_noisy_signal"], "S", True),
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_figure(path: Path, config: dict[str, Any], stage: dict[str, Any]) -> None:
    display = config["display"]
    axis_um = stage["axis_um"]
    mid = stage["mid"]

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "legend.fontsize": 8.5,
            "figure.dpi": 120,
        }
    )
    fig, axs = plt.subplots(2, 2, figsize=(11, 6.4))

    im0 = axs[0, 0].imshow(
        stage["dark"]["noisy_image"],
        extent=display["image_extent_um"],
        origin="lower",
        cmap="inferno",
        vmin=-0.005,
        vmax=0.05,
    )
    axs[0, 0].set_title("Dark-field camera image")
    plt.colorbar(im0, ax=axs[0, 0], fraction=0.03, label=DARK_FIELD_INTENSITY)

    im1 = axs[0, 1].imshow(
        stage["dual_port_noisy_signal"],
        extent=display["image_extent_um"],
        origin="lower",
        cmap="RdBu_r",
        vmin=-0.5,
        vmax=0.5,
    )
    axs[0, 1].set_title("Dual-port signal map")
    plt.colorbar(im1, ax=axs[0, 1], fraction=0.03, label=DUAL_PORT_SIGNAL)

    for ax in axs[0]:
        ax.set_xlim(-45, 45)
        ax.set_ylim(-12, 12)
        ax.set_xlabel(coordinate_label("y"))
        ax.set_ylabel(coordinate_label("z"))

    axs[1, 0].plot(axis_um, stage["dark"]["camera_image"][mid], "C0", lw=2, label="noiseless binned")
    axs[1, 0].plot(axis_um, stage["dark"]["noisy_image"][mid], "C0.", ms=4, alpha=0.7, label="noisy pixels")
    axs[1, 0].axhline(0, color="gray", ls=":", lw=1)
    axs[1, 0].set_xlim(-45, 45)
    axs[1, 0].set_xlabel(coordinate_label("y"))
    axs[1, 0].set_ylabel(DARK_FIELD_INTENSITY)
    axs[1, 0].set_title("Dark-field lineout")
    axs[1, 0].legend()
    axs[1, 0].grid(alpha=0.25)

    axs[1, 1].plot(axis_um, stage["dual_port_ideal_signal"][mid], "C0", lw=2, label="noiseless binned")
    axs[1, 1].plot(axis_um, stage["dual_port_noisy_signal"][mid], "C0.", ms=4, alpha=0.7, label="noisy pixels")
    axs[1, 1].axhline(0, color="gray", ls=":", lw=1)
    axs[1, 1].set_xlim(-45, 45)
    axs[1, 1].set_xlabel(coordinate_label("y"))
    axs[1, 1].set_ylabel("$S$")
    axs[1, 1].set_title("Dual-port lineout")
    axs[1, 1].legend()
    axs[1, 1].grid(alpha=0.25)

    fig.suptitle("Faraday camera reference", fontsize=12)
    plt.tight_layout()
    fig.savefig(path, format="svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def metadata(config: dict[str, Any], config_path: Path, stage: dict[str, Any], outputs: dict[str, Path]) -> dict[str, Any]:
    atom = config["atom"]
    condensate = config["condensate"]
    grid = config["grid"]
    geometry = config["imaging_geometry"]
    camera = config["camera_recovery"]
    faraday = stage["faraday_stage"]
    return {
        "label": config["label"],
        "figure_type": "notebook-aligned recovery",
        "status": "Notebook-aligned Faraday camera-level reference panel.",
        "notebook_counterpart": {
            "file": "notebook_sections/06_faraday.py",
            "primary_cell": 51,
            "section": "17.2 / 17.3",
            "context_cells": [85, 87, 89, 91],
        },
        "generation_script": "scripts/recover_notebook_faraday_camera_panel.py",
        "config_file_used": str(config_path),
        "git_commit_hash": git_commit(),
        "physical_parameters": {
            "species": atom["species"],
            "transition_wavelength_m": atom["transition_wavelength_m"],
            "natural_linewidth_rad_s": atom["natural_linewidth_rad_s"],
            "resonant_cross_section_m2": atom["resonant_cross_section_m2"],
            "atom_number": condensate["atom_number"],
            "scattering_length_bohr": condensate["scattering_length_bohr"],
            "trap_frequencies_hz": condensate["trap_frequencies_hz"],
            "temperature_k": condensate["temperature_k"],
            "detuning_hz": faraday["phase_stage"]["detuning_hz"],
            "kappa_F": faraday["kappa_F"],
            "kappa_F_status": "phenomenological Version 1 placeholder",
            "theta_f_peak_rad": faraday["theta_f_peak_rad"],
        },
        "optical_parameters": {
            "probe_power_mw": stage["photon_scale"]["probe_power_mw"],
            "probe_diameter_m": geometry["probe_diameter_m"],
            "numerical_aperture": faraday["pupil_stage"]["numerical_aperture"],
            "focal_length_1_m": geometry["focal_length_1_m"],
            "focal_length_2_m": geometry["focal_length_2_m"],
        },
        "grid_and_axis": {
            "ngrid": grid["ngrid"],
            "field_of_view_m": grid["field_of_view_m"],
            "imaging_axis": geometry["imaging_axis"],
            "display_plane": "y,z",
            "high_resolution_shape": list(faraday["notebook_dark_field_intensity"].shape),
            "camera_shape": list(stage["dark"]["camera_image"].shape),
        },
        "camera_noise_parameters": {
            "exposure_time_s": camera["default_exposure_s"],
            "quantum_efficiency": camera["quantum_efficiency"],
            "bin_size": stage["bin_size"],
            "read_noise_electrons": stage["read_noise_electrons"],
            "photons_per_pixel": stage["photon_scale"]["photons_per_pixel"],
            "rng_seed": stage["rng_seed"],
            "rng_assumptions": "Explicit np.random.default_rng(seed) replay of notebook cell 51 dark/u/v camera call order.",
        },
        "normalisation_conventions": {
            "dark_field": "camera-level I_dark/I0, noisy after Poisson photon noise plus Gaussian read noise",
            "dual_port_ports": "camera-level I/I0 for each port",
            "dual_port_signal": "S = (I_v - I_u) / (I_v + I_u)",
            "lineouts": "central camera row versus y coordinate in um",
        },
        "comparison_status": {
            "exact_explicit_seed_helper_replay": True,
            "exact_arbitrary_notebook_global_rng_reproduction": False,
            "claim": "Numerical replay is exact for the explicit-seed cell 51 call order; arbitrary interactive notebook RNG state is not claimed.",
        },
        "caveats": [
            "Version 1 notebook-aligned / uncalibrated output.",
            "kappa_F = 1.0 is a phenomenological placeholder.",
            "No microscopic Faraday model is introduced.",
            "No experimental RAI or absorption calibration has been applied.",
            "This is not a final optimised operating point.",
            "Simulation uses finite-rotation expressions; small-angle formulae are interpretive scaling language only.",
        ],
        "generated_outputs": {key: str(path) for key, path in outputs.items()},
    }


def generate(config: dict[str, Any], config_path: Path) -> dict[str, str]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    stage = build_faraday_camera_panel(config)
    outputs = {
        "figure": OUTPUT_DIR / "faraday_camera_panel.svg",
        "comparison_report": OUTPUT_DIR / "comparison_report.json",
        "summary": OUTPUT_DIR / "faraday_camera_panel_summary.json",
        "metadata": OUTPUT_DIR / "metadata.json",
        "lineouts": OUTPUT_DIR / "lineouts.csv",
        "frame_statistics": OUTPUT_DIR / "frame_statistics.csv",
    }

    write_figure(outputs["figure"], config, stage)
    write_json(outputs["comparison_report"], comparison_report(config, stage))
    write_json(outputs["summary"], summary(config, stage))
    write_lineouts(outputs["lineouts"], stage)
    write_frame_statistics(outputs["frame_statistics"], stage)
    write_json(outputs["metadata"], metadata(config, config_path, stage, outputs))
    return {key: str(path) for key, path in outputs.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config = load_config(args.config)
    outputs = generate(config, args.config)
    print("Recovered notebook Faraday camera-level reference panel outputs:")
    for key, path in outputs.items():
        print(f"- {key}: {path}")


if __name__ == "__main__":
    main()
