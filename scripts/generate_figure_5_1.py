"""Generate Figure 5.1: Faraday SNR versus F with aligned noisy camera frames."""

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
from matplotlib.colors import TwoSlopeNorm
import numpy as np

from non_destructive_image import (
    centered_camera_shape,
    resample_to_camera_pixels,
    scalar_phase_shift,
    simulate_faraday_image,
    simulate_noisy_camera_image,
)
from scripts.recover_notebook_condensate_stage import load_config
from scripts.recover_notebook_multishot_stage import _basic_constants, intensity_at_atoms_notebook
from scripts.recover_notebook_pci_stage import build_pupil
from scripts.recover_notebook_phase_stage import build_phase_stage


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = REPO_ROOT / "configs" / "figure_5_1.json"
MODES = ("Faraday dark-field", "Faraday dual-port")
MODE_COLOURS = {
    "Faraday dark-field": "#C0842A",
    "Faraday dual-port": "#2B7A73",
}


def _load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _git_value(*args: str) -> str:
    for command in (["git", *args], [r"C:\Program Files\Git\cmd\git.exe", *args]):
        try:
            return subprocess.check_output(
                command,
                cwd=REPO_ROOT,
                text=True,
                stderr=subprocess.DEVNULL,
            ).strip()
        except (OSError, subprocess.CalledProcessError):
            continue
    return "unknown"


def _caption(config: dict[str, Any]) -> str:
    scan = config["scan"]
    fixed = config["fixed_context"]
    half_width = int(fixed["central_block_half_width_pixels"])
    block_side = 2 * half_width + 1
    minimum_snr = float(config["screening_guide"]["minimum_snr"])
    conditions = ", ".join(
        f"{chr(ord('A') + index)} (F={float(value):g} mW us)"
        for index, value in enumerate(scan["selected_fluence_mw_us"])
    )
    return (
        "Initial-frame Faraday SNR as a function of F=P tau at "
        f"|Delta|/2pi={float(scan['detuning_ghz']):g} GHz for the undepleted "
        "reference condensate. The estimator uses the absolute expected signal in "
        f"the central {block_side}x{block_side} camera-pixel block and includes "
        f"Poisson noise plus {float(fixed['read_noise_electrons_per_pixel_per_port']):g} "
        f"e- rms read noise per pixel and port. The dashed line marks the working "
        f"SNR={minimum_snr:g} camera-space screening threshold. The labelled conditions "
        f"are {conditions}. "
        "The lower panels are fixed noise realisations and do not determine the curves."
    )


def _fluence_axis(scan: dict[str, Any]) -> np.ndarray:
    minimum = float(scan["fluence_min_mw_us"])
    maximum = float(scan["fluence_max_mw_us"])
    points = int(scan["fluence_points"])
    if scan["spacing"] == "linear":
        values = np.linspace(minimum, maximum, points)
    elif scan["spacing"] == "log":
        values = np.geomspace(minimum, maximum, points)
    else:
        raise ValueError("Figure 5.1 spacing must be 'linear' or 'log'")
    required = [
        float(scan["reference_fluence_mw_us"]),
        *[float(value) for value in scan["selected_fluence_mw_us"]],
    ]
    return np.unique(np.sort(np.concatenate([values, np.asarray(required, dtype=float)])))


def _central_block(array: np.ndarray, half_width: int) -> np.ndarray:
    centre = tuple(size // 2 for size in array.shape)
    return array[
        centre[0] - half_width : centre[0] + half_width + 1,
        centre[1] - half_width : centre[1] + half_width + 1,
    ]


def _single_output_snr(
    image: np.ndarray,
    background: float,
    photoelectrons_per_i0_pixel: np.ndarray,
    read_noise_e: float,
    half_width: int,
) -> np.ndarray:
    block = _central_block(image, half_width)
    signal_per_i0 = float(np.sum(block - background))
    photon_variance_per_i0 = float(np.sum(np.clip(block, 0.0, None)))
    counts = np.asarray(photoelectrons_per_i0_pixel, dtype=float)
    signal_e = signal_per_i0 * counts
    variance_e2 = photon_variance_per_i0 * counts + block.size * read_noise_e**2
    return np.abs(signal_e) / np.sqrt(variance_e2)


def _dual_port_snr(
    port_h: np.ndarray,
    port_v: np.ndarray,
    photoelectrons_per_i0_pixel: np.ndarray,
    read_noise_e: float,
    half_width: int,
) -> np.ndarray:
    block_h = _central_block(port_h, half_width)
    block_v = _central_block(port_v, half_width)
    signal_per_i0 = float(np.sum(block_h - block_v))
    photon_variance_per_i0 = float(np.sum(np.clip(block_h + block_v, 0.0, None)))
    counts = np.asarray(photoelectrons_per_i0_pixel, dtype=float)
    signal_e = signal_per_i0 * counts
    variance_e2 = photon_variance_per_i0 * counts + block_h.size * 2 * read_noise_e**2
    return np.abs(signal_e) / np.sqrt(variance_e2)


def _reference_fields(
    model: dict[str, Any],
    config: dict[str, Any],
) -> dict[str, Any]:
    fixed = config["fixed_context"]
    constants = _basic_constants(model)
    phase_stage = build_phase_stage(model)
    reference_phase = np.asarray(phase_stage["notebook_phase_map_rad"], dtype=float)
    reference_peak = float(phase_stage["notebook_phase_peak_rad"])
    phase_profile = reference_phase / reference_peak
    pupil_stage = build_pupil(model)
    pupil = np.asarray(pupil_stage["pupil"], dtype=float)

    if int(phase_stage["condensate_stage"]["imaging_axis"]) != int(fixed["imaging_axis"]):
        raise RuntimeError("Figure 5.1 config does not match the reference condensate imaging axis")
    if not np.isclose(float(pupil_stage["numerical_aperture"]), float(fixed["numerical_aperture"])):
        raise RuntimeError("Figure 5.1 config does not match the notebook numerical aperture")

    detuning_ghz = float(config["scan"]["detuning_ghz"])
    phase_peak = scalar_phase_shift(
        detuning_ghz * 1e9,
        float(constants["column_density"][int(fixed["imaging_axis"])]),
        float(constants["resonant_cross_section"]),
        float(constants["gamma"]),
    )
    phase_map = phase_peak * phase_profile
    kappa_f = float(model["faraday_recovery"]["kappa_F"])
    if not np.isclose(kappa_f, float(fixed["kappa_F"]), rtol=0, atol=1e-15):
        raise RuntimeError("Figure 5.1 kappa_F does not match the active model")
    faraday = simulate_faraday_image(kappa_f * phase_map, pupil)
    requested_object_pixel_m = (
        float(fixed["camera_pixel_m"]) / float(fixed["effective_magnification"])
    )
    computational_spacing_m = (
        float(model["grid"]["field_of_view_m"]) / int(model["grid"]["ngrid"])
    )
    output_shape = tuple(int(value) for value in fixed["camera_output_shape"])
    expected_output_shape = centered_camera_shape(
        reference_phase.shape,
        computational_spacing_m,
        requested_object_pixel_m,
    )
    if output_shape != expected_output_shape:
        raise RuntimeError(
            "Figure 5.1 camera output shape is inconsistent with the numerical field, "
            "physical pixel pitch and effective magnification"
        )
    dark = resample_to_camera_pixels(
        faraday["dark_field_intensity"],
        computational_spacing_m,
        requested_object_pixel_m,
        output_shape,
    )
    port_h = resample_to_camera_pixels(
        faraday["analyser_h_intensity"],
        computational_spacing_m,
        requested_object_pixel_m,
        output_shape,
    )
    port_v = resample_to_camera_pixels(
        faraday["analyser_v_intensity"],
        computational_spacing_m,
        requested_object_pixel_m,
        output_shape,
    )

    field_of_view_m = float(model["grid"]["field_of_view_m"])
    ngrid = int(model["grid"]["ngrid"])
    object_pixel_um = requested_object_pixel_m * 1e6
    half_extent_um = object_pixel_um * dark.shape[0] / 2
    return {
        "phase_peak_rad": float(phase_peak),
        "dark": dark,
        "port_h": port_h,
        "port_v": port_v,
        "physical_object_pixel_m": float(requested_object_pixel_m),
        "object_pixel_um": float(object_pixel_um),
        "extent_um": [-half_extent_um, half_extent_um, -half_extent_um, half_extent_um],
    }


def _noisy_images(
    fields: dict[str, Any],
    selected_fluence: np.ndarray,
    selected_photoelectrons: np.ndarray,
    read_noise: float,
    base_seed: int,
) -> dict[str, list[np.ndarray]]:
    images: dict[str, list[np.ndarray]] = {mode: [] for mode in MODES}
    for index, (_, count_scale) in enumerate(zip(selected_fluence, selected_photoelectrons, strict=True)):
        dark_rng = np.random.default_rng(np.random.SeedSequence([base_seed, 0, index]))
        dark_counts = simulate_noisy_camera_image(
            fields["dark"],
            float(count_scale),
            dark_rng,
            read_noise,
            input_is_binned=True,
            normalize=False,
        )
        images["Faraday dark-field"].append(np.asarray(dark_counts, dtype=float) / count_scale)

        dual_rng = np.random.default_rng(np.random.SeedSequence([base_seed, 1, index]))
        h_counts = simulate_noisy_camera_image(
            fields["port_h"],
            float(count_scale),
            dual_rng,
            read_noise,
            input_is_binned=True,
            normalize=False,
        )
        v_counts = simulate_noisy_camera_image(
            fields["port_v"],
            float(count_scale),
            dual_rng,
            read_noise,
            input_is_binned=True,
            normalize=False,
        )
        h_counts = np.asarray(h_counts, dtype=float)
        v_counts = np.asarray(v_counts, dtype=float)
        total_counts = h_counts + v_counts
        images["Faraday dual-port"].append(
            np.divide(
                h_counts - v_counts,
                total_counts,
                out=np.zeros_like(total_counts),
                where=np.abs(total_counts) > np.finfo(float).eps,
            )
        )
    return images


def build_scan(config: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    model = load_config(REPO_ROOT / config["model_config"])
    fixed = config["fixed_context"]
    scan = config["scan"]
    fields = _reference_fields(model, config)
    fluence = _fluence_axis(scan)
    constants = _basic_constants(model)
    h_planck = 2 * np.pi * float(model["constants"]["hbar"])
    photon_energy = (
        h_planck
        * float(model["constants"]["speed_of_light"])
        / float(constants["wavelength"])
    )
    incident_intensity_for_1mw = intensity_at_atoms_notebook(model, 1.0)
    photoelectrons = (
        incident_intensity_for_1mw
        * float(fields["physical_object_pixel_m"]) ** 2
        * fluence
        * 1e-6
        * float(fixed["quantum_efficiency"])
        / photon_energy
    )
    half_width = int(fixed["central_block_half_width_pixels"])
    read_noise = float(model["camera_recovery"]["read_noise_electrons"])
    if not np.isclose(
        read_noise,
        float(fixed["read_noise_electrons_per_pixel_per_port"]),
        rtol=0,
        atol=1e-15,
    ):
        raise RuntimeError("Figure 5.1 read noise does not match the active model")
    snr = {
        "Faraday dark-field": _single_output_snr(
            fields["dark"],
            0.0,
            photoelectrons,
            read_noise,
            half_width,
        ),
        "Faraday dual-port": _dual_port_snr(
            fields["port_h"],
            fields["port_v"],
            photoelectrons,
            read_noise,
            half_width,
        ),
    }

    rows: list[dict[str, Any]] = []
    for mode in MODES:
        for index, fluence_value in enumerate(fluence):
            rows.append(
                {
                    "detuning_ghz": f"{float(scan['detuning_ghz']):.12g}",
                    "fluence_mw_us": f"{fluence_value:.12g}",
                    "mode": mode,
                    "central_block_snr": f"{snr[mode][index]:.12g}",
                    "phase_peak_rad": f"{fields['phase_peak_rad']:.12g}",
                    "photoelectrons_per_incident_I0_pixel": f"{photoelectrons[index]:.12g}",
                    "read_noise_e_per_pixel_per_port": f"{read_noise:.12g}",
                    "kappa_F": f"{float(fixed['kappa_F']):.12g}",
                }
            )

    selected_fluence = np.asarray(scan["selected_fluence_mw_us"], dtype=float)
    selected_indices = np.asarray(
        [
            int(np.flatnonzero(np.isclose(fluence, value, rtol=0, atol=1e-12))[0])
            for value in selected_fluence
        ],
        dtype=int,
    )
    selected_photoelectrons = photoelectrons[selected_indices]
    selected_snr = {
        mode: [float(snr[mode][index]) for index in selected_indices]
        for mode in MODES
    }
    noisy_images = _noisy_images(
        fields,
        selected_fluence,
        selected_photoelectrons,
        read_noise,
        int(config["noise_realisation"]["rng_seed"]),
    )

    reference_fluence = float(scan["reference_fluence_mw_us"])
    reference_index = int(
        np.flatnonzero(np.isclose(fluence, reference_fluence, rtol=0, atol=1e-12))[0]
    )
    canonical = {
        "detuning_ghz": float(scan["detuning_ghz"]),
        "fluence_mw_us": float(fluence[reference_index]),
        "phase_peak_rad": float(fields["phase_peak_rad"]),
        "photoelectrons_per_incident_I0_pixel": float(photoelectrons[reference_index]),
        "snr": {mode: float(snr[mode][reference_index]) for mode in MODES},
    }
    checks = config["canonical_checks"]
    for name in ("phase_peak_rad", "photoelectrons_per_incident_I0_pixel"):
        if not np.isclose(
            canonical[name],
            float(checks[name]),
            rtol=float(checks["relative_tolerance"]),
            atol=float(checks["absolute_tolerance"]),
        ):
            raise RuntimeError(f"Figure 5.1 canonical check failed for {name}")
    for mode, expected in checks["snr"].items():
        if not np.isclose(
            canonical["snr"][mode],
            float(expected),
            rtol=float(checks["relative_tolerance"]),
            atol=float(checks["absolute_tolerance"]),
        ):
            raise RuntimeError(f"Figure 5.1 canonical SNR check failed for {mode}")

    summary = {
        "canonical_gate": {"passed": True, "reference": canonical},
        "scan": {
            "detuning_ghz": float(scan["detuning_ghz"]),
            "fluence_mw_us": [float(value) for value in fluence],
            "snr_ranges": {
                mode: {"min": float(np.min(values)), "max": float(np.max(values))}
                for mode, values in snr.items()
            },
        },
        "selected_frames": {
            "fluence_mw_us": [float(value) for value in selected_fluence],
            "photoelectrons_per_incident_I0_pixel": [
                float(value) for value in selected_photoelectrons
            ],
            "snr": selected_snr,
            "rng_seed": int(config["noise_realisation"]["rng_seed"]),
            "interpretation": "fixed representative draws for camera appearance; not used to estimate SNR",
        },
        "estimator": {
            "name": (
                f"central {2 * half_width + 1}x{2 * half_width + 1} "
                "binned-camera-pixel block"
            ),
            "single_output": "abs(sum(signal_e))/sqrt(sum(total_count_e + read_variance_e2))",
            "dual_port": "abs(sum(N_H-N_V))/sqrt(sum(N_H+N_V+2*sigma_r^2))",
        },
        "screening_guide": config["screening_guide"],
        "fixed_context": fixed,
        "scope": config["scope"],
    }
    data = {
        "fluence_mw_us": fluence,
        "photoelectrons": photoelectrons,
        "snr": snr,
        "selected_fluence_mw_us": selected_fluence,
        "selected_snr": selected_snr,
        "noisy_images": noisy_images,
        "extent_um": fields["extent_um"],
        "rows": rows,
    }
    return data, summary


def _plot(data: dict[str, Any], config: dict[str, Any], paths: dict[str, Path]) -> None:
    figure = config["figure"]
    display = figure["image_display"]
    fluence = np.asarray(data["fluence_mw_us"], dtype=float)
    selected = np.asarray(data["selected_fluence_mw_us"], dtype=float)

    plt.style.use("default")
    matplotlib.rcParams["svg.hashsalt"] = "figure-5-1-faraday-snr-and-camera-frames"
    matplotlib.rcParams.update(
        {
            "font.size": 8.8,
            "axes.labelsize": 9.4,
            "axes.titlesize": 9.4,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    fig = plt.figure(figsize=tuple(figure["figsize_inches"]), constrained_layout=True)
    grid = fig.add_gridspec(3, 3, height_ratios=[1.35, 1.0, 1.0])
    curve_axis = fig.add_subplot(grid[0, :])
    image_axes = [
        [fig.add_subplot(grid[row, column]) for column in range(3)]
        for row in (1, 2)
    ]

    for mode in MODES:
        curve_axis.plot(
            fluence,
            np.asarray(data["snr"][mode], dtype=float),
            color=MODE_COLOURS[mode],
            linewidth=1.8,
            label=mode.replace("Faraday ", "").capitalize(),
        )
    markers = ("o", "s", "^")
    column_labels = ("A", "B", "C")
    for column, (label, marker, value) in enumerate(zip(column_labels, markers, selected, strict=True)):
        curve_axis.axvline(value, color="0.75", linewidth=0.8, linestyle=":")
        for mode in MODES:
            snr_value = float(data["selected_snr"][mode][column])
            curve_axis.plot(
                value,
                snr_value,
                marker=marker,
                markersize=5.2,
                markerfacecolor="white",
                markeredgecolor=MODE_COLOURS[mode],
                markeredgewidth=1.2,
                linestyle="none",
            )
        curve_axis.text(
            value,
            0.98,
            label,
            transform=curve_axis.get_xaxis_transform(),
            ha="center",
            va="top",
            fontweight="bold",
        )
    curve_axis.set_xlim(float(figure["xlim_mw_us"][0]), float(figure["xlim_mw_us"][1]))
    curve_axis.set_ylim(float(figure["ylim_snr"][0]), float(figure["ylim_snr"][1]))
    curve_axis.set_xlabel(r"Fluence coordinate $F=P\tau$ (mW $\mu$s)")
    half_width = int(config["fixed_context"]["central_block_half_width_pixels"])
    block_side = 2 * half_width + 1
    curve_axis.set_ylabel(
        rf"Initial-frame SNR (central ${block_side}\times{block_side}$ block)"
    )
    curve_axis.grid(color="0.88", linewidth=0.6)
    minimum_snr = float(config["screening_guide"]["minimum_snr"])
    curve_axis.axhline(
        minimum_snr,
        color="0.35",
        linewidth=1.0,
        linestyle="--",
        label=rf"screening threshold, SNR $={minimum_snr:g}$",
    )
    curve_axis.legend(frameon=False, loc="lower right")
    curve_axis.set_title(
        rf"(a) Fixed detuning $|\Delta|/2\pi={float(config['scan']['detuning_ghz']):g}$ GHz",
        loc="left",
        pad=5,
    )

    row_modes = ("Faraday dark-field", "Faraday dual-port")
    norms = {
        "Faraday dark-field": TwoSlopeNorm(
            vmin=float(display["dark_vmin"]),
            vcenter=0.0,
            vmax=float(display["dark_vmax"]),
        ),
        "Faraday dual-port": TwoSlopeNorm(
            vmin=float(display["dual_vmin"]),
            vcenter=0.0,
            vmax=float(display["dual_vmax"]),
        ),
    }
    row_labels = {
        "Faraday dark-field": "Dark-field",
        "Faraday dual-port": "Dual-port",
    }
    colourbar_labels = {
        "Faraday dark-field": r"$I_{\rm DF}/I_0$",
        "Faraday dual-port": r"$S$",
    }
    extent = data["extent_um"]
    for row, mode in enumerate(row_modes):
        image_handle = None
        for column, axis in enumerate(image_axes[row]):
            image_handle = axis.imshow(
                data["noisy_images"][mode][column],
                origin="lower",
                extent=extent,
                cmap="coolwarm",
                norm=norms[mode],
                interpolation="nearest",
                aspect="equal",
            )
            axis.set_xlim(*display["xlim_um"])
            axis.set_ylim(*display["ylim_um"])
            axis.set_title(
                (
                    f"{column_labels[column]}: "
                    rf"$F={selected[column]:g}$ mW $\mu$s"
                    "\n"
                    rf"SNR $={float(data['selected_snr'][mode][column]):.2f}$"
                ),
                pad=4,
            )
            if row == 1:
                axis.set_xlabel(r"$y$ ($\mu$m)")
            else:
                axis.tick_params(labelbottom=False)
            if column == 0:
                axis.set_ylabel(r"$z$ ($\mu$m)")
                axis.text(
                    -0.32,
                    0.5,
                    row_labels[mode],
                    transform=axis.transAxes,
                    ha="center",
                    va="center",
                    rotation=90,
                )
            else:
                axis.tick_params(labelleft=False)
        if image_handle is None:
            raise RuntimeError("Figure 5.1 produced no image panels")
        colourbar = fig.colorbar(
            image_handle,
            ax=image_axes[row],
            location="right",
            fraction=0.025,
            pad=0.02,
        )
        colourbar.set_label(colourbar_labels[mode])

    fig.savefig(paths["svg"], bbox_inches="tight")
    fig.savefig(paths["png"], dpi=int(figure["dpi"]), bbox_inches="tight")
    fig.savefig(paths["pdf"], bbox_inches="tight")
    plt.close(fig)


def generate(config_path: Path = DEFAULT_CONFIG, output_dir: Path | None = None) -> dict[str, Path]:
    config_path = config_path if config_path.is_absolute() else REPO_ROOT / config_path
    config = _load_json(config_path)
    destination = output_dir or REPO_ROOT / config["output_directory"]
    destination.mkdir(parents=True, exist_ok=True)
    paths = {
        "svg": destination / config["figure"]["svg_filename"],
        "png": destination / config["figure"]["png_filename"],
        "pdf": destination / config["figure"]["pdf_filename"],
        "csv": destination / "figure_5_1_data.csv",
        "values": destination / "figure_5_1_values.json",
        "metadata": destination / "metadata.json",
    }
    data, summary = build_scan(config)
    _plot(data, config, paths)
    _write_csv(paths["csv"], data["rows"])
    _write_json(paths["values"], summary)
    _write_json(
        paths["metadata"],
        {
            "label": config["label"],
            "config_path": str(config_path.relative_to(REPO_ROOT)),
            "git_branch": _git_value("branch", "--show-current"),
            "git_commit": _git_value("rev-parse", "HEAD"),
            "outputs": {
                key: str(path.relative_to(REPO_ROOT)) if path.is_relative_to(REPO_ROOT) else str(path)
                for key, path in paths.items()
            },
            "estimator": summary["estimator"],
            "screening_guide": summary["screening_guide"],
            "fixed_context": config["fixed_context"],
            "selected_frames": summary["selected_frames"],
            "scope": config["scope"],
            "faraday_boundary": config["fixed_context"]["faraday_calibration_status"],
            "noise_draw_boundary": config["noise_realisation"]["interpretation"],
            "caption": _caption(config),
        },
    )
    return paths


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    outputs = generate(args.config, args.output_dir)
    for label, path in outputs.items():
        print(f"- {label}: {path}")


if __name__ == "__main__":
    main()
