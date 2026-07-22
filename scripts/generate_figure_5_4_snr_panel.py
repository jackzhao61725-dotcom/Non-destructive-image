"""Generate Figure 5.4: fixed-detuning noisy dual-port multiframe sequences."""

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
from matplotlib.colors import TwoSlopeNorm, to_rgb
from matplotlib.patches import Patch
import numpy as np
from scipy.special import zeta

from non_destructive_image import (
    resample_to_camera_pixels,
    simulate_faraday_image,
    simulate_noisy_camera_image,
)
from scripts.generate_figure_5_2 import (
    REPO_ROOT,
    _basic_constants,
    _central_block,
    _load_json,
    _photoelectron_scale,
    load_config,
    reabsorption_fraction,
    scalar_phase_peak,
    scattered_photons_per_atom,
    self_consistent_total_atoms,
    tf_state_for_atoms,
)
from scripts.recover_notebook_pci_stage import build_pupil


DEFAULT_CONFIG = REPO_ROOT / "configs" / "figure_5_4.json"


def _json_value(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: _json_value(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_value(item) for item in value]
    return value


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(_json_value(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


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


def _grid_context(model: dict[str, Any], scan_config: dict[str, Any]) -> dict[str, Any]:
    ngrid = int(model["grid"]["ngrid"])
    field_of_view_m = float(model["grid"]["field_of_view_m"])
    spacing_m = field_of_view_m / ngrid
    coordinate = (np.arange(ngrid) - ngrid // 2) * spacing_m
    coordinate_a, coordinate_b = np.meshgrid(coordinate, coordinate)
    imaging_axis = int(scan_config["scan"]["imaging_axis"])
    plane = [index for index in range(3) if index != imaging_axis]
    camera_shape = tuple(
        int(value) for value in model["camera_recovery"]["camera_output_shape"]
    )
    object_pixel_m = float(model["camera_recovery"]["object_plane_pixel_m"])
    object_pixel_um = object_pixel_m * 1e6
    half_extent_um = object_pixel_um * camera_shape[0] / 2
    camera_coordinate_um = (
        np.arange(camera_shape[0]) - (camera_shape[0] - 1) / 2
    ) * object_pixel_um
    return {
        "pupil": np.asarray(build_pupil(model)["pupil"], dtype=float),
        "coordinate_a": coordinate_a,
        "coordinate_b": coordinate_b,
        "imaging_axis": imaging_axis,
        "plane": plane,
        "input_spacing_m": float(spacing_m),
        "object_pixel_m": float(object_pixel_m),
        "camera_shape": camera_shape,
        "object_pixel_um": float(object_pixel_um),
        "camera_coordinate_um": camera_coordinate_um,
        "extent_um": [
            -half_extent_um,
            half_extent_um,
            -half_extent_um,
            half_extent_um,
        ],
    }


def _dual_port_fields(
    condensate_atoms: float,
    detuning_hz: float,
    kappa_f: float,
    constants: dict[str, Any],
    context: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    if condensate_atoms > 0:
        state = tf_state_for_atoms(condensate_atoms, constants)
        radii = np.asarray(state["radii"], dtype=float)
        plane = context["plane"]
        profile = np.maximum(
            0.0,
            1.0
            - context["coordinate_a"] ** 2 / radii[plane[0]] ** 2
            - context["coordinate_b"] ** 2 / radii[plane[1]] ** 2,
        ) ** 1.5
        phase_peak = scalar_phase_peak(
            detuning_hz,
            float(np.asarray(state["column_density"])[context["imaging_axis"]]),
            constants,
        )
        theta_map = kappa_f * phase_peak * profile
    else:
        theta_map = np.zeros_like(context["coordinate_a"], dtype=float)

    faraday = simulate_faraday_image(theta_map, context["pupil"])
    port_h = resample_to_camera_pixels(
        faraday["analyser_h_intensity"],
        context["input_spacing_m"],
        context["object_pixel_m"],
        context["camera_shape"],
    )
    port_v = resample_to_camera_pixels(
        faraday["analyser_v_intensity"],
        context["input_spacing_m"],
        context["object_pixel_m"],
        context["camera_shape"],
    )
    return np.asarray(port_h, dtype=float), np.asarray(port_v, dtype=float)


def _dual_port_snr(
    port_h: np.ndarray,
    port_v: np.ndarray,
    count_scale: float,
    read_noise_e: float,
    half_width: int,
) -> float:
    block_h = _central_block(port_h, half_width)
    block_v = _central_block(port_v, half_width)
    signal_e = float(np.sum(block_h - block_v)) * count_scale
    variance_e2 = (
        float(np.sum(np.clip(block_h + block_v, 0.0, None))) * count_scale
        + block_h.size * 2 * read_noise_e**2
    )
    return float(abs(signal_e) / np.sqrt(variance_e2))


def _noisy_dual_port_signal(
    port_h: np.ndarray,
    port_v: np.ndarray,
    count_scale: float,
    read_noise_e: float,
    seed: np.random.SeedSequence,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    h_counts = np.asarray(
        simulate_noisy_camera_image(
            port_h,
            count_scale,
            rng,
            read_noise_e,
            input_is_binned=True,
            normalize=False,
        ),
        dtype=float,
    )
    v_counts = np.asarray(
        simulate_noisy_camera_image(
            port_v,
            count_scale,
            rng,
            read_noise_e,
            input_is_binned=True,
            normalize=False,
        ),
        dtype=float,
    )
    total = h_counts + v_counts
    return np.divide(
        h_counts - v_counts,
        total,
        out=np.zeros_like(total),
        where=np.abs(total) > np.finfo(float).eps,
    )


def _crop_and_rotate(
    image: np.ndarray,
    context: dict[str, Any],
    display: dict[str, Any],
) -> np.ndarray:
    coordinate = np.asarray(context["camera_coordinate_um"], dtype=float)
    y_mask = (coordinate >= float(display["xlim_um"][0])) & (
        coordinate <= float(display["xlim_um"][1])
    )
    z_mask = (coordinate >= float(display["ylim_um"][0])) & (
        coordinate <= float(display["ylim_um"][1])
    )
    cropped = image[np.ix_(z_mask, y_mask)]
    rotations = int(display["rotate_degrees_counterclockwise"]) // 90
    return np.rot90(cropped, k=rotations)


def _sequence_energy_context(
    model: dict[str, Any],
    constants: dict[str, Any],
    thermal: dict[str, float],
    detuning_hz: float,
    fluence_mw_us: float,
    power_mw: float,
) -> dict[str, float]:
    exposure_s = fluence_mw_us / power_mw * 1e-6
    n_gamma = scattered_photons_per_atom(
        model,
        constants,
        detuning_hz,
        power_mw,
        exposure_s,
    )
    reabsorption = reabsorption_fraction(detuning_hz, constants)
    recoil_multiplier = float(model["multishot_recovery"]["recoil_energy_multiplier"])
    deposited_energy = (
        recoil_multiplier * n_gamma * (1 + reabsorption) * constants["e_rec"]
    )
    critical_temperature = float(thermal["critical_temperature_k"])
    energy_coefficient = (
        3
        * float(zeta(4) / zeta(3))
        * float(constants["boltzmann_constant"])
        / critical_temperature**3
    )
    return {
        "n_gamma": float(n_gamma),
        "reabsorption_fraction": float(reabsorption),
        "deposited_energy_j_per_atom_per_frame": float(deposited_energy),
        "critical_temperature_k": critical_temperature,
        "energy_coefficient_j_per_k4": float(energy_coefficient),
    }


def _next_state(
    temperature: float,
    energy: dict[str, float],
    total_atoms: float,
) -> tuple[float, float]:
    next_temperature = (
        temperature**4
        + energy["deposited_energy_j_per_atom_per_frame"]
        / energy["energy_coefficient_j_per_k4"]
    ) ** 0.25
    raw_condensate_atoms = total_atoms * (
        1 - (next_temperature / energy["critical_temperature_k"]) ** 3
    )
    return float(next_temperature), float(max(0.0, raw_condensate_atoms))


def _depletion_limited_count(
    initial_temperature: float,
    initial_condensate_atoms: float,
    total_atoms: float,
    energy: dict[str, float],
    loss_limit: float,
    maximum_frames: int,
) -> int:
    temperature = initial_temperature
    accepted = 0
    for _ in range(maximum_frames):
        next_temperature, next_atoms = _next_state(temperature, energy, total_atoms)
        next_loss = 1 - next_atoms / initial_condensate_atoms
        if next_loss > loss_limit + 1e-12:
            break
        accepted += 1
        temperature = next_temperature
    return accepted


def build_sequences(config: dict[str, Any]) -> dict[str, Any]:
    figure_5_2_config = _load_json(REPO_ROOT / config["figure_5_2_config"])
    model = load_config(REPO_ROOT / figure_5_2_config["model_config"])
    constants = _basic_constants(model)
    thermal = self_consistent_total_atoms(model, constants)
    context = _grid_context(model, figure_5_2_config)
    sequence_config = config["sequence"]
    frames_to_show = int(sequence_config["displayed_frames"])
    detuning_hz = float(config["detuning_ghz"]) * 1e9
    power_mw = float(figure_5_2_config["scan"]["probe_power_mw"])
    kappa_f = float(model["faraday_recovery"]["kappa_F"])
    read_noise = float(model["camera_recovery"]["read_noise_electrons"])
    half_width = int(figure_5_2_config["frame_criteria"]["central_block_half_width_pixels"])
    initial_temperature = float(model["condensate"]["temperature_k"])
    initial_condensate_atoms = float(constants["atom_number"])
    total_atoms = float(thermal["total_atoms"])
    minimum_snr = float(sequence_config["minimum_snr"])
    loss_limit = float(sequence_config["maximum_post_exposure_condensate_loss"])
    base_seed = int(sequence_config["rng_seed"])

    sequences: dict[str, Any] = {}
    for row_index, fluence in enumerate(config["selected_fluence_mw_us"]):
        fluence_value = float(fluence)
        energy = _sequence_energy_context(
            model,
            constants,
            thermal,
            detuning_hz,
            fluence_value,
            power_mw,
        )
        count_scale = _photoelectron_scale(model, constants, fluence_value)
        temperature = initial_temperature
        frames: list[dict[str, Any]] = []
        quality_open = True
        usable_count = 0

        for frame_index in range(1, frames_to_show + 1):
            raw_condensate_atoms = total_atoms * (
                1 - (temperature / energy["critical_temperature_k"]) ** 3
            )
            condensate_atoms = float(max(0.0, raw_condensate_atoms))
            pre_loss = float(1 - condensate_atoms / initial_condensate_atoms)
            port_h, port_v = _dual_port_fields(
                condensate_atoms,
                detuning_hz,
                kappa_f,
                constants,
                context,
            )
            snr = _dual_port_snr(
                port_h,
                port_v,
                count_scale,
                read_noise,
                half_width,
            )
            noisy_signal = _noisy_dual_port_signal(
                port_h,
                port_v,
                count_scale,
                read_noise,
                np.random.SeedSequence([base_seed, row_index, frame_index]),
            )
            display_image = _crop_and_rotate(noisy_signal, context, config["image_display"])
            next_temperature, next_condensate_atoms = _next_state(
                temperature,
                energy,
                total_atoms,
            )
            post_loss = float(1 - next_condensate_atoms / initial_condensate_atoms)

            quality_failed = snr < minimum_snr
            depletion_failed = post_loss > loss_limit + 1e-12
            if quality_failed and depletion_failed:
                status = "both_limits"
            elif depletion_failed:
                status = "depletion_limited"
            elif quality_failed:
                status = "quality_limited"
            else:
                status = "usable"

            accepted = status == "usable" and quality_open
            if accepted:
                usable_count += 1
            else:
                quality_open = False

            frames.append(
                {
                    "image_number": frame_index,
                    "temperature_k": float(temperature),
                    "condensate_atoms": condensate_atoms,
                    "condensate_fraction": condensate_atoms / initial_condensate_atoms,
                    "pre_exposure_loss": pre_loss,
                    "post_exposure_loss": post_loss,
                    "snr_5x5": snr,
                    "status": status,
                    "counted_as_usable": accepted,
                    "camera_signal": display_image,
                }
            )
            temperature = next_temperature

        depletion_count = _depletion_limited_count(
            initial_temperature,
            initial_condensate_atoms,
            total_atoms,
            energy,
            loss_limit,
            int(model["multishot_recovery"]["max_shots"]),
        )
        sequences[f"{fluence_value:g}"] = {
            "fluence_mw_us": fluence_value,
            "row_label": config["row_labels"][f"{fluence_value:g}"],
            "photoelectrons_per_incident_i0_pixel": count_scale,
            "energy": energy,
            "depletion_limited_frames": depletion_count,
            "usable_frames": usable_count,
            "frames": frames,
        }

    expected_counts = config["canonical_checks"]["sequence_counts"]
    for key, expected in expected_counts.items():
        sequence = sequences[key]
        if sequence["depletion_limited_frames"] != expected["depletion"]:
            raise RuntimeError(f"Figure 5.4 depletion count drifted at F={key}")
        if sequence["usable_frames"] != expected["usable"]:
            raise RuntimeError(f"Figure 5.4 usable-frame count drifted at F={key}")

    return {
        "model": model,
        "figure_5_2_config": figure_5_2_config,
        "context": context,
        "sequences": sequences,
    }


def _rgb_strip(
    frames: list[dict[str, Any]],
    display: dict[str, Any],
) -> tuple[np.ndarray, list[float]]:
    norm = TwoSlopeNorm(
        vmin=float(display["vmin"]),
        vcenter=0.0,
        vmax=float(display["vmax"]),
    )
    cmap = plt.get_cmap(display["cmap"])
    band_height = int(display["status_band_pixels"])
    separator_width = int(display["separator_pixels"])
    status_colours = {
        key: np.asarray(to_rgb(value), dtype=float)
        for key, value in display["status_colours"].items()
        if key != "both_limits"
    }
    combined_colours = [
        np.asarray(to_rgb(value), dtype=float)
        for value in display["status_colours"]["both_limits"]
    ]
    strips: list[np.ndarray] = []
    centres: list[float] = []
    cursor = 0
    for frame in frames:
        image = np.asarray(frame["camera_signal"], dtype=float)
        rgb = np.asarray(cmap(norm(image))[..., :3], dtype=float)
        band = np.ones((band_height, rgb.shape[1], 3), dtype=float)
        if frame["status"] == "both_limits":
            yy, xx = np.indices((band_height, rgb.shape[1]))
            diagonal = ((xx + 2 * yy) // 3) % 2
            band[diagonal == 0] = combined_colours[0]
            band[diagonal == 1] = combined_colours[1]
        elif frame["status"] != "usable":
            band[:] = status_colours[frame["status"]]
        tile = np.concatenate([band, rgb], axis=0)
        strips.append(tile)
        centres.append(cursor + tile.shape[1] / 2)
        cursor += tile.shape[1]
        if separator_width > 0 and frame is not frames[-1]:
            strips.append(np.ones((tile.shape[0], separator_width, 3), dtype=float))
            cursor += separator_width
    return np.concatenate(strips, axis=1), centres


def _plot(
    config: dict[str, Any],
    result: dict[str, Any],
    destination: Path,
) -> dict[str, Path]:
    figure_config = config["figure"]
    display = config["image_display"]
    sequence_config = config["sequence"]
    selected = [float(value) for value in config["selected_fluence_mw_us"]]

    plt.style.use("default")
    matplotlib.rcParams["svg.hashsalt"] = "figure-5-4-dual-port-multiframe"
    matplotlib.rcParams.update(
        {
            "font.size": 8.4,
            "axes.labelsize": 9.2,
            "axes.titlesize": 9.2,
            "xtick.labelsize": 7.6,
            "ytick.labelsize": 7.6,
            "axes.linewidth": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    figure, axes = plt.subplots(
        len(selected),
        1,
        figsize=tuple(figure_config["figsize_inches"]),
        sharex=True,
    )
    figure.subplots_adjust(left=0.185, right=0.88, top=0.86, bottom=0.12, hspace=0.13)

    tick_numbers = [int(value) for value in figure_config["labelled_frame_numbers"]]
    strip_handles: list[np.ndarray] = []
    tick_positions: list[float] | None = None
    for axis, fluence in zip(axes, selected, strict=True):
        key = f"{fluence:g}"
        sequence = result["sequences"][key]
        strip, centres = _rgb_strip(sequence["frames"], display)
        strip_handles.append(strip)
        if tick_positions is None:
            tick_positions = [centres[number - 1] for number in tick_numbers]
        axis.imshow(strip, origin="upper", interpolation="nearest", aspect="auto")
        axis.set_yticks([])
        axis.tick_params(axis="x", length=2.5, pad=2)
        row_label = sequence["row_label"]
        axis.set_ylabel(
            (
                rf"{row_label}   $F={fluence:g}$"
                "\n"
                rf"$N_{{\rm use}}={sequence['usable_frames']}$, "
                rf"$N_{{\rm dep}}={sequence['depletion_limited_frames']}$"
            ),
            rotation=0,
            ha="right",
            va="center",
            labelpad=8,
        )
        for side in ("left", "right", "bottom"):
            axis.spines[side].set_visible(False)
        axis.spines["top"].set_color("0.75")
        axis.spines["top"].set_linewidth(0.5)

    if tick_positions is None:
        raise RuntimeError("Figure 5.4 produced no filmstrip rows")
    axes[0].set_xticks(tick_positions, [str(value) for value in tick_numbers])
    axes[0].xaxis.set_ticks_position("top")
    axes[0].xaxis.set_label_position("top")
    axes[0].set_xlabel("Image number", labelpad=5)
    for axis in axes[1:]:
        axis.tick_params(axis="x", labelbottom=False, bottom=False)

    normaliser = TwoSlopeNorm(
        vmin=float(display["vmin"]),
        vcenter=0.0,
        vmax=float(display["vmax"]),
    )
    scalar_mappable = matplotlib.cm.ScalarMappable(
        norm=normaliser,
        cmap=display["cmap"],
    )
    colourbar_axis = figure.add_axes([0.895, 0.22, 0.018, 0.52])
    colourbar = figure.colorbar(scalar_mappable, cax=colourbar_axis)
    colourbar.set_label(r"Dual-port signal $S$")
    colourbar.set_ticks(
        [float(display["vmin"]), 0.0, float(display["vmax"])]
    )

    legend = [
        Patch(
            facecolor=display["status_colours"]["quality_limited"],
            edgecolor="none",
            label=rf"SNR below {float(sequence_config['minimum_snr']):g} only",
        ),
        Patch(
            facecolor=display["status_colours"]["depletion_limited"],
            edgecolor="none",
            label="loss above 30% only",
        ),
        Patch(
            facecolor=display["status_colours"]["both_limits"][1],
            edgecolor=display["status_colours"]["both_limits"][0],
            hatch="///",
            label="both limits exceeded",
        ),
    ]
    figure.legend(
        handles=legend,
        loc="lower center",
        bbox_to_anchor=(0.53, 0.025),
        ncol=3,
        frameon=False,
        fontsize=7.7,
        handlelength=1.5,
        columnspacing=1.8,
    )

    destination.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, Path] = {}
    for extension in figure_config["formats"]:
        path = destination / f"{figure_config['output_stem']}.{extension}"
        kwargs: dict[str, Any] = {"bbox_inches": "tight"}
        if extension == "png":
            kwargs["dpi"] = int(figure_config["dpi"])
        figure.savefig(path, **kwargs)
        outputs[extension] = path
    plt.close(figure)
    return outputs


def _write_frame_csv(path: Path, result: dict[str, Any]) -> None:
    rows: list[dict[str, Any]] = []
    for sequence in result["sequences"].values():
        for frame in sequence["frames"]:
            rows.append(
                {
                    "fluence_mw_us": f"{sequence['fluence_mw_us']:.12g}",
                    "image_number": frame["image_number"],
                    "snr_5x5": f"{frame['snr_5x5']:.12g}",
                    "condensate_atoms": f"{frame['condensate_atoms']:.12g}",
                    "condensate_fraction": f"{frame['condensate_fraction']:.12g}",
                    "temperature_k": f"{frame['temperature_k']:.12g}",
                    "pre_exposure_loss": f"{frame['pre_exposure_loss']:.12g}",
                    "post_exposure_loss": f"{frame['post_exposure_loss']:.12g}",
                    "status": frame["status"],
                    "counted_as_usable": frame["counted_as_usable"],
                }
            )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _output_destination(config: dict[str, Any], output_dir: Path | None) -> Path:
    canonical = (REPO_ROOT / config["output_directory"]).resolve()
    destination = output_dir or canonical
    if not destination.is_absolute():
        destination = REPO_ROOT / destination
    destination = destination.resolve()
    lifecycle = config.get("lifecycle", {})
    if (
        not lifecycle.get("canonical_regeneration_allowed", True)
        and destination == canonical
    ):
        raise RuntimeError(
            "canonical Figure 5.4 is frozen pending the approved heating "
            "replacement; provide --output-dir for a non-canonical historical "
            "reproduction"
        )
    return destination


def generate(
    config_path: Path = DEFAULT_CONFIG,
    output_dir: Path | None = None,
) -> dict[str, Path]:
    config_path = config_path if config_path.is_absolute() else REPO_ROOT / config_path
    config = _load_json(config_path)
    destination = _output_destination(config, output_dir)
    result = build_sequences(config)
    outputs = _plot(config, result, destination)
    csv_path = destination / "figure_5_4_frame_data.csv"
    _write_frame_csv(csv_path, result)
    outputs["csv"] = csv_path

    metadata_path = destination / "metadata.json"
    _write_json(
        metadata_path,
        {
            "label": config["label"],
            "config_path": str(config_path.relative_to(REPO_ROOT)),
            "git_branch": _git_value("branch", "--show-current"),
            "git_commit": _git_value("rev-parse", "HEAD"),
            "result_status": config["lifecycle"],
            "detuning_ghz": config["detuning_ghz"],
            "selected_fluence_mw_us": config["selected_fluence_mw_us"],
            "displayed_frames": config["sequence"]["displayed_frames"],
            "criteria": {
                "minimum_snr": config["sequence"]["minimum_snr"],
                "maximum_post_exposure_condensate_loss": config["sequence"][
                    "maximum_post_exposure_condensate_loss"
                ],
                "status_encoding": config["sequence"]["status_encoding"],
            },
            "sequence_counts": {
                key: {
                    "depletion_limited_frames": sequence[
                        "depletion_limited_frames"
                    ],
                    "usable_frames": sequence["usable_frames"],
                    "photoelectrons_per_incident_i0_pixel": sequence[
                        "photoelectrons_per_incident_i0_pixel"
                    ],
                }
                for key, sequence in result["sequences"].items()
            },
            "camera": {
                "model": result["model"]["camera_recovery"]["camera_model"],
                "magnification": result["model"]["imaging_geometry"]["magnification"],
                "quantum_efficiency": result["model"]["camera_recovery"][
                    "quantum_efficiency"
                ],
                "read_noise_electrons_per_pixel_per_port": result["model"][
                    "camera_recovery"
                ]["read_noise_electrons"],
                "sampling_method": result["model"]["camera_recovery"]["sampling_method"],
                "camera_output_shape": result["model"]["camera_recovery"]["camera_output_shape"],
                "object_plane_pixel_m": result["model"]["camera_recovery"]["object_plane_pixel_m"],
                "rng_seed": config["sequence"]["rng_seed"],
                "readout_mode": result["model"]["camera_recovery"]["readout_mode"],
            },
            "atomic_response": {
                "kappa_F": result["model"]["faraday_recovery"]["kappa_F"],
                "coefficient_status": result["model"]["faraday_recovery"][
                    "kappa_F_status"
                ],
            },
            "optics": {
                "effective_numerical_aperture": float(
                    build_pupil(result["model"])["numerical_aperture"]
                )
            },
            "display": config["image_display"],
            "scope": config["scope"],
            "caption": (
                f"{int(config['sequence']['displayed_frames'])}-frame DPFI sequences "
                "at |Delta|/2pi=1.5 GHz for "
                "A, B and C (F=90, 150 and 300 mW us). Each tile shows one fixed-noise "
                "realisation of the signed signal S=(I_H-I_V)/(I_H+I_V). Row labels "
                "report N_use and N_dep. Orange marks SNR failure, red marks depletion "
                "above 30%, and striped bands mark simultaneous failure. Frames beyond "
                "the first failed criterion remain visible but are excluded from N_use."
            ),
            "outputs": {
                key: str(path.relative_to(REPO_ROOT))
                if path.is_relative_to(REPO_ROOT)
                else str(path)
                for key, path in outputs.items()
            },
        },
    )
    outputs["metadata"] = metadata_path
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    for label, path in generate(args.config, args.output_dir).items():
        print(f"- {label}: {path}")


if __name__ == "__main__":
    main()
