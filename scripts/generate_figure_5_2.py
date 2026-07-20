"""Generate Figure 5.2: usable-frame screening across detuning and fluence."""

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
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
from scipy.special import zeta

from non_destructive_image import bin_to_camera_pixels, simulate_faraday_image
from scripts.recover_notebook_condensate_stage import load_config
from scripts.recover_notebook_multishot_stage import (
    _basic_constants,
    intensity_at_atoms_notebook,
    reabsorption_fraction,
    scalar_phase_peak,
    scattered_photons_per_atom,
    self_consistent_total_atoms,
    tf_state_for_atoms,
)
from scripts.recover_notebook_pci_stage import build_pupil


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = REPO_ROOT / "configs" / "figure_5_2.json"
MODES = ("Faraday dark-field", "Faraday dual-port")
MODE_TITLES = {
    "Faraday dark-field": "(a) Dark-field Faraday",
    "Faraday dual-port": "(b) Dual-port Faraday",
}


def _load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


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


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("Figure 5.2 produced no scan rows")
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


def _scan_axes(config: dict[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    scan = config["scan"]
    detuning = np.linspace(
        float(scan["detuning_min_ghz"]),
        float(scan["detuning_max_ghz"]),
        int(scan["detuning_points"]),
    )
    detuning = np.unique(
        np.sort(
            np.concatenate(
                [detuning, np.asarray(scan["required_detuning_ghz"], dtype=float)]
            )
        )
    )
    fluence = np.arange(
        float(scan["fluence_min_mw_us"]),
        float(scan["fluence_max_mw_us"]) + 0.5 * float(scan["fluence_step_mw_us"]),
        float(scan["fluence_step_mw_us"]),
    )
    fluence = np.unique(
        np.sort(
            np.concatenate(
                [fluence, np.asarray(scan["required_fluence_mw_us"], dtype=float)]
            )
        )
    )
    return detuning, fluence


def _central_block(array: np.ndarray, half_width: int) -> np.ndarray:
    centre = tuple(size // 2 for size in array.shape)
    return array[
        centre[0] - half_width : centre[0] + half_width + 1,
        centre[1] - half_width : centre[1] + half_width + 1,
    ]


def _photoelectron_scale(
    model: dict[str, Any],
    constants: dict[str, Any],
    fluence_mw_us: float,
) -> float:
    geometry = model["imaging_geometry"]
    camera = model["camera_recovery"]
    h_planck = 2 * np.pi * float(model["constants"]["hbar"])
    photon_energy = (
        h_planck
        * float(model["constants"]["speed_of_light"])
        / float(constants["wavelength"])
    )
    object_pixel_m = float(geometry["camera_pixel_m"]) / float(geometry["magnification"])
    return float(
        intensity_at_atoms_notebook(model, 1.0)
        * object_pixel_m**2
        * fluence_mw_us
        * 1e-6
        * float(camera["quantum_efficiency"])
        / photon_energy
    )


def _response_table(
    model: dict[str, Any],
    config: dict[str, Any],
    detuning_ghz: np.ndarray,
    condensate_fraction: np.ndarray,
) -> dict[str, np.ndarray]:
    constants = _basic_constants(model)
    scan = config["scan"]
    criteria = config["frame_criteria"]
    pupil = np.asarray(build_pupil(model)["pupil"], dtype=float)
    ngrid = int(model["grid"]["ngrid"])
    field_of_view_m = float(model["grid"]["field_of_view_m"])
    grid_axis = (np.arange(ngrid) - ngrid // 2) * (field_of_view_m / ngrid)
    coordinate_a, coordinate_b = np.meshgrid(grid_axis, grid_axis)
    imaging_axis = int(scan["imaging_axis"])
    plane = [index for index in range(3) if index != imaging_axis]
    bin_size = int(model["camera_recovery"]["bin_size"])
    half_width = int(criteria["central_block_half_width_pixels"])
    initial_condensate_atoms = float(constants["atom_number"])

    shape = (detuning_ghz.size, condensate_fraction.size)
    dark_signal = np.empty(shape, dtype=float)
    dual_signal = np.empty(shape, dtype=float)
    dual_total = np.empty(shape, dtype=float)

    for detuning_index, detuning_value in enumerate(detuning_ghz):
        detuning_hz = float(detuning_value * 1e9)
        for fraction_index, fraction in enumerate(condensate_fraction):
            state = tf_state_for_atoms(float(fraction * initial_condensate_atoms), constants)
            radii = np.asarray(state["radii"], dtype=float)
            profile = np.maximum(
                0.0,
                1.0
                - coordinate_a**2 / radii[plane[0]] ** 2
                - coordinate_b**2 / radii[plane[1]] ** 2,
            ) ** 1.5
            phase_peak = scalar_phase_peak(
                detuning_hz,
                float(np.asarray(state["column_density"])[imaging_axis]),
                constants,
            )
            faraday = simulate_faraday_image(
                float(scan["kappa_F"]) * phase_peak * profile,
                pupil,
            )
            dark = bin_to_camera_pixels(faraday["dark_field_intensity"], bin_size)
            # Dissertation convention: H is the historical notebook v port and V is u.
            port_h = bin_to_camera_pixels(faraday["dual_port_v_intensity"], bin_size)
            port_v = bin_to_camera_pixels(faraday["dual_port_u_intensity"], bin_size)
            dark_block = _central_block(np.asarray(dark, dtype=float), half_width)
            h_block = _central_block(np.asarray(port_h, dtype=float), half_width)
            v_block = _central_block(np.asarray(port_v, dtype=float), half_width)
            dark_signal[detuning_index, fraction_index] = float(np.sum(dark_block))
            dual_signal[detuning_index, fraction_index] = float(np.sum(h_block - v_block))
            dual_total[detuning_index, fraction_index] = float(np.sum(h_block + v_block))

    return {
        "dark_signal_per_i0": dark_signal,
        "dual_signal_per_i0": dual_signal,
        "dual_total_per_i0": dual_total,
    }


def _interpolated_snr(
    response: dict[str, np.ndarray],
    detuning_index: int,
    condensate_fraction: float,
    fraction_axis: np.ndarray,
    photoelectrons_per_i0_pixel: float,
    read_noise_e: float,
    block_pixels: int,
) -> dict[str, float]:
    dark_signal = float(
        np.interp(
            condensate_fraction,
            fraction_axis,
            response["dark_signal_per_i0"][detuning_index],
        )
    )
    dual_signal = float(
        np.interp(
            condensate_fraction,
            fraction_axis,
            response["dual_signal_per_i0"][detuning_index],
        )
    )
    dual_total = float(
        np.interp(
            condensate_fraction,
            fraction_axis,
            response["dual_total_per_i0"][detuning_index],
        )
    )
    count_scale = float(photoelectrons_per_i0_pixel)
    dark_snr = abs(dark_signal * count_scale) / np.sqrt(
        max(dark_signal, 0.0) * count_scale + block_pixels * read_noise_e**2
    )
    dual_snr = abs(dual_signal * count_scale) / np.sqrt(
        max(dual_total, 0.0) * count_scale + block_pixels * 2 * read_noise_e**2
    )
    return {
        "Faraday dark-field": float(dark_snr),
        "Faraday dual-port": float(dual_snr),
    }


def _sequence_result(
    model: dict[str, Any],
    config: dict[str, Any],
    constants: dict[str, Any],
    thermal: dict[str, float],
    response: dict[str, np.ndarray],
    detuning_index: int,
    detuning_ghz: float,
    fluence_mw_us: float,
    fraction_axis: np.ndarray,
) -> dict[str, Any]:
    scan = config["scan"]
    criteria = config["frame_criteria"]
    power_mw = float(scan["probe_power_mw"])
    exposure_s = fluence_mw_us / power_mw * 1e-6
    detuning_hz = detuning_ghz * 1e9
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
    temperature = float(model["condensate"]["temperature_k"])
    initial_condensate_atoms = float(constants["atom_number"])
    total_atoms = float(thermal["total_atoms"])
    count_scale = _photoelectron_scale(model, constants, fluence_mw_us)
    read_noise = float(model["camera_recovery"]["read_noise_electrons"])
    half_width = int(criteria["central_block_half_width_pixels"])
    block_pixels = (2 * half_width + 1) ** 2
    minimum_snr = float(criteria["minimum_snr"])
    loss_limit = float(criteria["maximum_post_exposure_condensate_loss"])

    depletion_frames = 0
    usable_frames = {mode: 0 for mode in MODES}
    quality_open = {mode: True for mode in MODES}
    framewise_snr = {mode: [] for mode in MODES}
    pre_loss: list[float] = []
    post_loss: list[float] = []
    frame_index = 0
    while frame_index < int(model["multishot_recovery"]["max_shots"]):
        condensate_atoms = total_atoms * (1 - (temperature / critical_temperature) ** 3)
        condensate_fraction = condensate_atoms / initial_condensate_atoms
        if condensate_fraction <= 0:
            break
        snr = _interpolated_snr(
            response,
            detuning_index,
            float(condensate_fraction),
            fraction_axis,
            count_scale,
            read_noise,
            block_pixels,
        )
        next_temperature = (
            temperature**4 + deposited_energy / energy_coefficient
        ) ** 0.25
        next_condensate_atoms = total_atoms * (
            1 - (next_temperature / critical_temperature) ** 3
        )
        next_loss = 1 - next_condensate_atoms / initial_condensate_atoms
        if next_loss > loss_limit + 1e-12:
            break

        depletion_frames += 1
        pre_loss.append(float(1 - condensate_fraction))
        post_loss.append(float(next_loss))
        for mode in MODES:
            framewise_snr[mode].append(float(snr[mode]))
            if quality_open[mode] and snr[mode] >= minimum_snr:
                usable_frames[mode] += 1
            else:
                quality_open[mode] = False
        temperature = next_temperature
        frame_index += 1

    result: dict[str, Any] = {
        "n_gamma": float(n_gamma),
        "reabsorption_fraction": float(reabsorption),
        "depletion_limited_frames": int(depletion_frames),
        "post_sequence_loss": float(post_loss[-1]) if post_loss else 0.0,
        "next_pulse_loss": float(next_loss) if "next_loss" in locals() else np.nan,
        "initial_snr": {},
        "final_depletion_accepted_snr": {},
        "usable_frames": {},
        "stopping_reason": {},
        "framewise_snr": framewise_snr,
    }
    for mode in MODES:
        values = framewise_snr[mode]
        result["initial_snr"][mode] = float(values[0]) if values else np.nan
        result["final_depletion_accepted_snr"][mode] = (
            float(values[-1]) if values else np.nan
        )
        result["usable_frames"][mode] = int(usable_frames[mode])
        result["stopping_reason"][mode] = (
            "image quality" if usable_frames[mode] < depletion_frames else "depletion"
        )
    return result


def build_scan(config: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    model_path = REPO_ROOT / config["model_config"]
    model = load_config(model_path)
    constants = _basic_constants(model)
    thermal = self_consistent_total_atoms(model, constants)
    detuning, fluence = _scan_axes(config)
    interpolation = config["response_interpolation"]
    fraction_axis = np.linspace(
        float(interpolation["minimum_condensate_fraction"]),
        float(interpolation["maximum_condensate_fraction"]),
        int(interpolation["fraction_points"]),
    )
    response = _response_table(model, config, detuning, fraction_axis)

    shape = (detuning.size, fluence.size)
    depletion_frames = np.zeros(shape, dtype=int)
    usable = {mode: np.zeros(shape, dtype=int) for mode in MODES}
    reasons = {mode: np.empty(shape, dtype=object) for mode in MODES}
    initial_snr = {mode: np.zeros(shape, dtype=float) for mode in MODES}
    final_snr = {mode: np.zeros(shape, dtype=float) for mode in MODES}
    rows: list[dict[str, Any]] = []
    details: dict[tuple[float, float], dict[str, Any]] = {}

    for detuning_index, detuning_value in enumerate(detuning):
        for fluence_index, fluence_value in enumerate(fluence):
            result = _sequence_result(
                model,
                config,
                constants,
                thermal,
                response,
                detuning_index,
                float(detuning_value),
                float(fluence_value),
                fraction_axis,
            )
            details[(float(detuning_value), float(fluence_value))] = result
            depletion_frames[detuning_index, fluence_index] = result[
                "depletion_limited_frames"
            ]
            for mode in MODES:
                usable[mode][detuning_index, fluence_index] = result["usable_frames"][mode]
                reasons[mode][detuning_index, fluence_index] = result["stopping_reason"][mode]
                initial_snr[mode][detuning_index, fluence_index] = result["initial_snr"][mode]
                final_snr[mode][detuning_index, fluence_index] = result[
                    "final_depletion_accepted_snr"
                ][mode]
                rows.append(
                    {
                        "detuning_ghz": f"{detuning_value:.12g}",
                        "fluence_mw_us": f"{fluence_value:.12g}",
                        "mode": mode,
                        "minimum_snr": f"{float(config['frame_criteria']['minimum_snr']):.12g}",
                        "initial_snr": f"{result['initial_snr'][mode]:.12g}",
                        "final_depletion_accepted_snr": f"{result['final_depletion_accepted_snr'][mode]:.12g}",
                        "depletion_limited_frames": result["depletion_limited_frames"],
                        "usable_frames": result["usable_frames"][mode],
                        "stopping_reason": result["stopping_reason"][mode],
                        "post_sequence_loss": f"{result['post_sequence_loss']:.12g}",
                        "next_pulse_loss": f"{result['next_pulse_loss']:.12g}",
                        "n_gamma": f"{result['n_gamma']:.12g}",
                        "reabsorption_fraction": f"{result['reabsorption_fraction']:.12g}",
                    }
                )

    reference_key = (1.5, 90.0)
    if reference_key not in details:
        raise RuntimeError("Figure 5.2 scan does not contain the 1.5 GHz, 90 mW us reference")
    reference = details[reference_key]
    if reference["depletion_limited_frames"] != 10:
        raise RuntimeError("Figure 5.2 canonical depletion count is not 10 frames")
    expected_initial = {
        "Faraday dark-field": 6.22710517196034,
        "Faraday dual-port": 16.06033992848492,
    }
    for mode, expected in expected_initial.items():
        if not np.isclose(reference["initial_snr"][mode], expected, rtol=5e-5, atol=1e-8):
            raise RuntimeError(
                f"Figure 5.2 initial {mode} SNR does not reproduce Figure 5.1"
            )

    summary = {
        "canonical_gate": {"passed": True, "reference": reference},
        "scan": {
            "detuning_ghz": detuning,
            "fluence_mw_us": fluence,
            "usable_frame_range": {
                mode: [int(np.min(values)), int(np.max(values))]
                for mode, values in usable.items()
            },
            "depletion_frame_range": [
                int(np.min(depletion_frames)),
                int(np.max(depletion_frames)),
            ],
        },
        "criteria": config["frame_criteria"],
        "selected_conditions": {
            label: details[(1.5, fluence_value)]
            for label, fluence_value in zip(
                ("A", "B", "C"),
                (30.0, 90.0, 150.0),
                strict=True,
            )
        },
        "range_rationale": config["range_rationale"],
        "scope": config["scope"],
    }
    data = {
        "detuning_ghz": detuning,
        "fluence_mw_us": fluence,
        "depletion_frames": depletion_frames,
        "usable_frames": usable,
        "stopping_reason": reasons,
        "initial_snr": initial_snr,
        "final_snr": final_snr,
        "rows": rows,
    }
    return data, summary


def _frame_colormap(maximum_frames: int) -> tuple[ListedColormap, BoundaryNorm]:
    base = plt.get_cmap("viridis")
    colours = [(0.92, 0.92, 0.92, 1.0)]
    colours.extend(base(np.linspace(0.12, 0.96, maximum_frames)))
    cmap = ListedColormap(colours, name="usable_frames")
    boundaries = np.arange(-0.5, maximum_frames + 1.5, 1.0)
    return cmap, BoundaryNorm(boundaries, cmap.N, clip=True)


def _plot(data: dict[str, Any], config: dict[str, Any], paths: dict[str, Path]) -> None:
    figure = config["figure"]
    detuning = np.asarray(data["detuning_ghz"], dtype=float)
    fluence = np.asarray(data["fluence_mw_us"], dtype=float)
    maximum_frames = int(figure["maximum_displayed_frames"])
    if any(np.max(values) > maximum_frames for values in data["usable_frames"].values()):
        raise RuntimeError("Figure 5.2 colour scale clips the usable-frame count")
    cmap, norm = _frame_colormap(maximum_frames)

    plt.style.use("default")
    matplotlib.rcParams["svg.hashsalt"] = "figure-5-2-usable-frame-screen"
    matplotlib.rcParams.update(
        {
            "font.size": 8.8,
            "axes.labelsize": 9.2,
            "axes.titlesize": 9.3,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "hatch.linewidth": 0.45,
        }
    )
    fig, axes = plt.subplots(
        1,
        2,
        figsize=tuple(figure["figsize_inches"]),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    image = None
    for axis, mode in zip(axes, MODES, strict=True):
        values = np.asarray(data["usable_frames"][mode], dtype=int)
        image = axis.pcolormesh(
            fluence,
            detuning,
            values,
            shading="nearest",
            cmap=cmap,
            norm=norm,
            rasterized=False,
        )
        image_limited = (
            (np.asarray(data["stopping_reason"][mode], dtype=object) == "image quality")
            & (values > 0)
        )
        hatch = axis.contourf(
            fluence,
            detuning,
            image_limited.astype(float),
            levels=[0.5, 1.5],
            colors="none",
            hatches=["///"],
        )
        hatch.set_edgecolor("0.32")
        hatch.set_linewidth(0.0)
        axis.axhline(1.5, color="white", linewidth=0.9, linestyle=":", alpha=0.95)
        for label, marker, fluence_value in zip(
            ("A", "B", "C"),
            ("o", "s", "^"),
            (30.0, 90.0, 150.0),
            strict=True,
        ):
            axis.plot(
                fluence_value,
                1.5,
                marker=marker,
                markersize=4.7,
                markerfacecolor="white",
                markeredgecolor="black",
                markeredgewidth=0.8,
                linestyle="none",
                zorder=5,
            )
            axis.annotate(
                label,
                (fluence_value, 1.5),
                xytext=(0, 6),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=7.7,
                fontweight="bold",
                color="black",
                zorder=6,
            )
        axis.set_title(MODE_TITLES[mode], loc="left", pad=5)
        axis.set_xlabel(r"Fluence coordinate $F=P\tau$ (mW $\mu$s)")
        axis.set_xlim(float(fluence[0]), float(fluence[-1]))
        axis.set_ylim(float(detuning[0]), float(detuning[-1]))
        axis.set_xticks([30, 60, 90, 120, 150, 180])
        axis.set_yticks([0.75, 1.0, 1.5, 2.0, 2.5, 3.0])
    axes[0].set_ylabel(r"Absolute detuning $|\Delta|/2\pi$ (GHz)")
    if image is None:
        raise RuntimeError("Figure 5.2 produced no map")
    colourbar = fig.colorbar(image, ax=axes, fraction=0.035, pad=0.025)
    colourbar.set_label("Usable frames")
    colourbar.set_ticks([0, 5, 10, 15, 20])
    legend = [
        Patch(
            facecolor="none",
            edgecolor="0.32",
            hatch="///",
            label="image-quality limited",
        )
    ]
    axes[1].legend(
        handles=legend,
        frameon=True,
        framealpha=0.9,
        fontsize=7.4,
        loc="upper left",
        borderpad=0.35,
        handlelength=1.5,
    )

    fig.savefig(paths["svg"], bbox_inches="tight")
    fig.savefig(paths["png"], dpi=int(figure["dpi"]), bbox_inches="tight")
    fig.savefig(paths["pdf"], bbox_inches="tight")
    plt.close(fig)


def _plot_dual_port_heatmap(
    data: dict[str, Any],
    config: dict[str, Any],
    paths: dict[str, Path],
) -> None:
    figure = config["figure"]
    mode = "Faraday dual-port"
    detuning = np.asarray(data["detuning_ghz"], dtype=float)
    fluence = np.asarray(data["fluence_mw_us"], dtype=float)
    values = np.asarray(data["usable_frames"][mode], dtype=int)
    maximum_frames = int(figure["maximum_displayed_frames"])
    if np.max(values) > maximum_frames:
        raise RuntimeError("Dual-port heatmap colour scale clips the usable-frame count")
    cmap, norm = _frame_colormap(maximum_frames)

    relative_fraction = float(figure["high_n_relative_fraction"])
    row_maximum = np.max(values, axis=1, keepdims=True)
    high_n_threshold = np.ceil(relative_fraction * row_maximum)
    high_n_band = (values > 0) & (values >= high_n_threshold)

    plt.style.use("default")
    matplotlib.rcParams["svg.hashsalt"] = "figure-5-2-dual-port-heatmap"
    matplotlib.rcParams.update(
        {
            "font.size": 9.0,
            "axes.labelsize": 9.8,
            "axes.titlesize": 9.8,
            "xtick.labelsize": 8.4,
            "ytick.labelsize": 8.4,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    fig, axis = plt.subplots(
        figsize=tuple(figure["dual_port_heatmap_figsize_inches"]),
        constrained_layout=True,
    )
    image = axis.pcolormesh(
        fluence,
        detuning,
        values,
        shading="nearest",
        cmap=cmap,
        norm=norm,
        rasterized=False,
    )

    # A white underlay keeps the high-N boundary visible across the full colour scale.
    axis.contour(
        fluence,
        detuning,
        high_n_band.astype(float),
        levels=[0.5],
        colors="white",
        linewidths=3.0,
    )
    axis.contour(
        fluence,
        detuning,
        high_n_band.astype(float),
        levels=[0.5],
        colors="black",
        linewidths=1.15,
    )
    axis.axhline(1.5, color="white", linewidth=1.0, linestyle=":", alpha=0.95)
    for label, marker, fluence_value in zip(
        ("A", "B", "C"),
        ("o", "s", "^"),
        (30.0, 90.0, 150.0),
        strict=True,
    ):
        axis.plot(
            fluence_value,
            1.5,
            marker=marker,
            markersize=5.2,
            markerfacecolor="white",
            markeredgecolor="black",
            markeredgewidth=0.9,
            linestyle="none",
            zorder=6,
        )
        axis.annotate(
            label,
            (fluence_value, 1.5),
            xytext=(0, 7),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8.0,
            fontweight="bold",
            zorder=7,
        )

    axis.set_xlabel(r"Fluence coordinate $F=P\tau$ (mW $\mu$s)")
    axis.set_ylabel(r"Absolute detuning $|\Delta|/2\pi$ (GHz)")
    axis.set_xlim(float(fluence[0]), float(fluence[-1]))
    axis.set_ylim(float(detuning[0]), float(detuning[-1]))
    axis.set_xticks([30, 60, 90, 120, 150, 180])
    axis.set_yticks([0.75, 1.0, 1.5, 2.0, 2.5, 3.0])
    colourbar = fig.colorbar(image, ax=axis, fraction=0.045, pad=0.025)
    colourbar.set_label("Usable frames")
    colourbar.set_ticks([0, 5, 10, 15, 20])

    handles = [
        Line2D(
            [0],
            [0],
            color="black",
            linewidth=1.4,
            label=rf"high-$N$ band ($\geq {100 * relative_fraction:.0f}\%$ of maximum at each detuning)",
        )
    ]
    axis.legend(
        handles=handles,
        frameon=True,
        framealpha=0.92,
        fontsize=7.7,
        loc="upper left",
        borderpad=0.4,
        handlelength=2.2,
    )

    fig.savefig(paths["svg"], bbox_inches="tight")
    fig.savefig(paths["png"], dpi=int(figure["dpi"]), bbox_inches="tight")
    fig.savefig(paths["pdf"], bbox_inches="tight")
    plt.close(fig)


def _plot_operating_band(
    data: dict[str, Any],
    config: dict[str, Any],
    paths: dict[str, Path],
) -> None:
    figure = config["figure"]
    mode = "Faraday dual-port"
    detuning = np.asarray(data["detuning_ghz"], dtype=float)
    fluence = np.asarray(data["fluence_mw_us"], dtype=float)
    target_detuning = float(figure["operating_band_detuning_ghz"])
    detuning_index = int(np.argmin(np.abs(detuning - target_detuning)))
    if not np.isclose(detuning[detuning_index], target_detuning, atol=1e-12, rtol=0):
        raise RuntimeError("Operating-band detuning is absent from the Figure 5.2 scan")
    minimum_fluence = float(figure["operating_band_fluence_min_mw_us"])
    mask = fluence >= minimum_fluence
    x = fluence[mask]
    depletion = np.asarray(data["depletion_frames"], dtype=int)[detuning_index, mask]
    usable = np.asarray(data["usable_frames"][mode], dtype=int)[detuning_index, mask]

    relative_fraction = float(figure["high_n_relative_fraction"])
    high_threshold = int(np.ceil(relative_fraction * np.max(usable)))
    high_mask = usable >= high_threshold
    if not np.any(high_mask):
        raise RuntimeError("Operating-band cross-section contains no high-N interval")
    step = float(np.median(np.diff(x)))
    band_left = float(x[high_mask][0] - step / 2)
    band_right = float(x[high_mask][-1] + step / 2)

    plt.style.use("default")
    matplotlib.rcParams["svg.hashsalt"] = "figure-5-2-operating-band"
    matplotlib.rcParams.update(
        {
            "font.size": 9.0,
            "axes.labelsize": 9.8,
            "xtick.labelsize": 8.4,
            "ytick.labelsize": 8.4,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    fig, axis = plt.subplots(
        figsize=tuple(figure["operating_band_figsize_inches"]),
        constrained_layout=True,
    )
    axis.fill_between(
        x,
        usable,
        depletion,
        color="0.88",
        alpha=0.85,
        label="frames rejected by image-quality criterion",
    )
    axis.axvspan(
        band_left,
        band_right,
        color="#F2C14E",
        alpha=0.24,
        linewidth=0,
        label=rf"high-$N$ band ($N_{{\rm use}}\geq{high_threshold}$)",
    )
    axis.plot(
        x,
        depletion,
        color="0.35",
        linewidth=1.7,
        linestyle="--",
        label=r"strict depletion limit $N_{\rm dep}$",
    )
    axis.plot(
        x,
        usable,
        color="#2B7A73",
        linewidth=2.1,
        label=rf"usable frames ($\mathrm{{SNR}}_{{3\times3}}\geq{float(config['frame_criteria']['minimum_snr']):g}$)",
    )
    for label, marker, fluence_value in zip(
        ("A", "B", "C"),
        ("o", "s", "^"),
        (30.0, 90.0, 150.0),
        strict=True,
    ):
        index = int(np.flatnonzero(np.isclose(x, fluence_value, atol=1e-12, rtol=0))[0])
        axis.plot(
            fluence_value,
            usable[index],
            marker=marker,
            markersize=5.5,
            markerfacecolor="white",
            markeredgecolor="black",
            markeredgewidth=0.9,
            linestyle="none",
            zorder=6,
        )
        x_offset = 5 if label == "A" else 0
        axis.annotate(
            label,
            (fluence_value, usable[index]),
            xytext=(x_offset, 7),
            textcoords="offset points",
            ha="left" if label == "A" else "center",
            va="bottom",
            fontsize=8.0,
            fontweight="bold",
        )

    axis.set_xlim(float(x[0]) - 2, float(x[-1]))
    axis.set_ylim(0, max(35, int(np.max(depletion)) + 2))
    axis.set_xticks([30, 60, 90, 120, 150, 180])
    axis.set_xlabel(r"Fluence coordinate $F=P\tau$ (mW $\mu$s)")
    axis.set_ylabel("Consecutive frames")
    axis.grid(axis="y", color="0.9", linewidth=0.6)
    handles, labels = axis.get_legend_handles_labels()
    order = [3, 2, 0, 1]
    axis.legend(
        [handles[index] for index in order],
        [labels[index] for index in order],
        frameon=False,
        fontsize=7.8,
        loc="upper right",
        handlelength=2.3,
    )

    fig.savefig(paths["svg"], bbox_inches="tight")
    fig.savefig(paths["png"], dpi=int(figure["dpi"]), bbox_inches="tight")
    fig.savefig(paths["pdf"], bbox_inches="tight")
    plt.close(fig)


def generate(
    config_path: Path = DEFAULT_CONFIG,
    output_dir: Path | None = None,
) -> dict[str, Path]:
    config_path = config_path if config_path.is_absolute() else REPO_ROOT / config_path
    config = _load_json(config_path)
    destination = output_dir or REPO_ROOT / config["output_directory"]
    destination.mkdir(parents=True, exist_ok=True)
    paths = {
        "svg": destination / config["figure"]["svg_filename"],
        "png": destination / config["figure"]["png_filename"],
        "pdf": destination / config["figure"]["pdf_filename"],
        "dual_port_heatmap_svg": destination
        / config["figure"]["dual_port_heatmap_svg_filename"],
        "dual_port_heatmap_png": destination
        / config["figure"]["dual_port_heatmap_png_filename"],
        "dual_port_heatmap_pdf": destination
        / config["figure"]["dual_port_heatmap_pdf_filename"],
        "operating_band_svg": destination
        / config["figure"]["operating_band_svg_filename"],
        "operating_band_png": destination
        / config["figure"]["operating_band_png_filename"],
        "operating_band_pdf": destination
        / config["figure"]["operating_band_pdf_filename"],
        "csv": destination / "figure_5_2_data.csv",
        "values": destination / "figure_5_2_values.json",
        "metadata": destination / "metadata.json",
    }
    data, summary = build_scan(config)
    _write_csv(paths["csv"], data["rows"])
    _write_json(paths["values"], summary)
    _plot(data, config, paths)
    _plot_dual_port_heatmap(
        data,
        config,
        {
            "svg": paths["dual_port_heatmap_svg"],
            "png": paths["dual_port_heatmap_png"],
            "pdf": paths["dual_port_heatmap_pdf"],
        },
    )
    _plot_operating_band(
        data,
        config,
        {
            "svg": paths["operating_band_svg"],
            "png": paths["operating_band_png"],
            "pdf": paths["operating_band_pdf"],
        },
    )
    _write_json(
        paths["metadata"],
        {
            "label": config["label"],
            "config_path": str(config_path.relative_to(REPO_ROOT)),
            "model_config": config["model_config"],
            "git_branch": _git_value("branch", "--show-current"),
            "git_commit": _git_value("rev-parse", "HEAD"),
            "outputs": {
                key: str(path.relative_to(REPO_ROOT))
                if path.is_relative_to(REPO_ROOT)
                else str(path)
                for key, path in paths.items()
            },
            "criteria": summary["criteria"],
            "range_rationale": summary["range_rationale"],
            "scope": summary["scope"],
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
