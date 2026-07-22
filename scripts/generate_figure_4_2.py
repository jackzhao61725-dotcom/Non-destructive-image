"""Generate dissertation Figure 4.2 from the canonical condensate pipeline."""

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
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

from scripts.recover_notebook_condensate_stage import load_config
from scripts.recover_notebook_dgi_stage import build_dgi_stage
from scripts.recover_notebook_faraday_stage import build_faraday_stage
from scripts.recover_notebook_pci_stage import build_pci_stage


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = REPO_ROOT / "configs" / "figure_4_2.json"


def _load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _git_value(*args: str) -> str:
    commands = [
        ["git", *args],
        [r"C:\Program Files\Git\cmd\git.exe", *args],
    ]
    for command in commands:
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


def _stats(array: np.ndarray) -> dict[str, Any]:
    peak_index = np.unravel_index(int(np.argmax(array)), array.shape)
    centre = (array.shape[0] // 2, array.shape[1] // 2)
    return {
        "shape": [int(value) for value in array.shape],
        "min": float(np.min(array)),
        "max": float(np.max(array)),
        "mean": float(np.mean(array)),
        "centre_value": float(array[centre]),
        "peak_index": [int(value) for value in peak_index],
    }


def _assert_close(name: str, actual: float, expected: float, checks: dict[str, Any]) -> None:
    if not np.isclose(
        actual,
        expected,
        rtol=float(checks["relative_tolerance"]),
        atol=float(checks["absolute_tolerance"]),
    ):
        raise RuntimeError(f"Figure 4.2 canonical gate failed for {name}: {actual} != {expected}")


def build_figure_data(config: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Rebuild the cloud and propagate it through all four canonical readouts."""

    model_path = REPO_ROOT / config["model_config"]
    model = load_config(model_path)

    # Each recovery stage starts from the same condensate -> column density ->
    # scalar phase path and then applies the common pupil before its readout.
    pci_stage = build_pci_stage(model)
    dgi_stage = build_dgi_stage(model)
    faraday_stage = build_faraday_stage(model)

    pci_phase = np.asarray(pci_stage["phase_stage"]["notebook_phase_map_rad"])
    dgi_phase = np.asarray(dgi_stage["phase_stage"]["notebook_phase_map_rad"])
    faraday_phase = np.asarray(faraday_stage["phase_stage"]["notebook_phase_map_rad"])
    if not (np.array_equal(pci_phase, dgi_phase) and np.array_equal(pci_phase, faraday_phase)):
        raise RuntimeError("Figure 4.2 branches did not start from an identical scalar phase map")

    pci_pupil = np.asarray(pci_stage["notebook_pupil"])
    if not (
        np.array_equal(pci_pupil, dgi_stage["notebook_pupil"])
        and np.array_equal(pci_pupil, faraday_stage["notebook_pupil"])
    ):
        raise RuntimeError("Figure 4.2 branches did not use an identical finite-NA pupil")

    condensate = pci_stage["phase_stage"]["condensate_stage"]
    coordinate_um = np.asarray(condensate["coordinate_axis_m"]) * 1e6
    column_density = np.asarray(condensate["notebook_column_density_m2"])
    theta_f = np.asarray(faraday_stage["notebook_theta_f_map_rad"])
    pci = np.asarray(pci_stage["notebook_pci_image_intensity"])
    dgi = np.asarray(dgi_stage["notebook_dgi_image_intensity"])
    filtered_scattered = np.asarray(dgi_stage["notebook_propagated_scattered_field"])
    dark = np.asarray(faraday_stage["notebook_dark_field_intensity"])
    dual_s = np.asarray(faraday_stage["notebook_dual_port_signal"])

    # Dissertation notation: H is the historical notebook v port and V is u.
    port_h = np.asarray(faraday_stage["notebook_dual_port_v_intensity"])
    port_v = np.asarray(faraday_stage["notebook_dual_port_u_intensity"])

    scattered = np.exp(1j * pci_phase) - 1
    spectrum = np.fft.fft2(scattered)
    throughput = float(
        np.sum(np.abs(spectrum * pci_pupil) ** 2) / np.sum(np.abs(spectrum) ** 2)
    )
    surviving_rotation = float(np.arcsin(np.sqrt(np.max(dark))))
    centre = tuple(size // 2 for size in filtered_scattered.shape)
    phase_peak_index = np.unravel_index(int(np.argmax(pci_phase)), pci_phase.shape)
    theta_at_phase_peak = float(theta_f[phase_peak_index])
    dual_abs_peak_signed = float(dual_s.flat[int(np.argmax(np.abs(dual_s)))])
    centre_filtered_scattered = complex(filtered_scattered[centre])

    data = {
        "coordinate_um": coordinate_um,
        "column_density_m2": column_density,
        "phase_rad": pci_phase,
        "theta_f_rad": theta_f,
        "pupil": pci_pupil,
        "pci_I_I0": pci,
        "dgi_I_I0": dgi,
        "dark_I_I0": dark,
        "dual_S": dual_s,
        "port_H_I_I0": port_h,
        "port_V_I_I0": port_v,
    }

    metrics = {
        "grid": {
            "shape": [int(value) for value in pci.shape],
            "field_of_view_um": float(model["grid"]["field_of_view_m"]) * 1e6,
            "spacing_nm": float(model["grid"]["field_of_view_m"])
            / int(model["grid"]["ngrid"])
            * 1e9,
        },
        "condensate": {
            "imaging_axis": int(condensate["imaging_axis"]),
            "imaging_axis_label": model["imaging_geometry"]["axis_labels"][condensate["imaging_axis"]],
            "transverse_plane_labels": condensate["transverse_plane_labels"],
            "radii_um": [float(value) * 1e6 for value in condensate["state"]["radii"]],
            "peak_column_density_m2": float(pci_stage["phase_stage"]["column_density_peak_m2"]),
        },
        "optical_input": {
            "detuning_ghz": float(pci_stage["phase_stage"]["detuning_hz"]) / 1e9,
            "phase_peak_rad": float(np.max(pci_phase)),
            "kappa_F": float(faraday_stage["kappa_F"]),
            "theta_f_at_scalar_phase_peak_rad": theta_at_phase_peak,
            "theta_f_abs_peak_rad": float(np.max(np.abs(theta_f))),
            "numerical_aperture": float(pci_stage["pupil_stage"]["numerical_aperture"]),
            "scattered_power_throughput": throughput,
            "surviving_rotation_abs_rad": surviving_rotation,
            "surviving_rotation_abs_deg": float(np.degrees(surviving_rotation)),
            "centre_filtered_scattered_field": {
                "real": float(centre_filtered_scattered.real),
                "imag": float(centre_filtered_scattered.imag),
                "magnitude_squared": float(abs(centre_filtered_scattered) ** 2),
                "source": "dgi_stage.notebook_propagated_scattered_field at the image centre",
            },
        },
        "readouts": {
            "pci": {
                "background_I_I0": float(pci_stage["plate_background_intensity"]),
                "peak_signal_change_I_I0": float(np.max(pci) - pci_stage["plate_background_intensity"]),
                "stats": _stats(pci),
            },
            "dgi": {
                "background_I_I0": float(dgi_stage["dgi_reference_intensity"]),
                "peak_signal_change_I_I0": float(np.max(dgi) - dgi_stage["dgi_reference_intensity"]),
                "stats": _stats(dgi),
            },
            "dark_field_faraday": {
                "background_I_I0": 0.0,
                "stats": _stats(dark),
            },
            "dual_port_faraday": {
                "no_atom_background_H_I_I0": 0.5,
                "no_atom_background_V_I_I0": 0.5,
                "port_mapping": "I_H = historical notebook I_v; I_V = historical notebook I_u",
                "signal_definition": "S = (I_H - I_V)/(I_H + I_V)",
                "port_H_stats": _stats(port_h),
                "port_V_stats": _stats(port_v),
                "signal_stats": _stats(dual_s),
            },
        },
    }

    checks = config["canonical_checks"]
    actual_checks = {
        "phase_peak_rad": metrics["optical_input"]["phase_peak_rad"],
        "aperture_throughput": throughput,
        "pci_background_I_I0": metrics["readouts"]["pci"]["background_I_I0"],
        "pci_peak_I_I0": metrics["readouts"]["pci"]["stats"]["max"],
        "dgi_background_I_I0": metrics["readouts"]["dgi"]["background_I_I0"],
        "dgi_peak_I_I0": metrics["readouts"]["dgi"]["stats"]["max"],
        "faraday_dark_peak_I_I0": metrics["readouts"]["dark_field_faraday"]["stats"]["max"],
        "theta_f_at_phase_peak_rad": theta_at_phase_peak,
        "dual_port_centre_S": metrics["readouts"]["dual_port_faraday"]["signal_stats"]["centre_value"],
        "dual_port_min_S": metrics["readouts"]["dual_port_faraday"]["signal_stats"]["min"],
        "dual_port_max_S": metrics["readouts"]["dual_port_faraday"]["signal_stats"]["max"],
        "dual_port_abs_peak_signed_S": dual_abs_peak_signed,
        "dual_port_centre_I_H_I0": metrics["readouts"]["dual_port_faraday"]["port_H_stats"]["centre_value"],
        "dual_port_centre_I_V_I0": metrics["readouts"]["dual_port_faraday"]["port_V_stats"]["centre_value"],
        "surviving_rotation_abs_rad": surviving_rotation,
        "centre_filtered_scattered_field_real": centre_filtered_scattered.real,
        "centre_filtered_scattered_field_imag": centre_filtered_scattered.imag,
    }
    for name, actual in actual_checks.items():
        _assert_close(name, float(actual), float(checks[name]), checks)
    metrics["canonical_gate"] = {"passed": True, "values": actual_checks}
    return data, metrics


def _plot(
    data: dict[str, Any],
    config: dict[str, Any],
    svg_path: Path,
    png_path: Path,
    pdf_path: Path,
) -> None:
    figure_config = config["figure"]
    coordinate = np.asarray(data["coordinate_um"])
    extent = [float(coordinate[0]), float(coordinate[-1]), float(coordinate[0]), float(coordinate[-1])]

    plt.style.use("default")
    matplotlib.rcParams["svg.hashsalt"] = "figure-4-2-canonical-readouts"
    matplotlib.rcParams.update(
        {
            "font.size": 8.8,
            "axes.labelsize": 9.5,
            "axes.titlesize": 10.0,
            "xtick.labelsize": 8.2,
            "ytick.labelsize": 8.2,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )

    fig, axes = plt.subplots(
        2,
        2,
        figsize=tuple(figure_config["figsize_inches"]),
        sharex=True,
        sharey=True,
    )
    fig.subplots_adjust(left=0.085, right=0.91, bottom=0.14, top=0.93, wspace=0.48, hspace=0.46)
    panels = [
        (
            axes[0, 0],
            data["pci_I_I0"],
            "(a) PCI",
            "RdBu_r",
            TwoSlopeNorm(vmin=0.85, vcenter=0.9025, vmax=1.35),
            r"$I/I_0$",
            [0.85, 0.90, 1.10, 1.35],
            "%.2f",
        ),
        (
            axes[0, 1],
            data["dgi_I_I0"],
            "(b) DGI",
            "magma",
            matplotlib.colors.Normalize(vmin=0.0, vmax=0.045),
            r"$I/I_0$",
            [0.0, 0.02, 0.04],
            "%.3f",
        ),
        (
            axes[1, 0],
            data["dark_I_I0"],
            "(c) Dark-field Faraday",
            "magma",
            matplotlib.colors.Normalize(vmin=0.0, vmax=0.045),
            r"$I/I_0$",
            [0.0, 0.02, 0.04],
            "%.3f",
        ),
        (
            axes[1, 1],
            data["dual_S"],
            "(d) Dual-port Faraday",
            "RdBu_r",
            TwoSlopeNorm(vmin=-0.22, vcenter=0.0, vmax=0.22),
            r"$S$",
            [-0.2, 0.0, 0.2],
            "%.1f",
        ),
    ]

    for axis, image, title, cmap, norm, colorbar_label, ticks, tick_format in panels:
        plotted = axis.imshow(
            image,
            extent=extent,
            origin="lower",
            cmap=cmap,
            norm=norm,
            interpolation="nearest",
            rasterized=True,
        )
        axis.set_title(title, pad=5)
        axis.set_xlim(*figure_config["crop_y_um"])
        axis.set_ylim(*figure_config["crop_z_um"])
        axis.set_aspect("equal")
        axis.set_xticks([-30, -15, 0, 15, 30])
        axis.set_yticks([-5, 0, 5])
        divider = make_axes_locatable(axis)
        colorbar_axis = divider.append_axes("right", size="4%", pad=0.05)
        colorbar = fig.colorbar(plotted, cax=colorbar_axis, ticks=ticks)
        colorbar.set_label(colorbar_label, labelpad=4)
        colorbar.ax.yaxis.set_major_formatter(FormatStrFormatter(tick_format))

    axes[1, 0].set_xlabel(r"$y$ ($\mu\mathrm{m}$)")
    axes[1, 1].set_xlabel(r"$y$ ($\mu\mathrm{m}$)")
    axes[0, 0].set_ylabel(r"$z$ ($\mu\mathrm{m}$)")
    axes[1, 0].set_ylabel(r"$z$ ($\mu\mathrm{m}$)")

    fig.savefig(svg_path, format="svg", metadata={"Date": None})
    fig.savefig(
        png_path,
        format="png",
        dpi=int(figure_config["dpi"]),
        metadata={"Software": "matplotlib"},
    )
    fig.savefig(pdf_path, format="pdf", metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)

    svg_text = svg_path.read_text(encoding="utf-8")
    svg_path.write_text("\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n", encoding="utf-8")


def _write_lineouts(path: Path, data: dict[str, Any]) -> None:
    coordinate = np.asarray(data["coordinate_um"])
    middle = len(coordinate) // 2
    columns = {
        "y_um": coordinate,
        "column_density_m2": np.asarray(data["column_density_m2"])[middle],
        "phase_rad": np.asarray(data["phase_rad"])[middle],
        "theta_F_rad": np.asarray(data["theta_f_rad"])[middle],
        "pci_I_I0": np.asarray(data["pci_I_I0"])[middle],
        "dgi_I_I0": np.asarray(data["dgi_I_I0"])[middle],
        "dark_I_I0": np.asarray(data["dark_I_I0"])[middle],
        "dual_S": np.asarray(data["dual_S"])[middle],
        "port_H_I_I0": np.asarray(data["port_H_I_I0"])[middle],
        "port_V_I_I0": np.asarray(data["port_V_I_I0"])[middle],
    }
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for row in zip(*columns.values()):
            writer.writerow([f"{float(value):.12g}" for value in row])


def generate(config_path: Path = DEFAULT_CONFIG, output_dir: Path | None = None) -> dict[str, Path]:
    config = _load_json(config_path)
    data, metrics = build_figure_data(config)
    destination = output_dir or (REPO_ROOT / config["output_directory"])
    destination.mkdir(parents=True, exist_ok=True)

    figure_config = config["figure"]
    outputs = {
        "svg": destination / figure_config["svg_filename"],
        "png": destination / figure_config["png_filename"],
        "pdf": destination / figure_config["pdf_filename"],
        "lineouts": destination / "central_lineouts.csv",
        "values": destination / "figure_4_2_values.json",
        "metadata": destination / "metadata.json",
    }
    _plot(data, config, outputs["svg"], outputs["png"], outputs["pdf"])
    _write_lineouts(outputs["lineouts"], data)
    _write_json(outputs["values"], metrics)

    metadata = {
        "figure": "Figure 4.2",
        "figure_type": "canonical noiseless finite-aperture single-frame readout comparison",
        "pipeline": [
            "contact-only Thomas-Fermi condensate",
            "x-integrated column density",
            "scalar phase and ideal 166Er Faraday rotation maps",
            f"common effective NA={metrics['optical_input']['numerical_aperture']:.3f} pupil propagation",
            "PCI, DGI, dark-field Faraday, and dual-port Faraday readouts",
        ],
        "script": str(Path(__file__).resolve().relative_to(REPO_ROOT)).replace("\\", "/"),
        "config": str(config_path.resolve().relative_to(REPO_ROOT)).replace("\\", "/"),
        "model_config": config["model_config"],
        "git_branch": _git_value("branch", "--show-current"),
        "git_commit": _git_value("rev-parse", "HEAD"),
        "canonical_context": config["canonical_context"],
        "plotted_quantities": {
            "PCI": "noiseless image-plane intensity I/I0",
            "DGI": "noiseless image-plane intensity I/I0 with OD_s=4 residual carrier",
            "dark_field_faraday": "noiseless crossed-analyser intensity I/I0",
            "dual_port_faraday": "normalised difference S=(I_H-I_V)/(I_H+I_V)",
        },
        "port_mapping": "I_H = historical notebook I_v; I_V = historical notebook I_u",
        "camera_noise_included": False,
        "camera_binning_included": False,
        "multishot_evolution_included": False,
        "atomic_response_status": "ideal fully spin-polarised axial 166Er response with kappa_F=-45/91",
        "apparatus_calibration_status": "effective pupil, analyser-port gains, background and registration remain to be measured",
        "canonical_gate": metrics["canonical_gate"],
        "caption": (
            "Representative noiseless finite-aperture readouts of the reference condensate at "
            "|Delta|/2pi=1.5 GHz, imaged along x through the common effective NA=0.130 pupil. "
            "Panels (a)-(c) show image-plane intensity normalised to I0; panels (b) and (c) "
            "use the same intensity scale. Panel (d) shows the signed dual-port signal "
            "S=(I_H-I_V)/(I_H+I_V), whose negative central value follows the stated H/V "
            "convention. Camera sampling, detector noise and condensate evolution are omitted."
        ),
    }
    _write_json(outputs["metadata"], metadata)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    config_path = args.config if args.config.is_absolute() else REPO_ROOT / args.config
    output_dir = args.output_dir
    if output_dir is not None and not output_dir.is_absolute():
        output_dir = REPO_ROOT / output_dir
    outputs = generate(config_path, output_dir)
    for name, path in outputs.items():
        print(f"- {name}: {path.relative_to(REPO_ROOT) if path.is_relative_to(REPO_ROOT) else path}")


if __name__ == "__main__":
    main()
