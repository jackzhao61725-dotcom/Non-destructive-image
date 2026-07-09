"""Generate configurable Version 1 dissertation figure gallery outputs."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from html import escape
from pathlib import Path
from typing import Any

import numpy as np

from non_destructive_image import (
    accumulate_snr,
    compute_optical_density,
    estimate_cloud_moments,
    reabsorption_fraction,
    scalar_phase_shift,
    scattered_photons_per_atom,
    simulate_camera_image,
    simulate_dgi_image,
    simulate_faraday_image,
    simulate_multishot_sequence,
    simulate_noisy_camera_image,
    simulate_pci_image,
    summarise_faraday_sweep,
    sweep_faraday_detuning,
    sweep_faraday_exposure_time,
    sweep_faraday_intensity,
    thomas_fermi_profile_2d,
)


HELPERS_USED = [
    "thomas_fermi_profile_2d",
    "scalar_phase_shift",
    "scattered_photons_per_atom",
    "reabsorption_fraction",
    "simulate_pci_image",
    "simulate_dgi_image",
    "simulate_faraday_image",
    "simulate_camera_image",
    "simulate_noisy_camera_image",
    "simulate_multishot_sequence",
    "accumulate_snr",
    "sweep_faraday_detuning",
    "sweep_faraday_intensity",
    "sweep_faraday_exposure_time",
    "summarise_faraday_sweep",
    "compute_optical_density",
    "estimate_cloud_moments",
]


def _load_config(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _git_commit() -> str:
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


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    def convert(value: Any) -> Any:
        if isinstance(value, np.ndarray):
            return value.tolist()
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


def _write_rows(path: Path, rows: list[dict[str, float]]) -> None:
    if not rows:
        raise ValueError(f"no rows to write for {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: f"{value:.12g}" if isinstance(value, (float, int, np.floating, np.integer)) else value
                    for key, value in row.items()
                }
            )


def _hex_to_rgb(hex_colour: str) -> tuple[int, int, int]:
    hex_colour = hex_colour.lstrip("#")
    return tuple(int(hex_colour[i : i + 2], 16) for i in (0, 2, 4))


def _rgb_to_hex(rgb: tuple[int, int, int]) -> str:
    return "#" + "".join(f"{channel:02x}" for channel in rgb)


def _interpolate_colour(value: float, colours: list[str]) -> str:
    value = float(np.clip(value, 0.0, 1.0))
    if value == 1:
        return colours[-1]
    position = value * (len(colours) - 1)
    index = int(position)
    fraction = position - index
    left = _hex_to_rgb(colours[index])
    right = _hex_to_rgb(colours[index + 1])
    return _rgb_to_hex(tuple(int(left[i] + fraction * (right[i] - left[i])) for i in range(3)))


def _normalise(values: np.ndarray) -> np.ndarray:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return np.zeros_like(values, dtype=float)
    minimum = float(np.min(finite))
    maximum = float(np.max(finite))
    if minimum == maximum:
        return np.ones_like(values, dtype=float)
    return (values - minimum) / (maximum - minimum)


def _polyline(points: list[tuple[float, float]], colour: str, width: float = 2.4) -> str:
    return (
        '<polyline points="'
        + " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
        + f'" fill="none" stroke="{colour}" stroke-width="{width}" />'
    )


def _line_points(
    x_values: np.ndarray,
    y_values: np.ndarray,
    *,
    left: float,
    top: float,
    width: float,
    height: float,
) -> list[tuple[float, float]]:
    x_min = float(np.min(x_values))
    x_max = float(np.max(x_values))
    y_min = float(np.min(y_values))
    y_max = float(np.max(y_values))
    x_span = x_max - x_min if x_max != x_min else 1.0
    y_span = y_max - y_min if y_max != y_min else 1.0
    return [
        (
            left + ((float(x) - x_min) / x_span) * width,
            top + (1 - (float(y) - y_min) / y_span) * height,
        )
        for x, y in zip(x_values, y_values)
    ]


def _svg_header(width: int, height: int, title: str, label: str) -> str:
    return f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
  <title>{escape(title)}</title>
  <desc>Version 1 representative / uncalibrated figure placeholder generated for {escape(label)}.</desc>
  <rect width="100%" height="100%" fill="white" />
"""


def _svg_footer() -> str:
    return "</svg>\n"


def _heatmap_rects(
    data: np.ndarray,
    *,
    x: float,
    y: float,
    width: float,
    height: float,
    colours: list[str],
    max_cells: int = 48,
    vmin: float | None = None,
    vmax: float | None = None,
) -> str:
    array = np.asarray(data, dtype=float)
    step_y = max(1, int(np.ceil(array.shape[0] / max_cells)))
    step_x = max(1, int(np.ceil(array.shape[1] / max_cells)))
    reduced = array[::step_y, ::step_x]
    low = float(np.nanmin(reduced)) if vmin is None else vmin
    high = float(np.nanmax(reduced)) if vmax is None else vmax
    span = high - low if high != low else 1.0
    cell_w = width / reduced.shape[1]
    cell_h = height / reduced.shape[0]
    parts = []
    for row in range(reduced.shape[0]):
        for col in range(reduced.shape[1]):
            norm = (float(reduced[row, col]) - low) / span
            colour = _interpolate_colour(norm, colours)
            parts.append(
                f'<rect x="{x + col * cell_w:.2f}" y="{y + row * cell_h:.2f}" '
                f'width="{cell_w + 0.2:.2f}" height="{cell_h + 0.2:.2f}" fill="{colour}" />'
            )
    return "\n".join(parts)


def _write_closed_loop(path: Path, label: str) -> None:
    width = 960
    height = 360
    boxes = [
        ("Notebook-equivalent\nphysical model", 42, 120),
        ("Design decision\ndetuning, power,\nexposure, frames", 210, 120),
        ("Faraday imaging\nexperiment", 400, 120),
        ("Absorption / RAI\nfeedback", 570, 120),
        ("Calibration\nOD, density,\nkappa_F", 730, 120),
        ("Updated\noptimisation", 830, 250),
    ]
    svg = [_svg_header(width, height, "Closed-loop calibration architecture", label)]
    svg.append('<text x="480" y="34" text-anchor="middle" font-family="Arial" font-size="20">Closed-loop calibration architecture</text>')
    svg.append('<text x="480" y="58" text-anchor="middle" font-family="Arial" font-size="12">Version 1 representative / uncalibrated planning schematic</text>')
    for text, x, y in boxes:
        svg.append(f'<rect x="{x}" y="{y}" width="130" height="78" rx="7" fill="#f7fbff" stroke="#2d5c88" stroke-width="1.5" />')
        for i, line in enumerate(text.split("\n")):
            svg.append(f'<text x="{x+65}" y="{y+24+i*16}" text-anchor="middle" font-family="Arial" font-size="12">{escape(line)}</text>')
    arrows = [
        (172, 159, 210, 159),
        (340, 159, 400, 159),
        (530, 159, 570, 159),
        (700, 159, 730, 159),
        (795, 198, 860, 250),
        (830, 289, 107, 198),
    ]
    for x1, y1, x2, y2 in arrows:
        svg.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="#333" stroke-width="1.6" marker-end="url(#arrow)" />')
    svg.insert(
        1,
        '<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 Z" fill="#333" /></marker></defs>',
    )
    svg.append('<text x="480" y="330" text-anchor="middle" font-family="Arial" font-size="12">The simulator is a calibration-aware optimisation loop, not a one-way predictor.</text>')
    svg.append(_svg_footer())
    path.write_text("\n".join(svg), encoding="utf-8")


def _write_density_figure(path: Path, axis_m: np.ndarray, density: np.ndarray, label: str) -> None:
    width = 900
    height = 420
    mid = density.shape[0] // 2
    x_um = axis_m * 1e6
    line = density[mid] / np.max(density)
    points = _line_points(x_um, line, left=520, top=110, width=320, height=210)
    svg = [_svg_header(width, height, "Thomas-Fermi density model", label)]
    svg.append('<text x="450" y="32" text-anchor="middle" font-family="Arial" font-size="20">Representative Thomas-Fermi column-density model</text>')
    svg.append('<text x="450" y="55" text-anchor="middle" font-family="Arial" font-size="12">Version 1 representative / uncalibrated density object</text>')
    svg.append(_heatmap_rects(density, x=70, y=90, width=300, height=300, colours=["#081d58", "#225ea8", "#41b6c4", "#c7e9b4", "#ffffd9"]))
    svg.append('<rect x="70" y="90" width="300" height="300" fill="none" stroke="#222" />')
    svg.append('<text x="220" y="410" text-anchor="middle" font-family="Arial" font-size="13">x / y field of view</text>')
    svg.append('<line x1="520" y1="320" x2="840" y2="320" stroke="black" />')
    svg.append('<line x1="520" y1="110" x2="520" y2="320" stroke="black" />')
    svg.append(_polyline(points, "#1f77b4"))
    svg.append('<text x="680" y="360" text-anchor="middle" font-family="Arial" font-size="13">position (um)</text>')
    svg.append('<text x="455" y="220" text-anchor="middle" font-family="Arial" font-size="13" transform="rotate(-90 455 220)">normalised centre line</text>')
    svg.append(_svg_footer())
    path.write_text("\n".join(svg), encoding="utf-8")


def _write_profile_lineouts(path: Path, axis_m: np.ndarray, profile: np.ndarray, label: str) -> None:
    width = 900
    height = 430
    mid = profile.shape[0] // 2
    x_um = axis_m * 1e6
    x_line = profile[mid, :]
    y_line = profile[:, mid]
    left, top, plot_w, plot_h = 100, 92, 700, 255
    svg = [_svg_header(width, height, "Thomas-Fermi profile lineouts", label)]
    svg.append('<text x="450" y="32" text-anchor="middle" font-family="Arial" font-size="20">Representative cloud profile lineouts</text>')
    svg.append('<text x="450" y="55" text-anchor="middle" font-family="Arial" font-size="12">Version 1 representative / uncalibrated Thomas-Fermi object</text>')
    svg.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="black" />')
    svg.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="black" />')
    svg.append(_polyline(_line_points(x_um, x_line, left=left, top=top, width=plot_w, height=plot_h), "#1f77b4"))
    svg.append(_polyline(_line_points(x_um, y_line, left=left, top=top, width=plot_w, height=plot_h), "#ff7f0e", 2.0))
    svg.append('<line x1="610" y1="112" x2="652" y2="112" stroke="#1f77b4" stroke-width="2.5" />')
    svg.append('<text x="660" y="116" font-family="Arial" font-size="12">x centre line</text>')
    svg.append('<line x1="610" y1="134" x2="652" y2="134" stroke="#ff7f0e" stroke-width="2.5" />')
    svg.append('<text x="660" y="138" font-family="Arial" font-size="12">y centre line</text>')
    svg.append('<text x="450" y="390" text-anchor="middle" font-family="Arial" font-size="13">position (um)</text>')
    svg.append('<text x="38" y="220" text-anchor="middle" font-family="Arial" font-size="13" transform="rotate(-90 38 220)">normalised column density</text>')
    svg.append(_svg_footer())
    path.write_text("\n".join(svg), encoding="utf-8")


def _write_tradeoff_figure(path: Path, rows: list[dict[str, float]], label: str) -> None:
    width = 900
    height = 460
    x = np.asarray([row["detuning_hz"] / 1e9 for row in rows])
    signal = _normalise(np.asarray([row["signal_scale"] for row in rows]))
    scatter = _normalise(np.asarray([row["scattered_photons_per_atom"] for row in rows]))
    efficiency = _normalise(np.asarray([row["signal_per_scattered_photon"] for row in rows]))
    left, top, plot_w, plot_h = 90, 90, 720, 280
    svg = [_svg_header(width, height, "Light-atom signal/destructiveness trade-off", label)]
    svg.append('<text x="450" y="32" text-anchor="middle" font-family="Arial" font-size="20">Light-atom signal / destructiveness trade-off</text>')
    svg.append('<text x="450" y="55" text-anchor="middle" font-family="Arial" font-size="12">Normalised curves; Version 1 representative / uncalibrated</text>')
    svg.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="black" />')
    svg.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="black" />')
    for values, colour, name, y in [
        (signal, "#1f77b4", "Faraday/dispersive signal", 106),
        (scatter, "#d62728", "scattered photons / destructiveness", 128),
        (efficiency, "#2ca02c", "signal per scattered photon", 150),
    ]:
        svg.append(_polyline(_line_points(x, values, left=left, top=top, width=plot_w, height=plot_h), colour))
        svg.append(f'<line x1="610" y1="{y}" x2="652" y2="{y}" stroke="{colour}" stroke-width="2.5" />')
        svg.append(f'<text x="660" y="{y+4}" font-family="Arial" font-size="12">{escape(name)}</text>')
    svg.append('<text x="450" y="420" text-anchor="middle" font-family="Arial" font-size="13">detuning (GHz)</text>')
    svg.append('<text x="35" y="230" text-anchor="middle" font-family="Arial" font-size="13" transform="rotate(-90 35 230)">normalised value</text>')
    svg.append(_svg_footer())
    path.write_text("\n".join(svg), encoding="utf-8")


def _write_phase_faraday_maps(
    path: Path,
    phase_map: np.ndarray,
    theta_map: np.ndarray,
    label: str,
    kappa_f: float,
) -> None:
    width = 920
    height = 390
    colours = ["#081d58", "#225ea8", "#41b6c4", "#c7e9b4", "#ffffd9"]
    svg = [_svg_header(width, height, "Scalar phase and Faraday rotation maps", label)]
    svg.append('<text x="460" y="32" text-anchor="middle" font-family="Arial" font-size="20">Scalar phase and Faraday rotation maps</text>')
    svg.append(f'<text x="460" y="55" text-anchor="middle" font-family="Arial" font-size="12">Version 1 representative / uncalibrated; kappa_F = {kappa_f:g} placeholder</text>')
    for title, data, x in [
        ("scalar phase map", phase_map, 115),
        ("theta_F map", theta_map, 525),
    ]:
        svg.append(_heatmap_rects(data, x=x, y=88, width=260, height=260, colours=colours))
        svg.append(f'<rect x="{x}" y="88" width="260" height="260" fill="none" stroke="#222" />')
        svg.append(f'<text x="{x+130}" y="372" text-anchor="middle" font-family="Arial" font-size="13">{escape(title)}</text>')
    svg.append(_svg_footer())
    path.write_text("\n".join(svg), encoding="utf-8")


def _write_image_comparison(path: Path, images: dict[str, np.ndarray], label: str) -> None:
    width = 980
    height = 350
    names = list(images)
    colours = ["#000004", "#3b0f70", "#8c2981", "#de4968", "#fe9f6d", "#fcfdbf"]
    svg = [_svg_header(width, height, "Representative imaging-mode comparison", label)]
    svg.append('<text x="490" y="32" text-anchor="middle" font-family="Arial" font-size="20">Representative PCI / DGI / Faraday image comparison</text>')
    svg.append('<text x="490" y="55" text-anchor="middle" font-family="Arial" font-size="12">Same representative cloud; Version 1 uncalibrated; kappa_F = 1.0 placeholder</text>')
    for i, name in enumerate(names):
        x = 45 + i * 235
        svg.append(_heatmap_rects(images[name], x=x, y=85, width=190, height=190, colours=colours))
        svg.append(f'<rect x="{x}" y="85" width="190" height="190" fill="none" stroke="#222" />')
        svg.append(f'<text x="{x+95}" y="305" text-anchor="middle" font-family="Arial" font-size="13">{escape(name)}</text>')
    svg.append(_svg_footer())
    path.write_text("\n".join(svg), encoding="utf-8")


def _write_image_lineouts(path: Path, axis_m: np.ndarray, images: dict[str, np.ndarray], label: str) -> None:
    width = 980
    height = 470
    x_um = axis_m * 1e6
    left, top, plot_w, plot_h = 95, 92, 735, 290
    colours = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd"]
    svg = [_svg_header(width, height, "Imaging-mode centre lineouts", label)]
    svg.append('<text x="490" y="32" text-anchor="middle" font-family="Arial" font-size="20">PCI / DGI / Faraday centre-line comparison</text>')
    svg.append('<text x="490" y="55" text-anchor="middle" font-family="Arial" font-size="12">Normalised image lineouts; Version 1 representative / uncalibrated</text>')
    svg.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="black" />')
    svg.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="black" />')
    for index, (name, data) in enumerate(images.items()):
        line = np.asarray(data)[np.asarray(data).shape[0] // 2]
        values = _normalise(line)
        y = 112 + index * 22
        colour = colours[index % len(colours)]
        svg.append(_polyline(_line_points(x_um, values, left=left, top=top, width=plot_w, height=plot_h), colour, 2.0))
        svg.append(f'<line x1="660" y1="{y}" x2="702" y2="{y}" stroke="{colour}" stroke-width="2.5" />')
        svg.append(f'<text x="710" y="{y+4}" font-family="Arial" font-size="12">{escape(name)}</text>')
    svg.append('<text x="490" y="435" text-anchor="middle" font-family="Arial" font-size="13">position (um)</text>')
    svg.append('<text x="38" y="237" text-anchor="middle" font-family="Arial" font-size="13" transform="rotate(-90 38 237)">normalised centre-line signal</text>')
    svg.append(_svg_footer())
    path.write_text("\n".join(svg), encoding="utf-8")


def _write_camera_realism(path: Path, ideal: np.ndarray, noisy: np.ndarray, label: str) -> None:
    width = 920
    height = 430
    mid = ideal.shape[0] // 2
    x = np.arange(ideal.shape[1])
    colours = ["#000004", "#3b0f70", "#8c2981", "#de4968", "#fe9f6d", "#fcfdbf"]
    svg = [_svg_header(width, height, "Camera realism", label)]
    svg.append('<text x="460" y="32" text-anchor="middle" font-family="Arial" font-size="20">Camera-level measurement realism</text>')
    svg.append('<text x="460" y="55" text-anchor="middle" font-family="Arial" font-size="12">Noiseless binned image versus one fixed-seed noisy frame</text>')
    svg.append(_heatmap_rects(ideal, x=55, y=85, width=190, height=190, colours=colours))
    svg.append(_heatmap_rects(noisy, x=285, y=85, width=190, height=190, colours=colours, vmin=float(np.min(ideal)), vmax=float(np.max(ideal))))
    svg.append('<rect x="55" y="85" width="190" height="190" fill="none" stroke="#222" />')
    svg.append('<rect x="285" y="85" width="190" height="190" fill="none" stroke="#222" />')
    svg.append('<text x="150" y="302" text-anchor="middle" font-family="Arial" font-size="13">Noiseless camera image</text>')
    svg.append('<text x="380" y="302" text-anchor="middle" font-family="Arial" font-size="13">Noisy camera image</text>')
    left, top, plot_w, plot_h = 550, 95, 310, 190
    svg.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="black" />')
    svg.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="black" />')
    svg.append(_polyline(_line_points(x, ideal[mid], left=left, top=top, width=plot_w, height=plot_h), "#1f77b4"))
    svg.append(_polyline(_line_points(x, noisy[mid], left=left, top=top, width=plot_w, height=plot_h), "#d62728", 1.8))
    svg.append('<text x="705" y="322" text-anchor="middle" font-family="Arial" font-size="13">camera pixel</text>')
    svg.append(_svg_footer())
    path.write_text("\n".join(svg), encoding="utf-8")


def _write_optimisation_summary(path: Path, rows: list[dict[str, float]], label: str, metric_key: str) -> None:
    width = 1000
    height = 500
    panels = [
        ("detuning", "detuning (GHz)", 70),
        ("intensity", "probe power (mW)", 365),
        ("exposure", "exposure (us)", 660),
    ]
    svg = [_svg_header(width, height, "Faraday single-variable optimisation summary", label)]
    svg.append('<text x="500" y="32" text-anchor="middle" font-family="Arial" font-size="20">Faraday single-variable optimisation summary</text>')
    svg.append(f'<text x="500" y="55" text-anchor="middle" font-family="Arial" font-size="12">Metric: {escape(metric_key)}; Version 1 representative / uncalibrated</text>')
    for sweep_name, x_label, left in panels:
        panel_rows = [row for row in rows if row["sweep_id"] == sweep_name]
        x_values = np.asarray([row["parameter_display"] for row in panel_rows])
        values = _normalise(np.asarray([row["metric_value"] for row in panel_rows]))
        top, plot_w, plot_h = 102, 245, 260
        svg.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="black" />')
        svg.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="black" />')
        svg.append(_polyline(_line_points(x_values, values, left=left, top=top, width=plot_w, height=plot_h), "#2ca02c", 2.2))
        best_index = int(np.nanargmax(np.asarray([row["metric_value"] for row in panel_rows], dtype=float)))
        x_span = float(np.max(x_values) - np.min(x_values)) or 1.0
        bx = left + ((float(x_values[best_index]) - float(np.min(x_values))) / x_span) * plot_w
        by = top + (1 - float(values[best_index])) * plot_h
        svg.append(f'<circle cx="{bx:.2f}" cy="{by:.2f}" r="4" fill="#d62728" />')
        svg.append(f'<text x="{left + plot_w / 2:.2f}" y="395" text-anchor="middle" font-family="Arial" font-size="12">{escape(x_label)}</text>')
        svg.append(f'<text x="{left + plot_w / 2:.2f}" y="86" text-anchor="middle" font-family="Arial" font-size="13">{escape(sweep_name)}</text>')
    svg.append('<text x="500" y="455" text-anchor="middle" font-family="Arial" font-size="12">Each panel is a separate one-dimensional sweep; this is not a 2D/3D optimisation.</text>')
    svg.append(_svg_footer())
    path.write_text("\n".join(svg), encoding="utf-8")


def _write_multishot(path: Path, rows: list[dict[str, float]], label: str) -> None:
    width = 940
    height = 460
    shot = np.asarray([row["shot"] for row in rows])
    survival = np.asarray([row["condensate_fraction"] for row in rows])
    loss = np.asarray([row["loss_fraction"] for row in rows])
    acc = _normalise(np.asarray([row["accumulated_snr"] for row in rows]))
    left, top, plot_w, plot_h = 90, 90, 720, 280
    svg = [_svg_header(width, height, "Deterministic multi-shot evolution", label)]
    svg.append('<text x="470" y="32" text-anchor="middle" font-family="Arial" font-size="20">Deterministic multi-shot evolution</text>')
    svg.append('<text x="470" y="55" text-anchor="middle" font-family="Arial" font-size="12">Representative sequence bookkeeping; no noisy frame rendering</text>')
    svg.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="black" />')
    svg.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="black" />')
    for values, colour, name, y in [
        (survival, "#1f77b4", "condensate fraction", 106),
        (loss, "#d62728", "loss fraction", 128),
        (acc, "#2ca02c", "normalised accumulated SNR", 150),
    ]:
        svg.append(_polyline(_line_points(shot, values, left=left, top=top, width=plot_w, height=plot_h), colour))
        svg.append(f'<line x1="590" y1="{y}" x2="632" y2="{y}" stroke="{colour}" stroke-width="2.5" />')
        svg.append(f'<text x="640" y="{y+4}" font-family="Arial" font-size="12">{escape(name)}</text>')
    svg.append('<text x="470" y="420" text-anchor="middle" font-family="Arial" font-size="13">frame number</text>')
    svg.append('<text x="35" y="230" text-anchor="middle" font-family="Arial" font-size="13" transform="rotate(-90 35 230)">fraction / normalised value</text>')
    svg.append(_svg_footer())
    path.write_text("\n".join(svg), encoding="utf-8")


def _write_absorption_observables(path: Path, od: np.ndarray, lineout: np.ndarray, moments: dict[str, float], label: str) -> None:
    width = 940
    height = 430
    x = np.arange(lineout.size)
    colours = ["#ffffd9", "#c7e9b4", "#41b6c4", "#225ea8", "#081d58"]
    svg = [_svg_header(width, height, "Synthetic absorption calibration observables", label)]
    svg.append('<text x="470" y="32" text-anchor="middle" font-family="Arial" font-size="20">Absorption / RAI calibration-readiness observables</text>')
    svg.append('<text x="470" y="55" text-anchor="middle" font-family="Arial" font-size="12">Synthetic OD example only; no experimental calibration applied</text>')
    svg.append(_heatmap_rects(od, x=65, y=90, width=230, height=230, colours=colours))
    svg.append('<rect x="65" y="90" width="230" height="230" fill="none" stroke="#222" />')
    svg.append('<text x="180" y="352" text-anchor="middle" font-family="Arial" font-size="13">synthetic optical density</text>')
    left, top, plot_w, plot_h = 375, 105, 300, 190
    svg.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="black" />')
    svg.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="black" />')
    svg.append(_polyline(_line_points(x, lineout, left=left, top=top, width=plot_w, height=plot_h), "#1f77b4", 2.2))
    svg.append('<text x="525" y="328" text-anchor="middle" font-family="Arial" font-size="13">camera column</text>')
    svg.append('<text x="335" y="200" text-anchor="middle" font-family="Arial" font-size="13" transform="rotate(-90 335 200)">OD lineout</text>')
    stats = [
        ("peak OD", moments["peak_od"]),
        ("integrated OD", moments["integrated_od"]),
        ("centre x", moments["centre_x"]),
        ("centre y", moments["centre_y"]),
        ("width x", moments["width_x"]),
        ("width y", moments["width_y"]),
    ]
    svg.append('<rect x="725" y="105" width="165" height="190" rx="6" fill="#f7fbff" stroke="#2d5c88" />')
    svg.append('<text x="807" y="132" text-anchor="middle" font-family="Arial" font-size="13">extracted observables</text>')
    for index, (name, value) in enumerate(stats):
        svg.append(f'<text x="742" y="{160 + index * 22}" font-family="Arial" font-size="12">{escape(name)}: {value:.4g}</text>')
    svg.append(_svg_footer())
    path.write_text("\n".join(svg), encoding="utf-8")


def _write_figure_index(path: Path, records: list[dict[str, Any]]) -> None:
    _write_json(
        path,
        {
            "status": "Version 1 representative / uncalibrated figure index.",
            "figures": records,
        },
    )


def _build_cloud(config: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    cloud = config["cloud"]
    points = int(cloud["grid_points"])
    axis = np.linspace(-cloud["field_of_view_m"] / 2, cloud["field_of_view_m"] / 2, points)
    x_grid, y_grid = np.meshgrid(axis, axis)
    profile = thomas_fermi_profile_2d(x_grid, y_grid, cloud["radius_x_m"], cloud["radius_y_m"])
    density = profile * config["physical_parameters"]["column_density_peak"]
    return axis, profile, density


def _build_pupil(shape: tuple[int, int], radius_fraction: float) -> np.ndarray:
    fy = np.fft.fftfreq(shape[0])
    fx = np.fft.fftfreq(shape[1])
    fx_grid, fy_grid = np.meshgrid(fx, fy)
    return (fx_grid**2 + fy_grid**2 <= radius_fraction**2).astype(float)


def _generate(config: dict[str, Any], config_path: Path) -> dict[str, str]:
    output_dir = Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    label = config["label"]
    physical = config["physical_parameters"]
    imaging = config["imaging"]
    camera = config["camera"]

    axis, profile, density = _build_cloud(config)
    phase_peak = scalar_phase_shift(
        imaging["detuning_hz"],
        physical["column_density_peak"],
        physical["resonant_cross_section"],
        physical["gamma_rad_per_s"],
    )
    phase_map = phase_peak * profile
    theta_map = physical["kappa_F"] * phase_map
    pupil = _build_pupil(profile.shape, imaging["pupil_radius_fraction"])

    paths = {
        "closed_loop_architecture": output_dir / "closed_loop_architecture.svg",
        "tf_density_model": output_dir / "tf_density_model.svg",
        "tf_density_model_data": output_dir / "tf_density_model_data.csv",
        "tf_profile_lineouts": output_dir / "tf_profile_lineouts.svg",
        "tf_profile_lineouts_data": output_dir / "tf_profile_lineouts_data.csv",
        "light_atom_tradeoff": output_dir / "light_atom_tradeoff.svg",
        "light_atom_tradeoff_data": output_dir / "light_atom_tradeoff_data.csv",
        "phase_and_faraday_maps": output_dir / "phase_and_faraday_maps.svg",
        "phase_and_faraday_maps_metadata": output_dir / "phase_and_faraday_maps_metadata.json",
        "imaging_mode_comparison": output_dir / "imaging_mode_comparison.svg",
        "imaging_mode_comparison_metadata": output_dir / "imaging_mode_comparison_metadata.json",
        "imaging_mode_lineouts": output_dir / "imaging_mode_lineouts.svg",
        "imaging_mode_lineouts_data": output_dir / "imaging_mode_lineouts_data.csv",
        "camera_realism": output_dir / "camera_realism.svg",
        "camera_realism_metadata": output_dir / "camera_realism_metadata.json",
        "multishot_evolution": output_dir / "multishot_evolution.svg",
        "multishot_evolution_data": output_dir / "multishot_evolution_data.csv",
        "optimisation_sweep_summary": output_dir / "optimisation_sweep_summary.svg",
        "optimisation_sweep_summary_data": output_dir / "optimisation_sweep_summary_data.csv",
        "optimisation_sweep_summary_metadata": output_dir / "optimisation_sweep_summary_metadata.json",
        "absorption_calibration_observables": output_dir / "absorption_calibration_observables.svg",
        "absorption_calibration_observables_data": output_dir / "absorption_calibration_observables_data.csv",
        "absorption_calibration_observables_metadata": output_dir / "absorption_calibration_observables_metadata.json",
        "figure_index": output_dir / "figure_index.json",
        "metadata": output_dir / "metadata.json",
        "summary": output_dir / "summary.json",
    }

    _write_closed_loop(paths["closed_loop_architecture"], label)
    _write_density_figure(paths["tf_density_model"], axis, density, label)
    mid = density.shape[0] // 2
    _write_rows(
        paths["tf_density_model_data"],
        [
            {
                "position_m": float(axis[index]),
                "position_um": float(axis[index] * 1e6),
                "normalised_centre_line": float(profile[mid, index]),
                "column_density": float(density[mid, index]),
            }
            for index in range(axis.size)
        ],
    )
    _write_profile_lineouts(paths["tf_profile_lineouts"], axis, profile, label)
    _write_rows(
        paths["tf_profile_lineouts_data"],
        [
            {
                "position_m": float(axis[index]),
                "position_um": float(axis[index] * 1e6),
                "x_centre_line_normalised": float(profile[mid, index]),
                "y_centre_line_normalised": float(profile[index, mid]),
            }
            for index in range(axis.size)
        ],
    )

    trade_rows = []
    for detuning in config["light_atom_tradeoff"]["detuning_values_hz"]:
        signal = abs(
            physical["kappa_F"]
            * scalar_phase_shift(
                detuning,
                physical["column_density_peak"],
                physical["resonant_cross_section"],
                physical["gamma_rad_per_s"],
            )
        )
        scattered = scattered_photons_per_atom(
            detuning,
            imaging["probe_power_mw"],
            imaging["pulse_duration_s"],
            physical["saturation_intensity"],
            physical["gamma_rad_per_s"],
            physical["probe_diameter_m"],
            use_peak_intensity=physical.get("use_peak_intensity", True),
        )
        reabs = reabsorption_fraction(
            detuning,
            np.asarray(physical["column_densities_for_reabsorption"]),
            physical["resonant_cross_section"],
            physical["gamma_rad_per_s"],
        )
        trade_rows.append(
            {
                "detuning_hz": float(detuning),
                "detuning_ghz": float(detuning / 1e9),
                "signal_scale": float(signal),
                "scattered_photons_per_atom": float(scattered),
                "reabsorption_fraction": float(reabs),
                "destructiveness_metric": float(scattered * (1 + reabs)),
                "signal_per_scattered_photon": float(signal / scattered),
            }
        )
    _write_rows(paths["light_atom_tradeoff_data"], trade_rows)
    _write_tradeoff_figure(paths["light_atom_tradeoff"], trade_rows, label)
    _write_phase_faraday_maps(paths["phase_and_faraday_maps"], phase_map, theta_map, label, physical["kappa_F"])
    _write_json(
        paths["phase_and_faraday_maps_metadata"],
        {
            "label": label,
            "phase_peak_rad": phase_peak,
            "theta_peak_rad": float(np.max(theta_map)),
            "kappa_F": physical["kappa_F"],
            "note": "Scalar phase and phenomenological Faraday rotation maps are Version 1 representative / uncalibrated.",
        },
    )

    pci = simulate_pci_image(
        phase_map,
        pupil,
        phase_plate_transmittance=imaging["phase_plate_transmittance"],
        phase_plate_phase=imaging["phase_plate_phase_rad"],
    )
    dgi = simulate_dgi_image(
        phase_map,
        pupil,
        stop_optical_depth=imaging["dgi_stop_optical_depth"],
    )
    faraday = simulate_faraday_image(theta_map, pupil)
    images = {
        "PCI intensity": np.asarray(pci),
        "DGI intensity": np.asarray(dgi),
        "Faraday dark field": np.asarray(faraday["dark_field_intensity"]),
        "Faraday dual-port S": np.asarray(faraday["dual_port_signal"]),
    }
    _write_image_comparison(paths["imaging_mode_comparison"], images, label)
    _write_image_lineouts(paths["imaging_mode_lineouts"], axis, images, label)
    _write_rows(
        paths["imaging_mode_lineouts_data"],
        [
            {
                "position_m": float(axis[index]),
                "position_um": float(axis[index] * 1e6),
                "pci_intensity": float(images["PCI intensity"][mid, index]),
                "dgi_intensity": float(images["DGI intensity"][mid, index]),
                "faraday_dark_field": float(images["Faraday dark field"][mid, index]),
                "faraday_dual_port_signal": float(images["Faraday dual-port S"][mid, index]),
            }
            for index in range(axis.size)
        ],
    )
    _write_json(
        paths["imaging_mode_comparison_metadata"],
        {
            "label": label,
            "kappa_F": physical["kappa_F"],
            "phase_peak_rad": phase_peak,
            "theta_peak_rad": float(np.max(theta_map)),
            "image_min_max": {key: [float(np.min(value)), float(np.max(value))] for key, value in images.items()},
            "note": "Version 1 representative / uncalibrated imaging-mode comparison.",
        },
    )

    ideal_camera = simulate_camera_image(
        pci,
        bin_size=camera["bin_size"],
        photons_per_pixel=camera["photons_per_pixel"],
    )
    rng = np.random.default_rng(camera["rng_seed"])
    noisy_camera = simulate_noisy_camera_image(
        pci,
        camera["photons_per_pixel"],
        rng,
        camera["read_noise_electrons"],
        bin_size=camera["bin_size"],
    )
    _write_camera_realism(paths["camera_realism"], np.asarray(ideal_camera), np.asarray(noisy_camera), label)
    _write_json(
        paths["camera_realism_metadata"],
        {
            "label": label,
            "rng_seed": camera["rng_seed"],
            "photons_per_pixel": camera["photons_per_pixel"],
            "read_noise_electrons": camera["read_noise_electrons"],
            "bin_size": camera["bin_size"],
            "note": "One fixed-seed noisy frame for representative measurement realism, not noise averaging.",
        },
    )

    multi = config["multishot"]
    ref_phase = float(multi["reference_phase_rad"])
    ref_atoms = float(multi["initial_condensate_atoms"])
    sequence = simulate_multishot_sequence(
        multi["photons_scattered_per_atom_per_shot"],
        ref_atoms,
        loss_fraction_limit=multi["loss_fraction_limit"],
        max_shots=multi["max_shots"],
        model=multi["model"],
        eta_coll=multi["eta_coll"],
        phase_from_n0=lambda n0: ref_phase * n0 / ref_atoms,
        snr_from_phi=lambda phi: multi["reference_snr"] * abs(phi) / ref_phase,
    )
    accumulated = accumulate_snr(sequence["snr"])
    multi_rows = [
        {
            "shot": float(sequence["shot"][index]),
            "N0": float(sequence["N0"][index]),
            "condensate_fraction": float(sequence["condensate_fraction"][index]),
            "loss_fraction": float(sequence["loss_fraction"][index]),
            "phi": float(sequence["phi"][index]),
            "snr": float(sequence["snr"][index]),
            "accumulated_snr": float(accumulated[index]),
        }
        for index in range(sequence["shot"].size)
    ]
    _write_rows(paths["multishot_evolution_data"], multi_rows)
    _write_multishot(paths["multishot_evolution"], multi_rows, label)

    optimisation = config["optimisation_sweeps"]
    metric_key = optimisation["metric_key"]
    sweep_kwargs = {
        "column_density_peak": physical["column_density_peak"],
        "resonant_cross_section": physical["resonant_cross_section"],
        "gamma_rad_per_s": physical["gamma_rad_per_s"],
        "saturation_intensity": physical["saturation_intensity"],
        "probe_diameter_m": physical["probe_diameter_m"],
        "kappa_f": physical["kappa_F"],
        "column_densities_for_reabsorption": np.asarray(physical["column_densities_for_reabsorption"]),
        "use_peak_intensity": physical.get("use_peak_intensity", True),
        "photons_per_camera_pixel": camera["photons_per_pixel"],
        "objective_key": metric_key,
    }
    detuning_sweep = sweep_faraday_detuning(
        optimisation["detuning_values_hz"],
        probe_power_mw=imaging["probe_power_mw"],
        pulse_duration_s=imaging["pulse_duration_s"],
        **sweep_kwargs,
    )
    intensity_sweep = sweep_faraday_intensity(
        optimisation["probe_power_mw_values"],
        detuning_hz=imaging["detuning_hz"],
        pulse_duration_s=imaging["pulse_duration_s"],
        **sweep_kwargs,
    )
    exposure_sweep = sweep_faraday_exposure_time(
        optimisation["pulse_duration_s_values"],
        detuning_hz=imaging["detuning_hz"],
        probe_power_mw=imaging["probe_power_mw"],
        **sweep_kwargs,
    )
    optimisation_rows: list[dict[str, Any]] = []
    for sweep_id, parameter_key, display_scale, sweep in [
        ("detuning", "detuning_hz", 1e-9, detuning_sweep),
        ("intensity", "probe_power_mw", 1.0, intensity_sweep),
        ("exposure", "pulse_duration_s", 1e6, exposure_sweep),
    ]:
        parameters = np.asarray(sweep[parameter_key], dtype=float)
        metrics = np.asarray(sweep[metric_key], dtype=float)
        signal = np.asarray(sweep["faraday_signal_scale"], dtype=float)
        scattered = np.asarray(sweep["scattered_photons_per_atom"], dtype=float)
        destructiveness = np.asarray(sweep["destructiveness_metric"], dtype=float)
        for index in range(parameters.size):
            optimisation_rows.append(
                {
                    "sweep_id": sweep_id,
                    "parameter_raw": float(parameters[index]),
                    "parameter_display": float(parameters[index] * display_scale),
                    "metric_value": float(metrics[index]),
                    "faraday_signal_scale": float(signal[index]),
                    "scattered_photons_per_atom": float(scattered[index]),
                    "destructiveness_metric": float(destructiveness[index]),
                }
            )
    _write_rows(paths["optimisation_sweep_summary_data"], optimisation_rows)
    _write_optimisation_summary(paths["optimisation_sweep_summary"], optimisation_rows, label, metric_key)
    _write_json(
        paths["optimisation_sweep_summary_metadata"],
        {
            "label": label,
            "metric_key": metric_key,
            "detuning_summary": summarise_faraday_sweep(detuning_sweep, metric_key),
            "intensity_summary": summarise_faraday_sweep(intensity_sweep, metric_key),
            "exposure_summary": summarise_faraday_sweep(exposure_sweep, metric_key),
            "note": "Single-variable deterministic sweeps only; no 2D/3D optimisation or stochastic averaging.",
        },
    )

    absorption = config["absorption_calibration_demo"]
    peak_od = float(absorption["peak_optical_density"])
    probe_counts = float(absorption["probe_counts"])
    dark_counts = float(absorption["dark_counts"])
    od_model = peak_od * profile
    probe_image = np.full_like(od_model, probe_counts)
    dark_image = np.full_like(od_model, dark_counts)
    atom_image = dark_image + (probe_image - dark_image) * np.exp(-od_model)
    od = compute_optical_density(atom_image, probe_image, dark_image, epsilon=float(absorption["epsilon"]))
    od_moments = estimate_cloud_moments(od, x=axis, y=axis)
    od_lineout = od[mid, :]
    _write_rows(
        paths["absorption_calibration_observables_data"],
        [
            {
                "position_m": float(axis[index]),
                "position_um": float(axis[index] * 1e6),
                "synthetic_od_centre_line": float(od_lineout[index]),
            }
            for index in range(axis.size)
        ],
    )
    _write_absorption_observables(paths["absorption_calibration_observables"], od, od_lineout, od_moments, label)
    _write_json(
        paths["absorption_calibration_observables_metadata"],
        {
            "label": label,
            "source": "Synthetic absorption/RAI calibration-readiness example generated from the representative profile.",
            "probe_counts": probe_counts,
            "dark_counts": dark_counts,
            "peak_optical_density": peak_od,
            "observables": od_moments,
            "calibration_status": "No experimental absorption/RAI data or fitted calibration has been applied.",
        },
    )

    figure_records = [
        {
            "filename": path.name,
            "title": title,
            "what_it_shows": description,
            "helpers_or_logic_used": helpers,
            "deterministic": True,
            "depends_on_placeholder_parameters": placeholder,
            "v1_representative_uncalibrated": True,
        }
        for path, title, description, helpers, placeholder in [
            (paths["closed_loop_architecture"], "Closed-loop calibration architecture", "Planned calibration-aware simulation workflow.", "script schematic", True),
            (paths["tf_density_model"], "Thomas-Fermi density model", "Representative column-density object used by the gallery.", "thomas_fermi_profile_2d", False),
            (paths["tf_profile_lineouts"], "Thomas-Fermi profile lineouts", "Centre-line views of the representative cloud object.", "thomas_fermi_profile_2d", False),
            (paths["light_atom_tradeoff"], "Light-atom trade-off", "Signal, scattering, and signal-per-scattered-photon trends versus detuning.", "scalar_phase_shift, scattered_photons_per_atom, reabsorption_fraction", True),
            (paths["phase_and_faraday_maps"], "Phase and Faraday maps", "Scalar phase and theta_F maps for the representative cloud.", "scalar_phase_shift with kappa_F placeholder", True),
            (paths["imaging_mode_comparison"], "Imaging mode comparison", "PCI, DGI, Faraday dark-field, and Faraday dual-port representative images.", "simulate_pci_image, simulate_dgi_image, simulate_faraday_image", True),
            (paths["imaging_mode_lineouts"], "Imaging mode lineouts", "Normalised centre-line signals from the representative imaging modes.", "simulate_pci_image, simulate_dgi_image, simulate_faraday_image", True),
            (paths["camera_realism"], "Camera realism", "Noiseless binned camera image and one fixed-seed noisy camera image.", "simulate_camera_image, simulate_noisy_camera_image", True),
            (paths["multishot_evolution"], "Multi-shot evolution", "Deterministic sequence loss and accumulated SNR bookkeeping.", "simulate_multishot_sequence, accumulate_snr", True),
            (paths["optimisation_sweep_summary"], "Optimisation sweep summary", "One-dimensional detuning, intensity, and exposure-time objective trends.", "sweep_faraday_detuning, sweep_faraday_intensity, sweep_faraday_exposure_time, summarise_faraday_sweep", True),
            (paths["absorption_calibration_observables"], "Absorption calibration observables", "Synthetic OD image and extracted calibration-readiness observables.", "compute_optical_density, estimate_cloud_moments", True),
        ]
    ]
    _write_figure_index(paths["figure_index"], figure_records)

    summary = {
        "label": label,
        "figures": {
            "closed_loop_architecture": str(paths["closed_loop_architecture"]),
            "tf_density_model": str(paths["tf_density_model"]),
            "tf_profile_lineouts": str(paths["tf_profile_lineouts"]),
            "light_atom_tradeoff": str(paths["light_atom_tradeoff"]),
            "phase_and_faraday_maps": str(paths["phase_and_faraday_maps"]),
            "imaging_mode_comparison": str(paths["imaging_mode_comparison"]),
            "imaging_mode_lineouts": str(paths["imaging_mode_lineouts"]),
            "camera_realism": str(paths["camera_realism"]),
            "multishot_evolution": str(paths["multishot_evolution"]),
            "optimisation_sweep_summary": str(paths["optimisation_sweep_summary"]),
            "absorption_calibration_observables": str(paths["absorption_calibration_observables"]),
        },
        "main_text_candidates": [
            "closed_loop_architecture.svg",
            "tf_density_model.svg",
            "tf_profile_lineouts.svg",
            "light_atom_tradeoff.svg",
            "phase_and_faraday_maps.svg",
            "imaging_mode_comparison.svg",
            "imaging_mode_lineouts.svg",
            "camera_realism.svg",
            "multishot_evolution.svg",
            "optimisation_sweep_summary.svg",
            "absorption_calibration_observables.svg",
        ],
        "index": str(paths["figure_index"]),
        "status": "Version 1 representative / uncalibrated figure gallery.",
    }
    _write_json(paths["summary"], summary)

    metadata = {
        "label": label,
        "status": "Version 1 representative / uncalibrated figures",
        "config_file_used": str(config_path),
        "git_commit_hash": _git_commit(),
        "helper_functions_used": HELPERS_USED,
        "generated_figure_filenames": [str(path) for key, path in paths.items() if str(path).endswith(".svg")],
        "kappa_F": physical["kappa_F"],
        "kappa_F_note": "kappa_F remains the Version 1 placeholder and has not been experimentally fitted.",
        "rng_seed": camera["rng_seed"],
        "validation_status": "pytest and notebook section validation are expected to be run by the caller.",
        "calibration_status": "No experimental RAI / absorption calibration has been applied.",
        "overclaiming_boundary": "These outputs are figure placeholders, not final experimental predictions.",
        "outputs": {key: str(path) for key, path in paths.items()},
    }
    _write_json(paths["metadata"], metadata)

    return {key: str(path) for key, path in paths.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()

    config = _load_config(args.config)
    outputs = _generate(config, args.config)
    print("Generated dissertation figure gallery outputs:")
    for key, path in outputs.items():
        print(f"- {key}: {path}")


if __name__ == "__main__":
    main()
