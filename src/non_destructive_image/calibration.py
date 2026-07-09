"""Deterministic absorption-image calibration readiness helpers."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


def _as_matching_array(value: ArrayLike, reference_shape: tuple[int, ...], name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.shape != reference_shape:
        raise ValueError(f"{name} must have shape {reference_shape}, got {array.shape}")
    return array


def compute_optical_density(
    atom_image: ArrayLike,
    probe_image: ArrayLike,
    dark_image: ArrayLike | None = None,
    *,
    epsilon: float = 1e-12,
) -> np.ndarray:
    """Compute absorption optical density from atom, probe, and dark images."""

    if epsilon <= 0:
        raise ValueError("epsilon must be positive")

    atom = np.asarray(atom_image, dtype=float)
    probe = _as_matching_array(probe_image, atom.shape, "probe_image")
    dark = np.zeros_like(atom) if dark_image is None else _as_matching_array(dark_image, atom.shape, "dark_image")

    atom_corrected = np.clip(atom - dark, epsilon, None)
    probe_corrected = np.clip(probe - dark, epsilon, None)
    ratio = np.clip(atom_corrected / probe_corrected, epsilon, None)
    return -np.log(ratio)


def integrate_optical_density(od: ArrayLike, pixel_area: float | None = None) -> float:
    """Integrate an optical-density map, optionally scaling by pixel area."""

    od_array = np.asarray(od, dtype=float)
    if od_array.size == 0:
        raise ValueError("od must be non-empty")
    if not np.all(np.isfinite(od_array)):
        raise ValueError("od must contain only finite values")

    integral = float(np.sum(od_array))
    if pixel_area is not None:
        if pixel_area <= 0:
            raise ValueError("pixel_area must be positive")
        integral *= float(pixel_area)
    return integral


def _coordinate_grids(
    shape: tuple[int, int],
    x: ArrayLike | None,
    y: ArrayLike | None,
) -> tuple[np.ndarray, np.ndarray]:
    if len(shape) != 2:
        raise ValueError("od must be a 2D array")

    rows, columns = shape
    if x is None:
        x_grid = np.broadcast_to(np.arange(columns, dtype=float), shape)
    else:
        x_array = np.asarray(x, dtype=float)
        if x_array.ndim == 1:
            if x_array.size != columns:
                raise ValueError("x must have one value per OD column")
            x_grid = np.broadcast_to(x_array, shape)
        elif x_array.shape == shape:
            x_grid = x_array
        else:
            raise ValueError(f"x must be 1D with length {columns} or 2D with shape {shape}")

    if y is None:
        y_grid = np.broadcast_to(np.arange(rows, dtype=float)[:, np.newaxis], shape)
    else:
        y_array = np.asarray(y, dtype=float)
        if y_array.ndim == 1:
            if y_array.size != rows:
                raise ValueError("y must have one value per OD row")
            y_grid = np.broadcast_to(y_array[:, np.newaxis], shape)
        elif y_array.shape == shape:
            y_grid = y_array
        else:
            raise ValueError(f"y must be 1D with length {rows} or 2D with shape {shape}")

    return x_grid, y_grid


def estimate_cloud_moments(
    od: ArrayLike,
    x: ArrayLike | None = None,
    y: ArrayLike | None = None,
) -> dict[str, float]:
    """Estimate centre, RMS widths, peak OD, and integrated OD from an OD map."""

    od_array = np.asarray(od, dtype=float)
    if od_array.size == 0:
        raise ValueError("od must be non-empty")
    if not np.all(np.isfinite(od_array)):
        raise ValueError("od must contain only finite values")

    x_grid, y_grid = _coordinate_grids(od_array.shape, x, y)
    weights = np.clip(od_array, 0.0, None)
    total_weight = float(np.sum(weights))
    if total_weight <= 0:
        raise ValueError("od must contain positive weight for moment estimation")

    centre_x = float(np.sum(weights * x_grid) / total_weight)
    centre_y = float(np.sum(weights * y_grid) / total_weight)
    width_x = float(np.sqrt(np.sum(weights * (x_grid - centre_x) ** 2) / total_weight))
    width_y = float(np.sqrt(np.sum(weights * (y_grid - centre_y) ** 2) / total_weight))

    return {
        "peak_od": float(np.max(od_array)),
        "integrated_od": integrate_optical_density(od_array),
        "centre_x": centre_x,
        "centre_y": centre_y,
        "width_x": width_x,
        "width_y": width_y,
    }


def extract_absorption_observables(
    atom_image: ArrayLike,
    probe_image: ArrayLike,
    dark_image: ArrayLike | None = None,
    *,
    x: ArrayLike | None = None,
    y: ArrayLike | None = None,
    pixel_area: float | None = None,
    epsilon: float = 1e-12,
) -> dict[str, float]:
    """Compute OD and extract deterministic absorption-image observables."""

    od = compute_optical_density(atom_image, probe_image, dark_image, epsilon=epsilon)
    observables = estimate_cloud_moments(od, x=x, y=y)
    if pixel_area is not None:
        observables["integrated_od"] = integrate_optical_density(od, pixel_area=pixel_area)
    return observables
