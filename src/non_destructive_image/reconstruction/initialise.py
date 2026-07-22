"""Blind, modality-specific initialisation for smooth-TF fits."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
from numpy.typing import ArrayLike
from scipy.ndimage import gaussian_filter

from .measurements import DarkFieldFaradayMeasurement, DualPortFaradayMeasurement
from .parameters import SmoothTFBounds, SmoothTFParameters


def _clip_inside(value: float, lower: float, upper: float) -> float:
    margin = 1e-8 * max(abs(lower), abs(upper), 1.0)
    return float(np.clip(value, lower + margin, upper - margin))


def _fallback(bounds: SmoothTFBounds) -> SmoothTFParameters:
    return SmoothTFParameters(
        column_density_peak_m2=np.sqrt(
            bounds.lower.column_density_peak_m2 * bounds.upper.column_density_peak_m2
        ),
        y0_um=0.5 * (bounds.lower.y0_um + bounds.upper.y0_um),
        z0_um=0.5 * (bounds.lower.z0_um + bounds.upper.z0_um),
        radius_y_um=np.sqrt(bounds.lower.radius_y_um * bounds.upper.radius_y_um),
        radius_z_um=np.sqrt(bounds.lower.radius_z_um * bounds.upper.radius_z_um),
    )


def _moments_to_parameters(
    *,
    weights: np.ndarray,
    camera_y_um: np.ndarray,
    camera_z_um: np.ndarray,
    peak_rotation: float,
    rotation_per_column_density: float,
    bounds: SmoothTFBounds,
) -> SmoothTFParameters:
    weight_sum = float(np.sum(weights))
    if weight_sum <= np.finfo(float).eps or rotation_per_column_density == 0:
        return _fallback(bounds)
    y_grid, z_grid = np.meshgrid(camera_y_um, camera_z_um)
    y0 = float(np.sum(weights * y_grid) / weight_sum)
    z0 = float(np.sum(weights * z_grid) / weight_sum)
    variance_y = float(np.sum(weights * (y_grid - y0) ** 2) / weight_sum)
    variance_z = float(np.sum(weights * (z_grid - z0) ** 2) / weight_sum)
    # For the projected M0 profile, <y^2>=R_y^2/7 and <z^2>=R_z^2/7.
    radius_y = np.sqrt(max(7.0 * variance_y, np.finfo(float).eps))
    radius_z = np.sqrt(max(7.0 * variance_z, np.finfo(float).eps))
    peak_density = abs(peak_rotation / rotation_per_column_density)
    return SmoothTFParameters(
        column_density_peak_m2=_clip_inside(
            peak_density,
            bounds.lower.column_density_peak_m2,
            bounds.upper.column_density_peak_m2,
        ),
        y0_um=_clip_inside(y0, bounds.lower.y0_um, bounds.upper.y0_um),
        z0_um=_clip_inside(z0, bounds.lower.z0_um, bounds.upper.z0_um),
        radius_y_um=_clip_inside(
            radius_y,
            bounds.lower.radius_y_um,
            bounds.upper.radius_y_um,
        ),
        radius_z_um=_clip_inside(
            radius_z,
            bounds.lower.radius_z_um,
            bounds.upper.radius_z_um,
        ),
    )


def estimate_dual_port_initial_parameters(
    measurement: DualPortFaradayMeasurement,
    observed_h_counts: ArrayLike,
    observed_v_counts: ArrayLike,
    bounds: SmoothTFBounds,
) -> SmoothTFParameters:
    """Initialise from a smoothed raw-port normalised difference."""

    h = np.asarray(observed_h_counts, dtype=float)
    v = np.asarray(observed_v_counts, dtype=float)
    if h.shape != measurement.camera_shape or v.shape != measurement.camera_shape:
        raise ValueError(f"dual-port observations must have shape {measurement.camera_shape}")
    denominator = h + v
    signal = np.divide(
        h - v,
        denominator,
        out=np.zeros_like(denominator),
        where=np.abs(denominator) > np.finfo(float).eps,
    )
    smoothed = gaussian_filter(signal, sigma=1.0)
    outside = ~measurement.grid.roi_mask
    background = float(np.median(smoothed[outside])) if np.any(outside) else 0.0
    response = measurement.response.rotation_per_column_density_rad_m2
    response_sign = float(np.sign(response))
    if response_sign == 0:
        return _fallback(bounds)
    oriented_signal = np.clip(response_sign * (smoothed - background), 0.0, None)
    weights = oriented_signal * measurement.grid.roi_mask
    # The ideal weak-angle dual-port signal is approximately 2 theta_F.
    peak_rotation = response_sign * 0.5 * max(float(np.max(oriented_signal)), 0.01)
    return _moments_to_parameters(
        weights=weights,
        camera_y_um=measurement.grid.camera_y_um,
        camera_z_um=measurement.grid.camera_z_um,
        peak_rotation=peak_rotation,
        rotation_per_column_density=measurement.response.rotation_per_column_density_rad_m2,
        bounds=bounds,
    )


def estimate_dark_field_initial_parameters(
    measurement: DarkFieldFaradayMeasurement,
    observed_counts: ArrayLike,
    bounds: SmoothTFBounds,
) -> SmoothTFParameters:
    """Initialise from the square root of background-subtracted intensity."""

    observed = np.asarray(observed_counts, dtype=float)
    if observed.shape != measurement.camera_shape:
        raise ValueError(f"dark-field observations must have shape {measurement.camera_shape}")
    normalised = observed / measurement.detector.photoelectrons_per_i0_pixel
    outside = ~measurement.grid.roi_mask
    background = float(np.median(normalised[outside])) if np.any(outside) else 0.0
    amplitude = np.sqrt(np.clip(normalised - background, 0.0, None))
    smoothed = gaussian_filter(amplitude, sigma=1.0)
    weights = np.clip(smoothed, 0.0, None) * measurement.grid.roi_mask
    peak_rotation = max(float(np.max(smoothed)), 0.005)
    return _moments_to_parameters(
        weights=weights,
        camera_y_um=measurement.grid.camera_y_um,
        camera_z_um=measurement.grid.camera_z_um,
        peak_rotation=peak_rotation,
        rotation_per_column_density=measurement.response.rotation_per_column_density_rad_m2,
        bounds=bounds,
    )


def make_multistart_parameters(
    initial: SmoothTFParameters,
    bounds: SmoothTFBounds,
    *,
    maximum_starts: int = 4,
) -> list[SmoothTFParameters]:
    """Return deterministic starts around a blind estimate."""

    if maximum_starts <= 0:
        raise ValueError("maximum_starts must be positive")
    variants: Iterable[tuple[float, float, float, float, float]] = (
        (1.0, 0.0, 0.0, 1.0, 1.0),
        (0.65, 0.0, 0.0, 0.8, 0.8),
        (1.45, 0.0, 0.0, 1.2, 1.4),
        (1.0, 0.75, -0.50, 1.1, 0.7),
        (1.0, -0.75, 0.50, 0.9, 1.3),
    )
    starts: list[SmoothTFParameters] = []
    for density_scale, y_shift, z_shift, ry_scale, rz_scale in list(variants)[:maximum_starts]:
        starts.append(
            SmoothTFParameters(
                column_density_peak_m2=_clip_inside(
                    initial.column_density_peak_m2 * density_scale,
                    bounds.lower.column_density_peak_m2,
                    bounds.upper.column_density_peak_m2,
                ),
                y0_um=_clip_inside(
                    initial.y0_um + y_shift,
                    bounds.lower.y0_um,
                    bounds.upper.y0_um,
                ),
                z0_um=_clip_inside(
                    initial.z0_um + z_shift,
                    bounds.lower.z0_um,
                    bounds.upper.z0_um,
                ),
                radius_y_um=_clip_inside(
                    initial.radius_y_um * ry_scale,
                    bounds.lower.radius_y_um,
                    bounds.upper.radius_y_um,
                ),
                radius_z_um=_clip_inside(
                    initial.radius_z_um * rz_scale,
                    bounds.lower.radius_z_um,
                    bounds.upper.radius_z_um,
                ),
            )
        )
    return starts
