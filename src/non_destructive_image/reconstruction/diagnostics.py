"""Fit diagnostics and simulation-only recovery metrics."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .measurements import DarkFieldFaradayMeasurement, DualPortFaradayMeasurement
from .parameters import PARAMETER_NAMES, SmoothTFParameters


FaradayMeasurement = DualPortFaradayMeasurement | DarkFieldFaradayMeasurement


@dataclass(frozen=True)
class FitDiagnostics:
    """Numerical and local-identifiability diagnostics from one fit."""

    success: bool
    message: str
    quasi_deviance: float
    weighted_chi_square: float
    reduced_chi_square: float
    degrees_of_freedom: int
    nfev: int
    irls_iterations: int
    jacobian_rank: int
    jacobian_condition: float
    local_standard_errors: dict[str, float]
    local_correlation: NDArray[np.floating]
    active_bounds: dict[str, bool]
    whitened_residual_vector: NDArray[np.floating]


@dataclass(frozen=True)
class SmoothTFFitResult:
    """Physical estimate, predicted raw channels and fit diagnostics."""

    parameters: SmoothTFParameters
    predicted_channels: tuple[NDArray[np.floating], ...]
    diagnostics: FitDiagnostics


def covariance_to_correlation(covariance: NDArray[np.floating]) -> NDArray[np.floating]:
    """Convert a covariance matrix to a finite correlation matrix."""

    standard = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
    denominator = np.outer(standard, standard)
    correlation = np.divide(
        covariance,
        denominator,
        out=np.zeros_like(covariance, dtype=float),
        where=denominator > np.finfo(float).eps,
    )
    np.fill_diagonal(correlation, np.where(standard > 0, 1.0, 0.0))
    return correlation


def synthetic_recovery_metrics(
    measurement: FaradayMeasurement,
    truth: SmoothTFParameters,
    fitted: SmoothTFParameters,
) -> dict[str, float]:
    """Compare a fit with known synthetic truth.

    These metrics are deliberately isolated from production fit construction so
    that experimental fitting never requires access to a truth object.
    """

    truth_map = measurement.column_density(truth)
    fitted_map = measurement.column_density(fitted)
    y_min = float(np.min(measurement.grid.camera_y_um[measurement.grid.roi_mask.any(axis=0)]))
    y_max = float(np.max(measurement.grid.camera_y_um[measurement.grid.roi_mask.any(axis=0)]))
    z_min = float(np.min(measurement.grid.camera_z_um[measurement.grid.roi_mask.any(axis=1)]))
    z_max = float(np.max(measurement.grid.camera_z_um[measurement.grid.roi_mask.any(axis=1)]))
    high_resolution_roi = (
        (measurement.grid.y_grid_m * 1e6 >= y_min)
        & (measurement.grid.y_grid_m * 1e6 <= y_max)
        & (measurement.grid.z_grid_m * 1e6 >= z_min)
        & (measurement.grid.z_grid_m * 1e6 <= z_max)
    )
    truth_values = truth_map[high_resolution_roi]
    fitted_values = fitted_map[high_resolution_roi]
    denominator = max(float(np.linalg.norm(truth_values)), np.finfo(float).eps)
    integrated_truth = float(np.sum(truth_values))
    integrated_fitted = float(np.sum(fitted_values))
    return {
        "map_relative_l2_error": float(
            np.linalg.norm(fitted_values - truth_values) / denominator
        ),
        "integrated_column_density_relative_error": (
            integrated_fitted / integrated_truth - 1.0
            if integrated_truth != 0
            else float("nan")
        ),
        "column_density_peak_relative_error": (
            fitted.column_density_peak_m2 / truth.column_density_peak_m2 - 1.0
        ),
        "y0_error_um": fitted.y0_um - truth.y0_um,
        "z0_error_um": fitted.z0_um - truth.z0_um,
        "radius_y_relative_error": fitted.radius_y_um / truth.radius_y_um - 1.0,
        "radius_z_relative_error": fitted.radius_z_um / truth.radius_z_um - 1.0,
    }


def strongest_local_parameter_correlation(result: SmoothTFFitResult) -> tuple[str, str, float]:
    """Return the strongest off-diagonal local correlation for reporting."""

    correlation = np.asarray(result.diagnostics.local_correlation, dtype=float)
    off_diagonal = np.abs(correlation - np.eye(correlation.shape[0]))
    row, col = np.unravel_index(np.argmax(off_diagonal), off_diagonal.shape)
    return PARAMETER_NAMES[row], PARAMETER_NAMES[col], float(correlation[row, col])
