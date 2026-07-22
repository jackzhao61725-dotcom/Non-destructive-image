"""Truth-independent initialisation for shape-flexible density fits."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import lsq_linear
from scipy.sparse.linalg import LinearOperator

from .initialise import estimate_dark_field_initial_parameters
from .measurements import DarkFieldFaradayMeasurement, DifferentiableDensityMeasurement
from .noise import poisson_gaussian_variance
from .object_models import (
    DifferentiableColumnDensityModel,
    NonnegativeBilinearDensityModel,
    smooth_tf_column_density,
)
from .parameters import SmoothTFBounds


@dataclass(frozen=True)
class DensityInitialisation:
    """Linearised non-negative estimate and its local data support."""

    coefficients: NDArray[np.floating]
    method: str
    data_jacobian_rank: int
    data_jacobian_condition: float
    weighted_residual_norm: float


def _upper_vector(value: float | ArrayLike, parameter_count: int) -> NDArray[np.floating]:
    upper = np.asarray(value, dtype=float)
    if upper.ndim == 0:
        upper = np.full(parameter_count, float(upper), dtype=float)
    if upper.shape != (parameter_count,):
        raise ValueError(f"coefficient upper bound must have shape ({parameter_count},)")
    if np.any(~np.isfinite(upper)) or np.any(upper <= 0):
        raise ValueError("coefficient upper bounds must be finite and positive")
    return upper


def _rank_condition(
    weighted_jacobian: NDArray[np.floating],
    parameter_count: int,
) -> tuple[int, float]:
    singular_values = np.linalg.svd(weighted_jacobian, compute_uv=False)
    threshold = (
        np.finfo(float).eps * max(weighted_jacobian.shape) * singular_values[0]
        if singular_values.size
        else float("inf")
    )
    rank = int(np.count_nonzero(singular_values > threshold))
    condition = (
        float(singular_values[0] / singular_values[-1])
        if rank == parameter_count and singular_values[-1] > threshold
        else float("inf")
    )
    return rank, condition


def project_density_to_nonnegative_basis(
    model: NonnegativeBilinearDensityModel,
    target_column_density_m2: ArrayLike,
    *,
    coefficient_upper: float | ArrayLike,
    ridge_strength: float = 0.0,
) -> NDArray[np.floating]:
    """Project a supplied density estimate onto the finite-support basis."""

    if not np.isfinite(ridge_strength) or ridge_strength < 0:
        raise ValueError("projection ridge strength must be finite and non-negative")
    target = np.asarray(target_column_density_m2, dtype=float)
    if target.shape != model.y_grid_m.shape:
        raise ValueError("projection target must match the density-model grid")
    if np.any(~np.isfinite(target)) or np.any(target < 0):
        raise ValueError("projection target must be finite and non-negative")
    target_scaled = target.ravel() / model.coefficient_scale_m2
    pixel_count = target_scaled.size
    if ridge_strength > 0:
        root_ridge = np.sqrt(ridge_strength)
        design = LinearOperator(
            shape=(pixel_count + model.parameter_count, model.parameter_count),
            matvec=lambda vector: np.concatenate(
                [model.basis_matvec(vector), root_ridge * vector]
            ),
            rmatvec=lambda vector: (
                model.basis_rmatvec(vector[:pixel_count])
                + root_ridge * vector[pixel_count:]
            ),
            dtype=float,
        )
        target_scaled = np.concatenate(
            [target_scaled, np.zeros(model.parameter_count, dtype=float)]
        )
    else:
        design = LinearOperator(
            shape=(pixel_count, model.parameter_count),
            matvec=model.basis_matvec,
            rmatvec=model.basis_rmatvec,
            dtype=float,
        )
    upper = _upper_vector(coefficient_upper, model.parameter_count)
    result = lsq_linear(
        design,
        target_scaled,
        bounds=(np.zeros(model.parameter_count), upper),
        method="trf",
        lsq_solver="lsmr",
        lsmr_tol="auto",
        max_iter=200,
    )
    if not result.success:
        raise RuntimeError(f"density-to-basis projection failed: {result.message}")
    return np.asarray(result.x, dtype=float)


def linearised_nonnegative_initialisation(
    measurement: DifferentiableDensityMeasurement,
    model: DifferentiableColumnDensityModel,
    observed_channels: tuple[ArrayLike, ...],
    *,
    coefficient_upper: float | ArrayLike,
    ridge_strength: float = 0.0,
) -> DensityInitialisation:
    """Estimate coefficients from the instrument Jacobian at zero density.

    This is blind to synthetic truth and uses the same raw channels, pupil and
    detector variance as the final fit.  It is suitable for a locally linear
    readout such as dual-port Faraday imaging.  A quadratic dark-field signal
    has a zero Jacobian at zero density and is rejected explicitly instead of
    returning an arbitrary seed.
    """

    if not np.isfinite(ridge_strength) or ridge_strength < 0:
        raise ValueError("ridge strength must be finite and non-negative")
    observed = measurement.flatten_observed(*observed_channels)
    zero = np.zeros(model.parameter_count, dtype=float)
    baseline, jacobian = measurement.expected_vector_and_jacobian_model(model, zero)
    standard_deviation = np.sqrt(
        poisson_gaussian_variance(baseline, measurement.read_noise_electrons)
    )
    weighted_jacobian = jacobian / standard_deviation[:, None]
    weighted_signal = (observed - baseline) / standard_deviation
    rank, condition = _rank_condition(weighted_jacobian, model.parameter_count)
    if rank == 0:
        raise ValueError(
            "the measurement has no first-order density information at zero; "
            "use a modality-specific nonlinear initialiser"
        )
    data_jacobian = weighted_jacobian
    data_signal = weighted_signal
    if ridge_strength > 0:
        weighted_jacobian = np.vstack(
            [weighted_jacobian, np.sqrt(ridge_strength) * np.eye(model.parameter_count)]
        )
        weighted_signal = np.concatenate(
            [weighted_signal, np.zeros(model.parameter_count, dtype=float)]
        )
    upper = _upper_vector(coefficient_upper, model.parameter_count)
    result = lsq_linear(
        weighted_jacobian,
        weighted_signal,
        bounds=(np.zeros(model.parameter_count), upper),
        method="trf",
        lsmr_tol="auto",
        max_iter=200,
    )
    if not result.success:
        raise RuntimeError(f"linearised density initialisation failed: {result.message}")
    return DensityInitialisation(
        coefficients=np.asarray(result.x, dtype=float),
        method="zero_density_linearisation",
        data_jacobian_rank=rank,
        data_jacobian_condition=condition,
        weighted_residual_norm=float(
            np.linalg.norm(data_jacobian @ result.x - data_signal)
        ),
    )


def dark_field_sqrt_moment_initialisation(
    measurement: DarkFieldFaradayMeasurement,
    model: NonnegativeBilinearDensityModel,
    observed_counts: ArrayLike,
    *,
    smooth_bounds: SmoothTFBounds,
    coefficient_upper: float | ArrayLike,
    projection_ridge_strength: float = 0.0,
) -> DensityInitialisation:
    """Build a blind dark-field basis seed from square-root image moments.

    The square root of the background-subtracted intensity supplies an
    amplitude-like camera map.  Its centre, widths and scale define only a
    smooth starting envelope, which is projected onto the flexible basis.  The
    subsequent nonlinear fit is not restricted to that envelope.
    """

    smooth_seed = estimate_dark_field_initial_parameters(
        measurement,
        observed_counts,
        smooth_bounds,
    )
    seed_density = smooth_tf_column_density(
        measurement.grid.y_grid_m,
        measurement.grid.z_grid_m,
        smooth_seed,
    )
    coefficients = project_density_to_nonnegative_basis(
        model,
        seed_density,
        coefficient_upper=coefficient_upper,
        ridge_strength=projection_ridge_strength,
    )
    observed = measurement.flatten_observed(observed_counts)
    prediction, jacobian = measurement.expected_vector_and_jacobian_model(
        model,
        coefficients,
    )
    standard_deviation = np.sqrt(
        poisson_gaussian_variance(prediction, measurement.read_noise_electrons)
    )
    weighted_jacobian = jacobian / standard_deviation[:, None]
    rank, condition = _rank_condition(weighted_jacobian, model.parameter_count)
    return DensityInitialisation(
        coefficients=coefficients,
        method="dark_field_sqrt_moment_envelope",
        data_jacobian_rank=rank,
        data_jacobian_condition=condition,
        weighted_residual_norm=float(np.linalg.norm((observed - prediction) / standard_deviation)),
    )
