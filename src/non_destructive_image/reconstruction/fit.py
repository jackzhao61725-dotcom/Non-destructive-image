"""Multi-start frozen-weight IRLS fits for raw Faraday count channels."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import least_squares

from .diagnostics import FitDiagnostics, FaradayMeasurement, SmoothTFFitResult, covariance_to_correlation
from .noise import gaussian_quasi_deviance, poisson_gaussian_variance
from .parameters import (
    PARAMETER_NAMES,
    SmoothTFBounds,
    SmoothTFParameters,
    covariance_to_physical,
    from_internal,
    to_internal,
)


@dataclass(frozen=True)
class FitOptions:
    """Numerical controls for a smooth-TF raw-channel fit."""

    irls_iterations: int = 2
    max_nfev: int = 60
    xtol: float = 1e-8
    ftol: float = 1e-8
    gtol: float = 1e-8

    def __post_init__(self) -> None:
        if self.irls_iterations <= 0:
            raise ValueError("irls_iterations must be positive")
        if self.max_nfev <= 0:
            raise ValueError("max_nfev must be positive")
        if min(self.xtol, self.ftol, self.gtol) <= 0:
            raise ValueError("fit tolerances must be positive")


def _active_bounds(
    vector: NDArray[np.floating],
    lower: NDArray[np.floating],
    upper: NDArray[np.floating],
) -> dict[str, bool]:
    tolerance = 1e-5 * np.maximum(np.abs(upper - lower), 1.0)
    return {
        name: bool(abs(value - lo) < tol or abs(hi - value) < tol)
        for name, value, lo, hi, tol in zip(
            PARAMETER_NAMES,
            vector,
            lower,
            upper,
            tolerance,
            strict=True,
        )
    }


def fit_smooth_tf(
    measurement: FaradayMeasurement,
    observed_channels: tuple[ArrayLike, ...],
    *,
    starts: Iterable[SmoothTFParameters],
    bounds: SmoothTFBounds,
    options: FitOptions | None = None,
) -> SmoothTFFitResult:
    """Fit a smooth-TF object to raw DPFI or DFFI camera channels.

    Count variances are frozen within each least-squares solve and updated
    between IRLS passes. The full Gaussian quasi-deviance ranks multi-start
    solutions. This approximation is intentionally explicit so it can later be
    benchmarked against a mixed Poisson-Gaussian likelihood.
    """

    fit_options = options or FitOptions()
    observed = measurement.flatten_observed(*observed_channels)
    lower = to_internal(bounds.lower)
    upper = to_internal(bounds.upper)
    starts_list = list(starts)
    if not starts_list:
        raise ValueError("at least one initial parameter set is required")

    best: tuple[float, object, NDArray[np.floating], int] | None = None
    for start in starts_list:
        current = to_internal(start)
        final_result = None
        completed_irls = 0
        for outer in range(fit_options.irls_iterations):
            initial_prediction, _ = measurement.expected_vector_and_jacobian_internal(current)
            standard_deviation = np.sqrt(
                poisson_gaussian_variance(
                    initial_prediction,
                    measurement.read_noise_electrons,
                )
            )
            cached_vector: NDArray[np.floating] | None = None
            cached_prediction: NDArray[np.floating] | None = None
            cached_jacobian: NDArray[np.floating] | None = None

            def evaluate(
                vector: NDArray[np.floating],
            ) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
                nonlocal cached_vector, cached_prediction, cached_jacobian
                if cached_vector is None or not np.array_equal(vector, cached_vector):
                    cached_prediction, cached_jacobian = (
                        measurement.expected_vector_and_jacobian_internal(vector)
                    )
                    cached_vector = np.array(vector, copy=True)
                if cached_prediction is None or cached_jacobian is None:
                    raise RuntimeError("prediction cache was not populated")
                return cached_prediction, cached_jacobian

            final_result = least_squares(
                lambda vector: (observed - evaluate(vector)[0]) / standard_deviation,
                current,
                jac=lambda vector: -evaluate(vector)[1] / standard_deviation[:, None],
                bounds=(lower, upper),
                x_scale="jac",
                max_nfev=fit_options.max_nfev,
                xtol=fit_options.xtol,
                ftol=fit_options.ftol,
                gtol=fit_options.gtol,
            )
            current = np.asarray(final_result.x, dtype=float)
            completed_irls = outer + 1

        if final_result is None:
            raise RuntimeError("smooth-TF fit did not execute")
        final_prediction, _ = measurement.expected_vector_and_jacobian_internal(current)
        objective = gaussian_quasi_deviance(
            observed,
            final_prediction,
            read_noise_electrons=measurement.read_noise_electrons,
        )
        candidate = (objective, final_result, current, completed_irls)
        if best is None or objective < best[0]:
            best = candidate

    if best is None:
        raise RuntimeError("smooth-TF fit did not produce a candidate")

    objective, result, vector, completed_irls = best
    parameters = from_internal(vector)
    prediction, jacobian_physical_coordinates = (
        measurement.expected_vector_and_jacobian_internal(vector)
    )
    variance = poisson_gaussian_variance(prediction, measurement.read_noise_electrons)
    standard_deviation = np.sqrt(variance)
    whitened = (observed - prediction) / standard_deviation
    weighted_jacobian = jacobian_physical_coordinates / standard_deviation[:, None]
    singular_values = np.linalg.svd(weighted_jacobian, compute_uv=False)
    threshold = (
        np.finfo(float).eps * max(weighted_jacobian.shape) * singular_values[0]
        if singular_values.size
        else float("inf")
    )
    rank = int(np.count_nonzero(singular_values > threshold))
    condition = (
        float(singular_values[0] / singular_values[-1])
        if singular_values.size and singular_values[-1] > threshold
        else float("inf")
    )
    covariance_internal = np.linalg.pinv(weighted_jacobian.T @ weighted_jacobian, rcond=1e-12)
    covariance_physical = covariance_to_physical(parameters, covariance_internal)
    standard_errors = {
        name: float(value)
        for name, value in zip(
            PARAMETER_NAMES,
            np.sqrt(np.clip(np.diag(covariance_physical), 0.0, None)),
            strict=True,
        )
    }
    degrees_of_freedom = max(observed.size - len(PARAMETER_NAMES), 1)
    chi_square = float(np.sum(whitened**2))
    diagnostics = FitDiagnostics(
        success=bool(result.success and np.all(np.isfinite(vector))),
        message=str(result.message),
        quasi_deviance=float(objective),
        weighted_chi_square=chi_square,
        reduced_chi_square=chi_square / degrees_of_freedom,
        degrees_of_freedom=degrees_of_freedom,
        nfev=int(result.nfev),
        irls_iterations=completed_irls,
        jacobian_rank=rank,
        jacobian_condition=condition,
        local_standard_errors=standard_errors,
        local_correlation=covariance_to_correlation(covariance_physical),
        active_bounds=_active_bounds(vector, lower, upper),
        whitened_residual_vector=whitened,
    )
    return SmoothTFFitResult(
        parameters=parameters,
        predicted_channels=measurement.expected_channels(parameters),
        diagnostics=diagnostics,
    )
