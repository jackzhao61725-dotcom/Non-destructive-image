"""Non-negative basis reconstruction from raw Faraday count channels."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import least_squares

from .measurements import DifferentiableDensityMeasurement
from .noise import gaussian_quasi_deviance, poisson_gaussian_variance
from .object_models import NonnegativeBilinearDensityModel
from .regularisation import CurvatureRegularisation


@dataclass(frozen=True)
class DensityFitOptions:
    """Numerical controls for the shape-flexible density reconstruction."""

    irls_iterations: int = 2
    max_nfev: int = 120
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


@dataclass(frozen=True)
class DensityFitDiagnostics:
    """Data support and numerical diagnostics for one density reconstruction."""

    success: bool
    message: str
    quasi_deviance: float
    weighted_chi_square: float
    reduced_chi_square: float
    degrees_of_freedom: int
    nfev: int
    irls_iterations: int
    data_jacobian_rank: int
    data_jacobian_condition: float
    curvature_weight_um2: float
    curvature_seminorm_per_um2: float
    regularisation_objective: float
    regularisation_density_scale_m2: float | None
    regularisation_boundary_policy: str | None
    regularisation_axis_weights: tuple[float, float, float] | None
    active_lower_coefficients: int
    active_upper_coefficients: int
    whitened_residual_vector: NDArray[np.floating]


@dataclass(frozen=True)
class DensityFitResult:
    """Recovered non-negative density, raw-channel prediction and diagnostics."""

    coefficients: NDArray[np.floating]
    column_density_m2: NDArray[np.floating]
    predicted_channels: tuple[NDArray[np.floating], ...]
    diagnostics: DensityFitDiagnostics


def _coefficient_upper_bound(
    coefficient_upper: float | ArrayLike,
    parameter_count: int,
) -> NDArray[np.floating]:
    upper = np.asarray(coefficient_upper, dtype=float)
    if upper.ndim == 0:
        upper = np.full(parameter_count, float(upper), dtype=float)
    if upper.shape != (parameter_count,):
        raise ValueError(f"coefficient upper bound must have shape ({parameter_count},)")
    if np.any(~np.isfinite(upper)) or np.any(upper <= 0):
        raise ValueError("coefficient upper bounds must be finite and positive")
    return upper


def fit_nonnegative_basis_density(
    measurement: DifferentiableDensityMeasurement,
    model: NonnegativeBilinearDensityModel,
    observed_channels: tuple[ArrayLike, ...],
    *,
    initial_coefficients: ArrayLike,
    coefficient_upper: float | ArrayLike,
    regularisation: CurvatureRegularisation | None,
    options: DensityFitOptions | None = None,
) -> DensityFitResult:
    """Fit a finite-support density without imposing a TF-like shape.

    The raw camera channels provide the data term.  Non-negativity and the
    fixed basis support are hard constraints.  An optional physically scaled
    curvature operator supplies the regularisation residual.  The rank
    diagnostic is calculated from the data Jacobian alone, so regularisation
    cannot masquerade as information supplied by the measurement.
    """

    fit_options = options or DensityFitOptions()
    observed = measurement.flatten_observed(*observed_channels)
    lower = np.zeros(model.parameter_count, dtype=float)
    upper = _coefficient_upper_bound(coefficient_upper, model.parameter_count)
    current = np.asarray(initial_coefficients, dtype=float)
    if current.shape != (model.parameter_count,):
        raise ValueError(
            f"initial coefficient vector must have shape ({model.parameter_count},)"
        )
    if np.any(~np.isfinite(current)) or np.any(current < lower) or np.any(current > upper):
        raise ValueError("initial coefficients must be finite and lie inside the bounds")

    if regularisation is not None:
        if regularisation.parameter_count != model.parameter_count:
            raise ValueError("regularisation parameter count does not match the model")
        if not np.array_equal(regularisation.knot_y_um, model.knot_y_um) or not np.array_equal(
            regularisation.knot_z_um,
            model.knot_z_um,
        ):
            raise ValueError("regularisation knots do not match the density model")
        penalty = regularisation.matrix_for_coefficient_scale(
            model.coefficient_scale_m2
        )
    else:
        penalty = np.zeros((0, model.parameter_count), dtype=float)
    final_result = None
    completed_irls = 0
    for outer in range(fit_options.irls_iterations):
        initial_prediction, _ = measurement.expected_vector_and_jacobian_model(model, current)
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
                    measurement.expected_vector_and_jacobian_model(model, vector)
                )
                cached_vector = np.array(vector, copy=True)
            if cached_prediction is None or cached_jacobian is None:
                raise RuntimeError("density prediction cache was not populated")
            return cached_prediction, cached_jacobian

        def residual(vector: NDArray[np.floating]) -> NDArray[np.floating]:
            data = (observed - evaluate(vector)[0]) / standard_deviation
            return np.concatenate([data, penalty @ vector]) if penalty.size else data

        def jacobian(vector: NDArray[np.floating]) -> NDArray[np.floating]:
            data = -evaluate(vector)[1] / standard_deviation[:, None]
            return np.vstack([data, penalty]) if penalty.size else data

        final_result = least_squares(
            residual,
            current,
            jac=jacobian,
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
        raise RuntimeError("density fit did not execute")

    prediction, data_jacobian = measurement.expected_vector_and_jacobian_model(model, current)
    variance = poisson_gaussian_variance(prediction, measurement.read_noise_electrons)
    standard_deviation = np.sqrt(variance)
    whitened = (observed - prediction) / standard_deviation
    weighted_data_jacobian = data_jacobian / standard_deviation[:, None]
    singular_values = np.linalg.svd(weighted_data_jacobian, compute_uv=False)
    threshold = (
        np.finfo(float).eps * max(weighted_data_jacobian.shape) * singular_values[0]
        if singular_values.size
        else float("inf")
    )
    rank = int(np.count_nonzero(singular_values > threshold))
    condition = (
        float(singular_values[0] / singular_values[-1])
        if singular_values.size
        and rank == model.parameter_count
        and singular_values[-1] > threshold
        else float("inf")
    )
    degrees_of_freedom = max(observed.size - rank, 1)
    weighted_chi_square = float(np.sum(whitened**2))
    regularisation_residual = penalty @ current
    regularisation_objective = 0.5 * float(
        regularisation_residual @ regularisation_residual
    )
    curvature_weight = 0.0 if regularisation is None else regularisation.weight_um2
    curvature_seminorm = (
        float(regularisation_residual @ regularisation_residual) / curvature_weight
        if curvature_weight > 0.0
        else 0.0
    )
    tolerance = 1e-6 * np.maximum(upper, 1.0)
    density = model.column_density(current)
    diagnostics = DensityFitDiagnostics(
        success=bool(final_result.success and np.all(np.isfinite(current))),
        message=str(final_result.message),
        quasi_deviance=gaussian_quasi_deviance(
            observed,
            prediction,
            read_noise_electrons=measurement.read_noise_electrons,
        ),
        weighted_chi_square=weighted_chi_square,
        reduced_chi_square=weighted_chi_square / degrees_of_freedom,
        degrees_of_freedom=degrees_of_freedom,
        nfev=int(final_result.nfev),
        irls_iterations=completed_irls,
        data_jacobian_rank=rank,
        data_jacobian_condition=condition,
        curvature_weight_um2=float(curvature_weight),
        curvature_seminorm_per_um2=curvature_seminorm,
        regularisation_objective=regularisation_objective,
        regularisation_density_scale_m2=(
            None if regularisation is None else regularisation.density_scale_m2
        ),
        regularisation_boundary_policy=(
            None if regularisation is None else regularisation.boundary_policy
        ),
        regularisation_axis_weights=(
            None
            if regularisation is None
            else (
                regularisation.axis_weights.yy,
                regularisation.axis_weights.yz,
                regularisation.axis_weights.zz,
            )
        ),
        active_lower_coefficients=int(np.count_nonzero(current <= tolerance)),
        active_upper_coefficients=int(np.count_nonzero(upper - current <= tolerance)),
        whitened_residual_vector=whitened,
    )
    return DensityFitResult(
        coefficients=current,
        column_density_m2=density,
        predicted_channels=measurement.expected_channels_from_density(density),
        diagnostics=diagnostics,
    )
