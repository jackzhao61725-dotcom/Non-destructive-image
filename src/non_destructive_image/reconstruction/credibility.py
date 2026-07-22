"""Data-support and uncertainty diagnostics for density reconstruction.

This module deliberately does not compare a reconstruction with a truth map.
It instead asks which fitted directions are constrained by the raw camera
counts, how strongly the declared regulariser contributes, whether the
forward residual is consistent with the detector model, and how detector
noise propagates into density and low-order observables.

The calculations are local to one fitted solution.  They are therefore a
credibility layer, not a proof that the optical forward model is complete.
Instrument mismatch, unmodelled backgrounds and alternative regularisation
choices must be assessed separately.
"""

from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Mapping
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .density_fit import (
    DensityFitOptions,
    DensityFitResult,
    fit_nonnegative_basis_density,
)
from .measurements import DifferentiableDensityMeasurement
from .noise import poisson_gaussian_variance, simulate_poisson_gaussian_counts
from .object_models import NonnegativeBilinearDensityModel
from .observables import (
    DensityObservableSummary,
    ObservableIntegrationSupport,
    extract_density_observables,
)
from .regularisation import CurvatureRegularisation


OBSERVABLE_PARAMETER_NAMES: tuple[str, ...] = (
    "integrated_response",
    "centroid_y_m",
    "centroid_z_m",
    "major_rms_width_m",
)
ObservableIntervalStatus = Literal["complete", "partial", "unsupported"]


@dataclass(frozen=True)
class ResidualChannelDiagnostics:
    """Standardised residual map and simple spatial checks for one channel."""

    channel_name: str
    standardised_residual_map: NDArray[np.floating]
    residual_mean: float
    residual_rms: float
    residual_standard_deviation: float
    lag_one_correlation_y: float
    lag_one_correlation_z: float


@dataclass(frozen=True)
class DensityFeatureSummary:
    """Model-independent low-order observables of a column-density map."""

    integrated_column_density: float
    peak_column_density_m2: float
    centroid_y_um: float
    centroid_z_um: float
    rms_y_um: float
    rms_z_um: float


@dataclass(frozen=True)
class FeatureInterval:
    """Point estimate and descriptive conditional-bootstrap spread."""

    estimate: float
    bootstrap_mean: float
    bootstrap_standard_deviation: float
    bootstrap_bias: float
    lower: float
    upper: float


@dataclass(frozen=True)
class ObservableBootstrapInterval:
    """Conditional interval status for one reported physical observable.

    Bounds are reported only when the point estimate and every requested
    bootstrap draw converge and support the observable.  A partial result
    retains its sample count and joint samples but does not silently condition
    an interval on only the finite draws.
    """

    estimate: float | None
    supported_draws: int
    successful_draws: int
    requested_draws: int
    status: ObservableIntervalStatus
    bootstrap_mean: float | None
    bootstrap_standard_deviation: float | None
    bootstrap_bias: float | None
    lower: float | None
    upper: float | None

    def __post_init__(self) -> None:
        if self.requested_draws <= 0:
            raise ValueError("requested_draws must be positive")
        if not 0 < self.successful_draws <= self.requested_draws:
            raise ValueError("successful_draws must lie within requested_draws")
        if not 0 <= self.supported_draws <= self.successful_draws:
            raise ValueError("supported_draws must lie within successful_draws")
        if self.estimate is not None and not np.isfinite(self.estimate):
            raise ValueError("observable point estimate must be finite or None")
        if self.status == "complete":
            if (
                self.estimate is None
                or self.supported_draws != self.successful_draws
                or self.successful_draws != self.requested_draws
            ):
                raise ValueError(
                    "complete interval requires every requested draw to support it"
                )
            reported = (
                self.bootstrap_mean,
                self.bootstrap_standard_deviation,
                self.bootstrap_bias,
                self.lower,
                self.upper,
            )
            if any(value is None or not np.isfinite(value) for value in reported):
                raise ValueError("complete interval statistics must be finite")
            assert self.bootstrap_standard_deviation is not None
            assert self.lower is not None
            assert self.upper is not None
            if self.bootstrap_standard_deviation < 0.0:
                raise ValueError("bootstrap standard deviation cannot be negative")
            if self.lower > self.upper:
                raise ValueError("bootstrap interval lower bound exceeds upper bound")
        else:
            if any(
                value is not None
                for value in (
                    self.bootstrap_mean,
                    self.bootstrap_standard_deviation,
                    self.bootstrap_bias,
                    self.lower,
                    self.upper,
                )
            ):
                raise ValueError(
                    "incomplete intervals cannot report conditional bounds"
                )
            if self.status == "partial" and (
                self.estimate is None
                or self.supported_draws == 0
                or (
                    self.supported_draws == self.successful_draws
                    and self.successful_draws == self.requested_draws
                )
            ):
                raise ValueError(
                    "partial interval requires missing or unsupported requested draws"
                )
            if self.status == "unsupported" and not (
                self.estimate is None or self.supported_draws == 0
            ):
                raise ValueError(
                    "unsupported interval must lack a point or sample support"
                )


@dataclass(frozen=True, eq=False)
class ObservableBootstrapSummary:
    """Joint conditional-bootstrap samples of the paper-facing observables.

    Rows correspond one-for-one to successful refits.  Unsupported moments
    remain in place as ``NaN`` with a false entry in ``supported_mask`` so the
    joint sampling distribution is never altered by dropping draws.
    """

    parameter_names: tuple[str, ...]
    requested_draws: int
    successful_draws: int
    point_estimate: DensityObservableSummary
    samples: NDArray[np.floating]
    supported_mask: NDArray[np.bool_]
    intervals: Mapping[str, ObservableBootstrapInterval]
    joint_supported_draw_count: int
    joint_covariance: NDArray[np.floating] | None

    def __post_init__(self) -> None:
        names = tuple(self.parameter_names)
        if names != OBSERVABLE_PARAMETER_NAMES:
            raise ValueError("observable bootstrap parameter order is not canonical")
        if not isinstance(self.point_estimate, DensityObservableSummary):
            raise TypeError("point_estimate must be DensityObservableSummary")
        if self.requested_draws <= 0:
            raise ValueError("requested_draws must be positive")
        if not 0 < self.successful_draws <= self.requested_draws:
            raise ValueError("successful_draws must lie within requested_draws")
        samples = np.array(self.samples, dtype=float, copy=True, order="C")
        supported = np.array(self.supported_mask, dtype=bool, copy=True, order="C")
        if samples.ndim != 2 or samples.shape[1] != len(names):
            raise ValueError("observable samples must have one canonical row per draw")
        if samples.shape[0] != self.successful_draws:
            raise ValueError("observable sample rows must equal successful_draws")
        if supported.shape != samples.shape:
            raise ValueError("observable supported mask must match the sample matrix")
        if not np.array_equal(supported, np.isfinite(samples)):
            raise ValueError(
                "observable support mask must identify finite samples exactly"
            )
        intervals = dict(self.intervals)
        if set(intervals) != set(names):
            raise ValueError("observable intervals must cover the canonical parameters")
        for index, name in enumerate(names):
            interval = intervals[name]
            if not isinstance(interval, ObservableBootstrapInterval):
                raise TypeError("observable intervals have the wrong value type")
            if (
                interval.requested_draws != self.requested_draws
                or interval.successful_draws != self.successful_draws
                or (
                    interval.supported_draws
                    != int(np.count_nonzero(supported[:, index]))
                )
            ):
                raise ValueError("observable interval counts do not match the samples")
        joint_count = int(np.count_nonzero(np.all(supported, axis=1)))
        if self.joint_supported_draw_count != joint_count:
            raise ValueError("joint supported draw count is inconsistent with samples")
        covariance: NDArray[np.floating] | None = None
        if self.joint_covariance is not None:
            covariance = np.array(
                self.joint_covariance,
                dtype=float,
                copy=True,
                order="C",
            )
            if covariance.shape != (len(names), len(names)):
                raise ValueError("joint covariance has the wrong shape")
            if (
                joint_count != self.requested_draws
                or self.successful_draws != self.requested_draws
                or samples.shape[0] < 2
            ):
                raise ValueError(
                    "joint covariance requires full support in at least two draws"
                )
            if np.any(~np.isfinite(covariance)):
                raise ValueError("joint covariance must be finite")
            if not np.allclose(covariance, covariance.T, rtol=1e-12, atol=0.0):
                raise ValueError("joint covariance must be symmetric")
        elif (
            joint_count == self.requested_draws
            and self.successful_draws == self.requested_draws
            and samples.shape[0] >= 2
        ):
            raise ValueError("fully supported samples require a joint covariance")
        samples.setflags(write=False)
        supported.setflags(write=False)
        if covariance is not None:
            covariance.setflags(write=False)
        object.__setattr__(self, "parameter_names", names)
        object.__setattr__(self, "samples", samples)
        object.__setattr__(self, "supported_mask", supported)
        object.__setattr__(self, "intervals", intervals)
        object.__setattr__(self, "joint_covariance", covariance)


@dataclass(frozen=True)
class LocalCredibilityDiagnostics:
    """Local data/prior support and conditional noise propagation.

    Active-bound coefficients are excluded from the local linear calculation.
    Their entries in coefficient-space arrays are ``NaN``.  The density
    uncertainty is conditional on the fitted active set, the declared
    measurement operator and the selected regulariser.
    """

    free_coefficient_mask: NDArray[np.bool_]
    data_singular_values: NDArray[np.floating]
    relative_data_singular_values: NDArray[np.floating]
    generalised_data_mode_fractions: NDArray[np.floating]
    combined_constrained_rank: int
    effective_data_degrees_of_freedom: float
    effective_prior_degrees_of_freedom: float
    locally_unconstrained_degrees_of_freedom: float
    coefficient_resolution_diagonal: NDArray[np.floating]
    coefficient_noise_standard_uncertainty: NDArray[np.floating]
    conditional_noise_covariance: NDArray[np.floating]
    density_standard_uncertainty_m2: NDArray[np.floating] | None
    residual_channels: tuple[ResidualChannelDiagnostics, ...]
    residual_port_cross_correlation: float | None
    assumptions: tuple[str, ...]


@dataclass(frozen=True)
class ParametricBootstrapResult:
    """Conditional detector-noise bootstrap around one fitted solution.

    The legacy coefficient, density-image and feature fields are populated by
    default for frozen-study compatibility.  Parameter-first callers can set
    ``retain_latent_artifacts=False`` and use ``observables`` exclusively.
    """

    requested_draws: int
    successful_draws: int
    seed: int
    confidence_level: float
    coefficient_samples: NDArray[np.floating] | None
    density_mean_m2: NDArray[np.floating] | None
    density_standard_uncertainty_m2: NDArray[np.floating] | None
    density_interval_lower_m2: NDArray[np.floating] | None
    density_interval_upper_m2: NDArray[np.floating] | None
    feature_intervals: Mapping[str, FeatureInterval]
    observables: ObservableBootstrapSummary | None
    assumptions: tuple[str, ...]


@dataclass(frozen=True)
class ReconstructionStabilitySummary:
    """Spread across explicitly different reconstruction choices."""

    variant_labels: tuple[str, ...]
    density_mean_m2: NDArray[np.floating]
    density_standard_deviation_m2: NDArray[np.floating]
    density_minimum_m2: NDArray[np.floating]
    density_maximum_m2: NDArray[np.floating]
    feature_ranges: Mapping[str, tuple[float, float]]


def _upper_bound_vector(
    coefficient_upper: float | ArrayLike,
    parameter_count: int,
) -> NDArray[np.floating]:
    upper = np.asarray(coefficient_upper, dtype=float)
    if upper.ndim == 0:
        upper = np.full(parameter_count, float(upper), dtype=float)
    if upper.shape != (parameter_count,):
        raise ValueError(
            f"coefficient upper bound must have shape ({parameter_count},)"
        )
    if np.any(~np.isfinite(upper)) or np.any(upper <= 0.0):
        raise ValueError("coefficient upper bounds must be finite and positive")
    return upper


def _regularisation_matrix(
    regularisation: CurvatureRegularisation | None,
    model: NonnegativeBilinearDensityModel,
) -> NDArray[np.floating]:
    if regularisation is None:
        return np.zeros((0, model.parameter_count), dtype=float)
    if regularisation.parameter_count != model.parameter_count:
        raise ValueError("regularisation parameter count does not match the model")
    if not np.array_equal(regularisation.knot_y_um, model.knot_y_um):
        raise ValueError("regularisation y knots do not match the density model")
    if not np.array_equal(regularisation.knot_z_um, model.knot_z_um):
        raise ValueError("regularisation z knots do not match the density model")
    return regularisation.matrix_for_coefficient_scale(model.coefficient_scale_m2)


def _symmetric_pseudoinverse_and_inverse_sqrt(
    matrix: NDArray[np.floating],
) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("matrix must be square")
    if matrix.size == 0:
        return matrix.copy(), matrix.copy()
    symmetrised = 0.5 * (matrix + matrix.T)
    eigenvalues, eigenvectors = np.linalg.eigh(symmetrised)
    largest = max(float(np.max(eigenvalues)), 0.0)
    threshold = np.finfo(float).eps * max(matrix.shape) * largest
    retained = eigenvalues > threshold
    inverse_values = np.zeros_like(eigenvalues)
    inverse_sqrt_values = np.zeros_like(eigenvalues)
    inverse_values[retained] = 1.0 / eigenvalues[retained]
    inverse_sqrt_values[retained] = 1.0 / np.sqrt(eigenvalues[retained])
    inverse = (eigenvectors * inverse_values) @ eigenvectors.T
    inverse_sqrt = (eigenvectors * inverse_sqrt_values) @ eigenvectors.T
    return inverse, inverse_sqrt


def _masked_lag_one_correlation(
    values: NDArray[np.floating],
    mask: NDArray[np.bool_],
    *,
    axis: int,
) -> float:
    if axis == 0:
        first, second = values[:-1, :], values[1:, :]
        valid = mask[:-1, :] & mask[1:, :]
    elif axis == 1:
        first, second = values[:, :-1], values[:, 1:]
        valid = mask[:, :-1] & mask[:, 1:]
    else:
        raise ValueError("axis must be zero or one")
    if np.count_nonzero(valid) < 2:
        return float("nan")
    a = first[valid]
    b = second[valid]
    a = a - float(np.mean(a))
    b = b - float(np.mean(b))
    denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(a @ b / denominator) if denominator > 0.0 else float("nan")


def _residual_diagnostics(
    measurement: DifferentiableDensityMeasurement,
    observed_channels: tuple[ArrayLike, ...],
    predicted_channels: tuple[NDArray[np.floating], ...],
) -> tuple[tuple[ResidualChannelDiagnostics, ...], float | None]:
    if len(observed_channels) != len(predicted_channels):
        raise ValueError("observed and predicted channel counts differ")
    names = tuple(
        str(name)
        for name in getattr(
            measurement,
            "channel_names",
            tuple(f"channel_{index}" for index in range(len(predicted_channels))),
        )
    )
    if len(names) != len(predicted_channels):
        raise ValueError("measurement channel names do not match its outputs")
    diagnostics: list[ResidualChannelDiagnostics] = []
    roi_vectors: list[NDArray[np.floating]] = []
    roi = measurement.grid.roi_mask
    for name, observed_raw, predicted_raw in zip(
        names,
        observed_channels,
        predicted_channels,
        strict=True,
    ):
        observed = np.asarray(observed_raw, dtype=float)
        predicted = np.asarray(predicted_raw, dtype=float)
        if observed.shape != predicted.shape or observed.shape != measurement.grid.camera_shape:
            raise ValueError("camera channels must share the declared camera shape")
        standard = np.sqrt(
            poisson_gaussian_variance(
                predicted,
                measurement.read_noise_electrons,
            )
        )
        residual = (observed - predicted) / standard
        residual_map = np.full(residual.shape, np.nan, dtype=float)
        residual_map[roi] = residual[roi]
        values = residual[roi]
        roi_vectors.append(values)
        diagnostics.append(
            ResidualChannelDiagnostics(
                channel_name=name,
                standardised_residual_map=residual_map,
                residual_mean=float(np.mean(values)),
                residual_rms=float(np.sqrt(np.mean(values**2))),
                residual_standard_deviation=float(np.std(values, ddof=1))
                if values.size > 1
                else 0.0,
                lag_one_correlation_y=_masked_lag_one_correlation(
                    residual,
                    roi,
                    axis=1,
                ),
                lag_one_correlation_z=_masked_lag_one_correlation(
                    residual,
                    roi,
                    axis=0,
                ),
            )
        )
    cross_correlation: float | None = None
    if len(roi_vectors) == 2:
        first = roi_vectors[0] - float(np.mean(roi_vectors[0]))
        second = roi_vectors[1] - float(np.mean(roi_vectors[1]))
        denominator = float(np.linalg.norm(first) * np.linalg.norm(second))
        cross_correlation = (
            float(first @ second / denominator)
            if denominator > 0.0
            else float("nan")
        )
    return tuple(diagnostics), cross_correlation


def summarise_density_features(
    column_density_m2: ArrayLike,
    measurement: DifferentiableDensityMeasurement,
) -> DensityFeatureSummary:
    """Return low-order observables without assuming a morphology family."""

    density = np.asarray(column_density_m2, dtype=float)
    grid = measurement.grid
    if density.shape != grid.y_grid_m.shape:
        raise ValueError("density map does not match the reconstruction grid")
    if np.any(~np.isfinite(density)) or np.any(density < 0.0):
        raise ValueError("density map must be finite and non-negative")
    y_axis = np.asarray(grid.y_grid_m[0], dtype=float)
    z_axis = np.asarray(grid.z_grid_m[:, 0], dtype=float)
    if y_axis.size < 2 or z_axis.size < 2:
        raise ValueError("density moments require at least two samples per axis")
    dy = float(np.median(np.diff(y_axis)))
    dz = float(np.median(np.diff(z_axis)))
    if dy <= 0.0 or dz <= 0.0:
        raise ValueError("reconstruction coordinates must increase along each axis")
    raw_total = float(np.sum(density))
    if raw_total <= 0.0:
        return DensityFeatureSummary(0.0, 0.0, *(float("nan"),) * 4)
    y_um = grid.y_grid_m * 1e6
    z_um = grid.z_grid_m * 1e6
    centroid_y = float(np.sum(density * y_um) / raw_total)
    centroid_z = float(np.sum(density * z_um) / raw_total)
    rms_y = float(np.sqrt(np.sum(density * (y_um - centroid_y) ** 2) / raw_total))
    rms_z = float(np.sqrt(np.sum(density * (z_um - centroid_z) ** 2) / raw_total))
    return DensityFeatureSummary(
        integrated_column_density=raw_total * dy * dz,
        peak_column_density_m2=float(np.max(density)),
        centroid_y_um=centroid_y,
        centroid_z_um=centroid_z,
        rms_y_um=rms_y,
        rms_z_um=rms_z,
    )


def analyse_local_credibility(
    measurement: DifferentiableDensityMeasurement,
    model: NonnegativeBilinearDensityModel,
    observed_channels: tuple[ArrayLike, ...],
    fit_result: DensityFitResult,
    *,
    coefficient_upper: float | ArrayLike,
    regularisation: CurvatureRegularisation | None,
    active_tolerance: float = 1e-6,
    retain_density_artifacts: bool = True,
) -> LocalCredibilityDiagnostics:
    """Quantify local data support without access to a truth density.

    ``generalised_data_mode_fractions`` are eigenvalues of the data Hessian in
    the metric of the combined data-plus-regularisation Hessian.  Values near
    one identify locally data-dominated modes; values near zero identify modes
    stabilised mainly by the declared regulariser.  They are diagnostics, not
    posterior probabilities.
    """

    if active_tolerance <= 0.0:
        raise ValueError("active tolerance must be positive")
    if not isinstance(retain_density_artifacts, bool):
        raise TypeError("retain_density_artifacts must be boolean")
    coefficients = np.asarray(fit_result.coefficients, dtype=float)
    if coefficients.shape != (model.parameter_count,):
        raise ValueError("fit coefficients do not match the density model")
    upper = _upper_bound_vector(coefficient_upper, model.parameter_count)
    tolerance = active_tolerance * np.maximum(upper, 1.0)
    active = (coefficients <= tolerance) | (upper - coefficients <= tolerance)
    free = ~active
    if not np.any(free):
        raise ValueError("no free coefficients remain for local credibility analysis")

    observed = measurement.flatten_observed(*observed_channels)
    prediction, jacobian = measurement.expected_vector_and_jacobian_model(
        model,
        coefficients,
    )
    if observed.shape != prediction.shape:
        raise ValueError("observed and predicted vectors do not match")
    standard = np.sqrt(
        poisson_gaussian_variance(prediction, measurement.read_noise_electrons)
    )
    weighted_jacobian = jacobian / standard[:, None]
    free_jacobian = weighted_jacobian[:, free]
    singular_values = np.linalg.svd(free_jacobian, compute_uv=False)
    relative_singular = (
        singular_values / singular_values[0]
        if singular_values.size and singular_values[0] > 0.0
        else np.zeros_like(singular_values)
    )

    penalty = _regularisation_matrix(regularisation, model)
    free_penalty = penalty[:, free]
    data_hessian = free_jacobian.T @ free_jacobian
    prior_hessian = free_penalty.T @ free_penalty
    combined_hessian = data_hessian + prior_hessian
    combined_inverse, combined_inverse_sqrt = (
        _symmetric_pseudoinverse_and_inverse_sqrt(combined_hessian)
    )
    support_metric = combined_inverse_sqrt @ data_hessian @ combined_inverse_sqrt
    support_metric = 0.5 * (support_metric + support_metric.T)
    mode_fractions = np.linalg.eigvalsh(support_metric)
    mode_fractions = np.clip(mode_fractions, 0.0, 1.0)[::-1]
    resolution_free = combined_inverse @ data_hessian
    effective_data = float(np.trace(resolution_free))
    prior_resolution_free = combined_inverse @ prior_hessian
    effective_prior = float(np.trace(prior_resolution_free))
    combined_projector = combined_inverse @ combined_hessian
    combined_rank = int(round(float(np.trace(combined_projector))))
    free_count = int(np.count_nonzero(free))
    effective_unconstrained = float(
        max(free_count - effective_data - effective_prior, 0.0)
    )
    noise_covariance_free = combined_inverse @ data_hessian @ combined_inverse
    noise_covariance_free = 0.5 * (
        noise_covariance_free + noise_covariance_free.T
    )

    full_covariance = np.full(
        (model.parameter_count, model.parameter_count),
        np.nan,
        dtype=float,
    )
    full_covariance[np.ix_(free, free)] = noise_covariance_free
    coefficient_standard = np.full(model.parameter_count, np.nan, dtype=float)
    coefficient_standard[free] = np.sqrt(
        np.clip(np.diag(noise_covariance_free), 0.0, None)
    )
    resolution_diagonal = np.full(model.parameter_count, np.nan, dtype=float)
    resolution_diagonal[free] = np.diag(resolution_free)

    density_standard: NDArray[np.floating] | None = None
    if retain_density_artifacts:
        _, density_derivatives = model.column_density_and_jacobian(coefficients)
        free_density_jacobian = density_derivatives[free].reshape(free_count, -1)
        density_variance = np.einsum(
            "ip,ij,jp->p",
            free_density_jacobian,
            noise_covariance_free,
            free_density_jacobian,
            optimize=True,
        )
        density_standard = np.sqrt(np.clip(density_variance, 0.0, None)).reshape(
            model.y_grid_m.shape
        )

    predicted_channels = measurement.expected_channels_from_density(
        fit_result.column_density_m2
    )
    residual_channels, port_cross = _residual_diagnostics(
        measurement,
        observed_channels,
        predicted_channels,
    )
    return LocalCredibilityDiagnostics(
        free_coefficient_mask=free,
        data_singular_values=singular_values,
        relative_data_singular_values=relative_singular,
        generalised_data_mode_fractions=mode_fractions,
        combined_constrained_rank=combined_rank,
        effective_data_degrees_of_freedom=effective_data,
        effective_prior_degrees_of_freedom=effective_prior,
        locally_unconstrained_degrees_of_freedom=effective_unconstrained,
        coefficient_resolution_diagonal=resolution_diagonal,
        coefficient_noise_standard_uncertainty=coefficient_standard,
        conditional_noise_covariance=full_covariance,
        density_standard_uncertainty_m2=density_standard,
        residual_channels=residual_channels,
        residual_port_cross_correlation=port_cross,
        assumptions=(
            "local linearisation at the fitted solution",
            "fixed active coefficient set",
            "declared Poisson-Gaussian variance model",
            "declared forward operator and curvature regulariser",
            "does not include instrument mismatch or morphology-model uncertainty",
        ),
    )


def _feature_mapping(summary: DensityFeatureSummary) -> dict[str, float]:
    return {
        "integrated_column_density": summary.integrated_column_density,
        "peak_column_density_m2": summary.peak_column_density_m2,
        "centroid_y_um": summary.centroid_y_um,
        "centroid_z_um": summary.centroid_z_um,
        "rms_y_um": summary.rms_y_um,
        "rms_z_um": summary.rms_z_um,
    }


def _observable_vector(summary: DensityObservableSummary) -> NDArray[np.floating]:
    """Return the canonical SI-unit observable vector, preserving missing moments."""

    return np.asarray(
        [
            summary.integrated_response,
            np.nan if summary.centroid_y_m is None else summary.centroid_y_m,
            np.nan if summary.centroid_z_m is None else summary.centroid_z_m,
            (
                np.nan
                if summary.major_rms_width_m is None
                else summary.major_rms_width_m
            ),
        ],
        dtype=float,
    )


def _summarise_observable_bootstrap(
    point_estimate: DensityObservableSummary,
    samples: NDArray[np.floating],
    *,
    confidence_level: float,
    requested_draws: int,
) -> ObservableBootstrapSummary:
    """Build a missingness-explicit joint summary from aligned bootstrap rows."""

    sample_array = np.asarray(samples, dtype=float)
    if sample_array.ndim != 2 or sample_array.shape[1] != len(
        OBSERVABLE_PARAMETER_NAMES
    ):
        raise ValueError("observable samples have the wrong shape")
    if sample_array.shape[0] == 0:
        raise ValueError("observable bootstrap requires at least one successful draw")
    if requested_draws < sample_array.shape[0]:
        raise ValueError("requested_draws cannot be smaller than successful rows")
    supported = np.isfinite(sample_array)
    point_values = _observable_vector(point_estimate)
    alpha = 0.5 * (1.0 - confidence_level)
    intervals: dict[str, ObservableBootstrapInterval] = {}
    for index, name in enumerate(OBSERVABLE_PARAMETER_NAMES):
        values = sample_array[:, index]
        supported_count = int(np.count_nonzero(supported[:, index]))
        estimate = (
            float(point_values[index]) if np.isfinite(point_values[index]) else None
        )
        if estimate is None or supported_count == 0:
            status: ObservableIntervalStatus = "unsupported"
        elif supported_count < requested_draws:
            status = "partial"
        else:
            status = "complete"
        if status == "complete":
            lower, upper = np.quantile(values, [alpha, 1.0 - alpha])
            bootstrap_mean = float(np.mean(values))
            bootstrap_standard = float(
                np.std(values, ddof=1 if values.size > 1 else 0)
            )
            bootstrap_bias = bootstrap_mean - estimate
            lower_value: float | None = float(lower)
            upper_value: float | None = float(upper)
        else:
            bootstrap_mean = None
            bootstrap_standard = None
            bootstrap_bias = None
            lower_value = None
            upper_value = None
        intervals[name] = ObservableBootstrapInterval(
            estimate=estimate,
            supported_draws=supported_count,
            successful_draws=int(sample_array.shape[0]),
            requested_draws=int(requested_draws),
            status=status,
            bootstrap_mean=bootstrap_mean,
            bootstrap_standard_deviation=bootstrap_standard,
            bootstrap_bias=bootstrap_bias,
            lower=lower_value,
            upper=upper_value,
        )
    joint_supported = int(np.count_nonzero(np.all(supported, axis=1)))
    joint_covariance = (
        np.cov(sample_array, rowvar=False, ddof=1)
        if (
            joint_supported == requested_draws
            and sample_array.shape[0] == requested_draws
            and sample_array.shape[0] > 1
        )
        else None
    )
    return ObservableBootstrapSummary(
        parameter_names=OBSERVABLE_PARAMETER_NAMES,
        requested_draws=int(requested_draws),
        successful_draws=int(sample_array.shape[0]),
        point_estimate=point_estimate,
        samples=sample_array,
        supported_mask=supported,
        intervals=intervals,
        joint_supported_draw_count=joint_supported,
        joint_covariance=joint_covariance,
    )


def parametric_bootstrap_reconstruction(
    measurement: DifferentiableDensityMeasurement,
    model: NonnegativeBilinearDensityModel,
    fit_result: DensityFitResult,
    *,
    coefficient_upper: float | ArrayLike,
    regularisation: CurvatureRegularisation | None,
    draws: int,
    seed: int,
    confidence_level: float = 0.68,
    options: DensityFitOptions | None = None,
    observable_integration_support: ObservableIntegrationSupport | None = None,
    minimum_integrated_response: float = 0.0,
    retain_latent_artifacts: bool = True,
) -> ParametricBootstrapResult:
    """Propagate detector noise through repeated fits around one solution.

    This is a *parametric* bootstrap: synthetic camera counts are drawn from
    the fitted forward prediction.  It quantifies conditional detector-noise
    uncertainty, not optical-model error or the effect of choosing a different
    regulariser.
    """

    if draws <= 0:
        raise ValueError("bootstrap draws must be positive")
    if not (0.0 < confidence_level < 1.0):
        raise ValueError("confidence level must lie between zero and one")
    if not isinstance(retain_latent_artifacts, bool):
        raise TypeError("retain_latent_artifacts must be boolean")
    if observable_integration_support is not None and not isinstance(
        observable_integration_support,
        ObservableIntegrationSupport,
    ):
        raise TypeError(
            "observable_integration_support must be ObservableIntegrationSupport"
        )
    point_observables: DensityObservableSummary | None = None
    if observable_integration_support is not None:
        grid = measurement.grid
        if not np.array_equal(model.y_grid_m, grid.y_grid_m) or not np.array_equal(
            model.z_grid_m,
            grid.z_grid_m,
        ):
            raise ValueError(
                "density model coordinates must match the measurement grid"
            )
        if not np.array_equal(
            observable_integration_support.y_grid_m,
            grid.y_grid_m,
        ) or not np.array_equal(
            observable_integration_support.z_grid_m,
            grid.z_grid_m,
        ):
            raise ValueError(
                "observable integration coordinates must exactly match the "
                "measurement grid"
            )
        if np.any(
            observable_integration_support.support_mask
            & ~np.asarray(model.support_mask, dtype=bool)
        ):
            raise ValueError(
                "observable integration mask cannot extend beyond latent model support"
            )
        point_observables = extract_density_observables(
            fit_result.column_density_m2,
            observable_integration_support,
            minimum_integrated_response=minimum_integrated_response,
        )
    upper = _upper_bound_vector(coefficient_upper, model.parameter_count)
    epsilon = 1e-9 * np.maximum(upper, 1.0)
    initial = np.clip(
        np.asarray(fit_result.coefficients, dtype=float),
        epsilon,
        upper - epsilon,
    )
    expected_channels = measurement.expected_channels_from_density(
        fit_result.column_density_m2
    )
    rng = np.random.default_rng(seed)
    coefficient_samples: list[NDArray[np.floating]] = []
    density_samples: list[NDArray[np.floating]] = []
    feature_samples: list[dict[str, float]] = []
    observable_samples: list[NDArray[np.floating]] = []
    for _ in range(draws):
        noisy = tuple(
            simulate_poisson_gaussian_counts(
                expected,
                read_noise_electrons=measurement.read_noise_electrons,
                rng=rng,
            )
            for expected in expected_channels
        )
        result = fit_nonnegative_basis_density(
            measurement,
            model,
            noisy,
            initial_coefficients=initial,
            coefficient_upper=upper,
            regularisation=regularisation,
            options=options,
        )
        if not result.diagnostics.success:
            continue
        coefficient_samples.append(np.asarray(result.coefficients, dtype=float))
        density_samples.append(np.asarray(result.column_density_m2, dtype=float))
        feature_samples.append(
            _feature_mapping(
                summarise_density_features(result.column_density_m2, measurement)
            )
        )
        if observable_integration_support is not None:
            observable_samples.append(
                _observable_vector(
                    extract_density_observables(
                        result.column_density_m2,
                        observable_integration_support,
                        minimum_integrated_response=minimum_integrated_response,
                    )
                )
            )
    if not density_samples:
        raise RuntimeError("no bootstrap reconstruction converged")
    coefficient_array = np.stack(coefficient_samples)
    density_array = np.stack(density_samples)
    alpha = 0.5 * (1.0 - confidence_level)
    point_features = _feature_mapping(
        summarise_density_features(fit_result.column_density_m2, measurement)
    )
    feature_intervals: dict[str, FeatureInterval] = {}
    for name, estimate in point_features.items():
        values = np.asarray([sample[name] for sample in feature_samples], dtype=float)
        finite = values[np.isfinite(values)]
        if finite.size:
            lower, upper_interval = np.quantile(finite, [alpha, 1.0 - alpha])
            bootstrap_mean = float(np.mean(finite))
            bootstrap_standard = float(
                np.std(finite, ddof=1 if finite.size > 1 else 0)
            )
        else:
            lower = upper_interval = float("nan")
            bootstrap_mean = bootstrap_standard = float("nan")
        feature_intervals[name] = FeatureInterval(
            estimate=float(estimate),
            bootstrap_mean=bootstrap_mean,
            bootstrap_standard_deviation=bootstrap_standard,
            bootstrap_bias=bootstrap_mean - float(estimate),
            lower=float(lower),
            upper=float(upper_interval),
        )
    ddof = 1 if density_array.shape[0] > 1 else 0
    observable_summary: ObservableBootstrapSummary | None = None
    if observable_integration_support is not None:
        assert point_observables is not None
        observable_summary = _summarise_observable_bootstrap(
            point_observables,
            np.stack(observable_samples),
            confidence_level=confidence_level,
            requested_draws=draws,
        )
    return ParametricBootstrapResult(
        requested_draws=int(draws),
        successful_draws=int(density_array.shape[0]),
        seed=int(seed),
        confidence_level=float(confidence_level),
        coefficient_samples=coefficient_array if retain_latent_artifacts else None,
        density_mean_m2=(
            np.mean(density_array, axis=0) if retain_latent_artifacts else None
        ),
        density_standard_uncertainty_m2=(
            np.std(density_array, axis=0, ddof=ddof)
            if retain_latent_artifacts
            else None
        ),
        density_interval_lower_m2=(
            np.quantile(density_array, alpha, axis=0)
            if retain_latent_artifacts
            else None
        ),
        density_interval_upper_m2=(
            np.quantile(density_array, 1.0 - alpha, axis=0)
            if retain_latent_artifacts
            else None
        ),
        feature_intervals=feature_intervals if retain_latent_artifacts else {},
        observables=observable_summary,
        assumptions=(
            "parametric resampling from the fitted camera prediction",
            "declared independent Poisson and Gaussian read noise",
            "fixed forward operator, basis, support and regulariser",
        )
        + (
            ("fixed declared observable integration support for all draws",)
            if observable_integration_support is not None
            else ()
        )
        + ("does not include calibration or instrument-model uncertainty",),
    )


def summarise_reconstruction_stability(
    variants: Mapping[str, DensityFitResult],
    measurement: DifferentiableDensityMeasurement,
) -> ReconstructionStabilitySummary:
    """Summarise spread across deliberately varied inverse assumptions."""

    if len(variants) < 2:
        raise ValueError("at least two reconstruction variants are required")
    labels = tuple(str(label) for label in variants)
    maps = []
    feature_maps: list[dict[str, float]] = []
    for result in variants.values():
        density = np.asarray(result.column_density_m2, dtype=float)
        if density.shape != measurement.grid.y_grid_m.shape:
            raise ValueError("variant density does not match the measurement grid")
        maps.append(density)
        feature_maps.append(
            _feature_mapping(summarise_density_features(density, measurement))
        )
    stack = np.stack(maps)
    ranges: dict[str, tuple[float, float]] = {}
    for name in feature_maps[0]:
        values = np.asarray([features[name] for features in feature_maps], dtype=float)
        finite = values[np.isfinite(values)]
        ranges[name] = (
            (float(np.min(finite)), float(np.max(finite)))
            if finite.size
            else (float("nan"), float("nan"))
        )
    return ReconstructionStabilitySummary(
        variant_labels=labels,
        density_mean_m2=np.mean(stack, axis=0),
        density_standard_deviation_m2=np.std(stack, axis=0, ddof=1),
        density_minimum_m2=np.min(stack, axis=0),
        density_maximum_m2=np.max(stack, axis=0),
        feature_ranges=ranges,
    )
