"""Synthetic model-mismatch benchmarks for density reconstruction."""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import fft as scipy_fft

from .contracts import ReconstructionGrid
from .density_fit import DensityFitOptions, DensityFitResult, fit_nonnegative_basis_density
from .density_initialise import (
    DensityInitialisation,
    dark_field_sqrt_moment_initialisation,
    linearised_nonnegative_initialisation,
)
from .measurements import DarkFieldFaradayMeasurement, DifferentiableDensityMeasurement
from .object_models import NonnegativeBilinearDensityModel
from .parameters import SmoothTFBounds
from .regularisation import CurvatureRegularisation
from .synthetic_morphologies import SyntheticMorphology


@dataclass(frozen=True)
class DensityRecoveryMetrics:
    """Map and low-order observable errors for known synthetic truth."""

    full_map_relative_l2_error: float
    supported_band_relative_l2_error: float
    integrated_density_relative_error: float
    centroid_y_error_um: float
    centroid_z_error_um: float
    rms_y_relative_error: float
    rms_z_relative_error: float


@dataclass(frozen=True)
class MorphologyBenchmarkResult:
    """One blind deterministic reconstruction of an independent truth map."""

    morphology: SyntheticMorphology
    initialisation: DensityInitialisation
    fit: DensityFitResult
    metrics: DensityRecoveryMetrics


def supported_band_projection(
    column_density_m2: ArrayLike,
    pupil: ArrayLike,
    *,
    relative_threshold: float = 1e-12,
) -> NDArray[np.floating]:
    """Project a density map onto spatial frequencies transmitted by the pupil.

    Only the support of the pupil is used.  Aberration phase is part of the
    forward measurement, not part of this comparison target.
    """

    density = np.asarray(column_density_m2, dtype=float)
    pupil_array = np.asarray(pupil)
    if density.ndim != 2 or density.shape != pupil_array.shape:
        raise ValueError("density and pupil must be same-shape 2D arrays")
    maximum = float(np.max(np.abs(pupil_array)))
    if maximum <= 0:
        raise ValueError("pupil transmits no spatial frequencies")
    support = np.abs(pupil_array) > relative_threshold * maximum
    return np.real(scipy_fft.ifft2(scipy_fft.fft2(density) * support))


def _density_moments(
    density: NDArray[np.floating],
    grid: ReconstructionGrid,
) -> tuple[float, float, float, float, float]:
    total = float(np.sum(density))
    if total <= 0:
        raise ValueError("density moments require a positive integrated density")
    y_um = grid.y_grid_m * 1e6
    z_um = grid.z_grid_m * 1e6
    centre_y = float(np.sum(density * y_um) / total)
    centre_z = float(np.sum(density * z_um) / total)
    rms_y = float(np.sqrt(np.sum(density * (y_um - centre_y) ** 2) / total))
    rms_z = float(np.sqrt(np.sum(density * (z_um - centre_z) ** 2) / total))
    return total, centre_y, centre_z, rms_y, rms_z


def density_recovery_metrics(
    truth_column_density_m2: ArrayLike,
    recovered_column_density_m2: ArrayLike,
    grid: ReconstructionGrid,
) -> DensityRecoveryMetrics:
    """Compare full maps and the density content supported by the aperture."""

    truth = np.asarray(truth_column_density_m2, dtype=float)
    recovered = np.asarray(recovered_column_density_m2, dtype=float)
    if truth.shape != grid.y_grid_m.shape or recovered.shape != truth.shape:
        raise ValueError("truth and recovered densities must match the reconstruction grid")
    if np.any(~np.isfinite(truth)) or np.any(~np.isfinite(recovered)):
        raise ValueError("density maps must be finite")
    if np.any(truth < 0) or np.any(recovered < 0):
        raise ValueError("density maps must be non-negative")

    full_denominator = max(float(np.linalg.norm(truth)), np.finfo(float).eps)
    truth_supported = supported_band_projection(truth, grid.pupil)
    recovered_supported = supported_band_projection(recovered, grid.pupil)
    supported_denominator = max(
        float(np.linalg.norm(truth_supported)),
        np.finfo(float).eps,
    )
    truth_moments = _density_moments(truth, grid)
    recovered_moments = _density_moments(recovered, grid)
    return DensityRecoveryMetrics(
        full_map_relative_l2_error=float(np.linalg.norm(recovered - truth) / full_denominator),
        supported_band_relative_l2_error=float(
            np.linalg.norm(recovered_supported - truth_supported) / supported_denominator
        ),
        integrated_density_relative_error=recovered_moments[0] / truth_moments[0] - 1.0,
        centroid_y_error_um=recovered_moments[1] - truth_moments[1],
        centroid_z_error_um=recovered_moments[2] - truth_moments[2],
        rms_y_relative_error=recovered_moments[3] / truth_moments[3] - 1.0,
        rms_z_relative_error=recovered_moments[4] / truth_moments[4] - 1.0,
    )


def run_linear_readout_morphology_benchmark(
    measurement: DifferentiableDensityMeasurement,
    model: NonnegativeBilinearDensityModel,
    morphologies: Iterable[SyntheticMorphology],
    *,
    coefficient_upper: float | ArrayLike,
    regularisation: CurvatureRegularisation | None,
    initialisation_ridge_strength: float = 0.0,
    fit_options: DensityFitOptions | None = None,
) -> tuple[MorphologyBenchmarkResult, ...]:
    """Blindly reconstruct deterministic independent truths through one readout.

    The initial estimate is obtained from the raw counts and the instrument
    Jacobian at zero density.  This benchmark therefore applies to readouts
    with first-order density information.  Quadratic readouts require a
    separate nonlinear initialisation and are rejected by the initializer.
    """

    return run_morphology_benchmark(
        measurement,
        model,
        morphologies,
        initialise=lambda observed: linearised_nonnegative_initialisation(
            measurement,
            model,
            observed,
            coefficient_upper=coefficient_upper,
            ridge_strength=initialisation_ridge_strength,
        ),
        coefficient_upper=coefficient_upper,
        regularisation=regularisation,
        fit_options=fit_options,
    )


def run_dark_field_morphology_benchmark(
    measurement: DarkFieldFaradayMeasurement,
    model: NonnegativeBilinearDensityModel,
    morphologies: Iterable[SyntheticMorphology],
    *,
    smooth_bounds: SmoothTFBounds,
    coefficient_upper: float | ArrayLike,
    regularisation: CurvatureRegularisation | None,
    projection_ridge_strength: float = 0.0,
    fit_options: DensityFitOptions | None = None,
) -> tuple[MorphologyBenchmarkResult, ...]:
    """Blind deterministic benchmark using the dark-field square-root seed."""

    return run_morphology_benchmark(
        measurement,
        model,
        morphologies,
        initialise=lambda observed: dark_field_sqrt_moment_initialisation(
            measurement,
            model,
            observed[0],
            smooth_bounds=smooth_bounds,
            coefficient_upper=coefficient_upper,
            projection_ridge_strength=projection_ridge_strength,
        ),
        coefficient_upper=coefficient_upper,
        regularisation=regularisation,
        fit_options=fit_options,
    )


def run_morphology_benchmark(
    measurement: DifferentiableDensityMeasurement,
    model: NonnegativeBilinearDensityModel,
    morphologies: Iterable[SyntheticMorphology],
    *,
    initialise: Callable[[tuple[NDArray[np.floating], ...]], DensityInitialisation],
    coefficient_upper: float | ArrayLike,
    regularisation: CurvatureRegularisation | None,
    fit_options: DensityFitOptions | None = None,
) -> tuple[MorphologyBenchmarkResult, ...]:
    """Run a deterministic model-mismatch benchmark with an explicit initializer."""

    results: list[MorphologyBenchmarkResult] = []
    for morphology in morphologies:
        truth = morphology.column_density_m2
        if truth.shape != measurement.grid.y_grid_m.shape:
            raise ValueError(
                f"morphology {morphology.name!r} does not match the measurement grid"
            )
        observed = measurement.expected_channels_from_density(truth)
        initialisation = initialise(observed)
        fit = fit_nonnegative_basis_density(
            measurement,
            model,
            observed,
            initial_coefficients=initialisation.coefficients,
            coefficient_upper=coefficient_upper,
            regularisation=regularisation,
            options=fit_options,
        )
        results.append(
            MorphologyBenchmarkResult(
                morphology=morphology,
                initialisation=initialisation,
                fit=fit,
                metrics=density_recovery_metrics(truth, fit.column_density_m2, measurement.grid),
            )
        )
    return tuple(results)
