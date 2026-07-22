"""Focused guards against interpreting prior-supported structure as measured."""

from __future__ import annotations

import numpy as np
import pytest

from non_destructive_image.camera import bin_to_camera_pixels
from non_destructive_image.reconstruction import (
    DarkFieldFaradayMeasurement,
    DensityFitOptions,
    DetectorContract,
    DualPortFaradayMeasurement,
    FaradayResponseContract,
    NonnegativeBilinearDensityModel,
    ReconstructionGrid,
    SmoothTFBounds,
    build_curvature_regularisation,
    dark_field_sqrt_moment_initialisation,
    fit_nonnegative_basis_density,
    linearised_nonnegative_initialisation,
    supported_band_projection,
)
from non_destructive_image.reconstruction.noise import poisson_gaussian_variance


def _small_grid() -> ReconstructionGrid:
    size = 48
    field_of_view_m = 60e-6
    spacing_m = field_of_view_m / size
    axis_m = (np.arange(size) - size / 2) * spacing_m
    y_grid_m, z_grid_m = np.meshgrid(axis_m, axis_m)
    frequency = np.fft.fftfreq(size, d=spacing_m)
    frequency_y, frequency_z = np.meshgrid(frequency, frequency)
    pupil = (
        np.sqrt(frequency_y**2 + frequency_z**2) <= 0.080 / 401e-9
    ).astype(float)
    bin_size = 3
    camera_y_um = bin_to_camera_pixels(y_grid_m * 1e6, bin_size)[0]
    camera_z_um = bin_to_camera_pixels(z_grid_m * 1e6, bin_size)[:, 0]
    camera_y, camera_z = np.meshgrid(camera_y_um, camera_z_um)
    roi = (np.abs(camera_y) <= 21.0) & (np.abs(camera_z) <= 9.0)
    return ReconstructionGrid.from_arrays(
        y_grid_m=y_grid_m,
        z_grid_m=z_grid_m,
        pupil=pupil,
        bin_size=bin_size,
        roi_mask=roi,
    )


@pytest.fixture(scope="module")
def guard_context() -> dict[str, object]:
    grid = _small_grid()
    detector = DetectorContract(
        photoelectrons_per_i0_pixel=700.0,
        read_noise_electrons_per_pixel_per_readout=3.0,
    )
    response = FaradayResponseContract(
        phase_per_column_density_rad_m2=4.0e-16,
        kappa_f=1.0,
    )
    model = NonnegativeBilinearDensityModel.from_grid(
        y_grid_m=grid.y_grid_m,
        z_grid_m=grid.z_grid_m,
        knot_y_um=[-18.0, -9.0, 0.0, 9.0, 18.0],
        knot_z_um=[-6.0, 0.0, 6.0],
        coefficient_scale_m2=5.0e14,
    )
    regularisation = build_curvature_regularisation(
        model.knot_y_um,
        model.knot_z_um,
        density_scale_m2=5.0e14,
        weight_um2=3.0,
    )
    bounds = SmoothTFBounds.from_mapping(
        {
            "column_density_peak_m2": [1.0e12, 1.5e15],
            "y0_um": [-15.0, 15.0],
            "z0_um": [-8.0, 8.0],
            "radius_y_um": [3.0, 30.0],
            "radius_z_um": [1.0, 12.0],
        }
    )
    return {
        "grid": grid,
        "model": model,
        "regularisation": regularisation,
        "bounds": bounds,
        "dual": DualPortFaradayMeasurement(
            grid=grid,
            detector=detector,
            response=response,
        ),
        "dark": DarkFieldFaradayMeasurement(
            grid=grid,
            detector=detector,
            response=response,
        ),
    }


def _fit_from_blind_initialisation(
    measurement: DualPortFaradayMeasurement | DarkFieldFaradayMeasurement,
    model: NonnegativeBilinearDensityModel,
    regularisation,
    bounds: SmoothTFBounds,
    observed: tuple[np.ndarray, ...],
):
    if isinstance(measurement, DualPortFaradayMeasurement):
        initial = linearised_nonnegative_initialisation(
            measurement,
            model,
            observed,
            coefficient_upper=1.5,
        ).coefficients
    else:
        initial = dark_field_sqrt_moment_initialisation(
            measurement,
            model,
            observed[0],
            smooth_bounds=bounds,
            coefficient_upper=1.5,
        ).coefficients
    return fit_nonnegative_basis_density(
        measurement,
        model,
        observed,
        initial_coefficients=initial,
        coefficient_upper=1.5,
        regularisation=regularisation,
        options=DensityFitOptions(irls_iterations=2, max_nfev=80),
    )


@pytest.mark.parametrize("readout", ["dual", "dark"])
def test_blank_ensemble_has_a_frozen_false_positive_baseline(
    guard_context: dict[str, object],
    readout: str,
) -> None:
    """Freeze the current non-negative blank bias for explicit review.

    These values are a regression baseline, not an acceptance threshold.  In
    particular, the much larger dark-field value documents why a real-data
    detection gate must be calibrated from blank frames before morphology is
    interpreted.
    """

    grid = guard_context["grid"]
    model = guard_context["model"]
    measurement = guard_context[readout]
    assert isinstance(grid, ReconstructionGrid)
    assert isinstance(model, NonnegativeBilinearDensityModel)
    assert isinstance(
        measurement,
        (DualPortFaradayMeasurement, DarkFieldFaradayMeasurement),
    )
    zero_density = np.zeros_like(grid.y_grid_m)
    reference_coefficients = np.asarray(
        [
            [0.0, 0.05, 0.10, 0.05, 0.0],
            [0.05, 0.55, 1.00, 0.55, 0.05],
            [0.0, 0.05, 0.10, 0.05, 0.0],
        ]
    ).ravel()
    reference_mass = float(np.sum(model.column_density(reference_coefficients)))
    false_mass_fractions = []
    peak_rotations = []
    for seed in (101, 202, 303, 404):
        observed = measurement.simulate_channels_from_density(
            zero_density,
            np.random.default_rng(seed),
        )
        fitted = _fit_from_blind_initialisation(
            measurement,
            model,
            guard_context["regularisation"],
            guard_context["bounds"],
            observed,
        )
        assert fitted.diagnostics.success
        false_mass_fractions.append(float(np.sum(fitted.column_density_m2)) / reference_mass)
        peak_rotations.append(
            measurement.response.rotation_per_column_density_rad_m2
            * float(np.max(fitted.column_density_m2))
        )

    expected_ranges = {
        "dual": ((0.06, 0.09), (0.04, 0.055)),
        "dark": ((0.40, 0.55), (0.11, 0.16)),
    }[readout]
    median_mass = float(np.median(false_mass_fractions))
    maximum_rotation = max(peak_rotations)
    assert expected_ranges[0][0] < median_mass < expected_ranges[0][1]
    assert expected_ranges[1][0] < maximum_rotation < expected_ranges[1][1]


def test_central_dual_port_reconstruction_is_stable_to_support_trimming(
    guard_context: dict[str, object],
) -> None:
    """A well-contained feature should not be created by the support boundary."""

    grid = guard_context["grid"]
    full_model = guard_context["model"]
    measurement = guard_context["dual"]
    assert isinstance(grid, ReconstructionGrid)
    assert isinstance(full_model, NonnegativeBilinearDensityModel)
    assert isinstance(measurement, DualPortFaradayMeasurement)
    y_um = grid.y_grid_m * 1e6
    z_um = grid.z_grid_m * 1e6
    trimmed_model = NonnegativeBilinearDensityModel.from_grid(
        y_grid_m=grid.y_grid_m,
        z_grid_m=grid.z_grid_m,
        knot_y_um=full_model.knot_y_um,
        knot_z_um=full_model.knot_z_um,
        coefficient_scale_m2=full_model.coefficient_scale_m2,
        support_mask=(np.abs(y_um) <= 15.0) & (np.abs(z_um) <= 5.0),
    )
    truth = 5.0e14 * np.exp(-0.5 * ((y_um / 7.0) ** 2 + (z_um / 2.4) ** 2))
    observed = measurement.simulate_channels_from_density(
        truth,
        np.random.default_rng(607),
    )

    fits = []
    for model in (full_model, trimmed_model):
        initial = linearised_nonnegative_initialisation(
            measurement,
            model,
            observed,
            coefficient_upper=1.5,
        ).coefficients
        regularisation = build_curvature_regularisation(
            model.knot_y_um,
            model.knot_z_um,
            density_scale_m2=5.0e14,
            weight_um2=3.0,
        )
        fits.append(
            fit_nonnegative_basis_density(
                measurement,
                model,
                observed,
                initial_coefficients=initial,
                coefficient_upper=1.5,
                regularisation=regularisation,
                options=DensityFitOptions(irls_iterations=2, max_nfev=80),
            )
        )
    assert all(fit.diagnostics.success for fit in fits)
    full_supported = supported_band_projection(fits[0].column_density_m2, grid.pupil)
    trimmed_supported = supported_band_projection(fits[1].column_density_m2, grid.pupil)
    relative_change = np.linalg.norm(trimmed_supported - full_supported) / np.linalg.norm(
        full_supported
    )
    assert relative_change < 0.12


def test_dark_field_solution_is_stable_across_distinct_feasible_starts(
    guard_context: dict[str, object],
) -> None:
    """A single dark-field envelope must not select a different local solution."""

    model = guard_context["model"]
    measurement = guard_context["dark"]
    assert isinstance(model, NonnegativeBilinearDensityModel)
    assert isinstance(measurement, DarkFieldFaradayMeasurement)
    truth_coefficients = np.asarray(
        [
            [0.01, 0.05, 0.09, 0.05, 0.01],
            [0.03, 0.75, 0.18, 0.62, 0.03],
            [0.01, 0.04, 0.08, 0.04, 0.01],
        ]
    ).ravel()
    observed = measurement.expected_channels_from_density(
        model.column_density(truth_coefficients)
    )
    blind = dark_field_sqrt_moment_initialisation(
        measurement,
        model,
        observed[0],
        smooth_bounds=guard_context["bounds"],
        coefficient_upper=1.5,
    ).coefficients
    starts = (
        blind,
        np.full(model.parameter_count, 0.18),
        np.asarray(
            [
                [0.02, 0.20, 0.05, 0.02, 0.01],
                [0.05, 0.15, 0.55, 0.25, 0.04],
                [0.01, 0.03, 0.08, 0.22, 0.02],
            ]
        ).ravel(),
    )
    fits = tuple(
        fit_nonnegative_basis_density(
            measurement,
            model,
            observed,
            initial_coefficients=start,
            coefficient_upper=1.5,
            regularisation=guard_context["regularisation"],
            options=DensityFitOptions(irls_iterations=2, max_nfev=100),
        )
        for start in starts
    )
    assert all(fit.diagnostics.success for fit in fits)
    reference_prediction = measurement.flatten_observed(*fits[0].predicted_channels)
    reference_density = supported_band_projection(fits[0].column_density_m2, measurement.grid.pupil)
    for fit in fits[1:]:
        prediction = measurement.flatten_observed(*fit.predicted_channels)
        prediction_change = np.linalg.norm(prediction - reference_prediction) / max(
            np.linalg.norm(reference_prediction),
            np.finfo(float).eps,
        )
        density = supported_band_projection(fit.column_density_m2, measurement.grid.pupil)
        density_change = np.linalg.norm(density - reference_density) / np.linalg.norm(
            reference_density
        )
        assert prediction_change < 1e-4
        assert density_change < 0.08


def test_out_of_band_density_pair_is_indistinguishable_at_the_camera(
    guard_context: dict[str, object],
) -> None:
    """Distinct null-band truths must not be treated as separately recovered."""

    grid = guard_context["grid"]
    assert isinstance(grid, ReconstructionGrid)
    index = np.arange(grid.y_grid_m.shape[1])
    out_of_band_mode = np.cos(2.0 * np.pi * 18.0 * index / index.size)[None, :]
    baseline = 4.0e14
    modulation = 0.30
    truth_a = baseline * (1.0 + modulation * out_of_band_mode)
    truth_b = baseline * (1.0 - modulation * out_of_band_mode)
    truth_a = np.broadcast_to(truth_a, grid.y_grid_m.shape).copy()
    truth_b = np.broadcast_to(truth_b, grid.y_grid_m.shape).copy()

    full_difference = np.linalg.norm(truth_a - truth_b) / np.linalg.norm(truth_a)
    supported_a = supported_band_projection(truth_a, grid.pupil)
    supported_b = supported_band_projection(truth_b, grid.pupil)
    supported_difference = np.linalg.norm(supported_a - supported_b) / np.linalg.norm(
        supported_a
    )
    assert full_difference > 0.35
    assert supported_difference < 1e-12

    for readout in ("dual", "dark"):
        measurement = guard_context[readout]
        assert isinstance(
            measurement,
            (DualPortFaradayMeasurement, DarkFieldFaradayMeasurement),
        )
        expected_a = measurement.expected_vector_from_density(truth_a)
        expected_b = measurement.expected_vector_from_density(truth_b)
        variance = poisson_gaussian_variance(
            0.5 * (expected_a + expected_b),
            measurement.read_noise_electrons,
        )
        squared_separation = float(np.sum((expected_a - expected_b) ** 2 / variance))
        assert squared_separation < 1.0
