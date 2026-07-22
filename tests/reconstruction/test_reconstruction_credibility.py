from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

import non_destructive_image.reconstruction.credibility as credibility_module
from non_destructive_image.camera import bin_to_camera_pixels
from non_destructive_image.reconstruction.contracts import (
    DetectorContract,
    FaradayResponseContract,
    ReconstructionGrid,
)
from non_destructive_image.reconstruction.credibility import (
    analyse_local_credibility,
    parametric_bootstrap_reconstruction,
    summarise_density_features,
    summarise_reconstruction_stability,
)
from non_destructive_image.reconstruction.density_fit import (
    DensityFitOptions,
    DensityFitResult,
    fit_nonnegative_basis_density,
)
from non_destructive_image.reconstruction.measurements import (
    DualPortFaradayMeasurement,
)
from non_destructive_image.reconstruction.object_models import (
    NonnegativeBilinearDensityModel,
)
from non_destructive_image.reconstruction.observables import (
    ObservableIntegrationSupport,
)
from non_destructive_image.reconstruction.regularisation import (
    build_curvature_regularisation,
)


@pytest.fixture(scope="module")
def credibility_case() -> dict[str, object]:
    size = 32
    field_of_view_m = 64e-6
    spacing_m = field_of_view_m / size
    axis_m = (np.arange(size) - size / 2) * spacing_m
    y_grid_m, z_grid_m = np.meshgrid(axis_m, axis_m)
    frequency = np.fft.fftfreq(size, d=spacing_m)
    fy, fz = np.meshgrid(frequency, frequency)
    pupil = ((fy**2 + fz**2) <= (2.0e5) ** 2).astype(float)
    bin_size = 2
    camera_y_um = bin_to_camera_pixels(y_grid_m * 1e6, bin_size)[0, :]
    camera_z_um = bin_to_camera_pixels(z_grid_m * 1e6, bin_size)[:, 0]
    camera_y, camera_z = np.meshgrid(camera_y_um, camera_z_um)
    roi_mask = (np.abs(camera_y) <= 24.0) & (np.abs(camera_z) <= 12.0)
    grid = ReconstructionGrid.from_arrays(
        y_grid_m=y_grid_m,
        z_grid_m=z_grid_m,
        pupil=pupil,
        bin_size=bin_size,
        roi_mask=roi_mask,
    )
    measurement = DualPortFaradayMeasurement(
        grid=grid,
        detector=DetectorContract(
            photoelectrons_per_i0_pixel=1200.0,
            read_noise_electrons_per_pixel_per_readout=3.0,
        ),
        response=FaradayResponseContract(
            phase_per_column_density_rad_m2=4.0e-16,
            kappa_f=1.0,
        ),
    )
    model = NonnegativeBilinearDensityModel.from_grid(
        y_grid_m=grid.y_grid_m,
        z_grid_m=grid.z_grid_m,
        knot_y_um=[-18.0, -9.0, 0.0, 9.0, 18.0],
        knot_z_um=[-6.0, 0.0, 6.0],
        coefficient_scale_m2=5.0e14,
    )
    truth = np.asarray(
        [
            [0.02, 0.04, 0.06, 0.04, 0.02],
            [0.04, 0.72, 0.20, 0.58, 0.03],
            [0.01, 0.05, 0.08, 0.04, 0.01],
        ],
        dtype=float,
    ).ravel()
    observed = measurement.expected_channels_from_density(model.column_density(truth))
    options = DensityFitOptions(irls_iterations=2, max_nfev=100)
    unregularised = fit_nonnegative_basis_density(
        measurement,
        model,
        observed,
        initial_coefficients=np.full(model.parameter_count, 0.15),
        coefficient_upper=1.5,
        regularisation=None,
        options=options,
    )
    regularisation = build_curvature_regularisation(
        model.knot_y_um,
        model.knot_z_um,
        density_scale_m2=5.0e14,
        weight_um2=20.0,
    )
    regularised = fit_nonnegative_basis_density(
        measurement,
        model,
        observed,
        initial_coefficients=np.full(model.parameter_count, 0.15),
        coefficient_upper=1.5,
        regularisation=regularisation,
        options=options,
    )
    return {
        "measurement": measurement,
        "model": model,
        "observed": observed,
        "options": options,
        "regularisation": regularisation,
        "unregularised": unregularised,
        "regularised": regularised,
    }


def test_local_credibility_separates_data_and_prior_support(
    credibility_case: dict[str, object],
) -> None:
    measurement = credibility_case["measurement"]
    model = credibility_case["model"]
    observed = credibility_case["observed"]
    unregularised = credibility_case["unregularised"]
    regularised = credibility_case["regularised"]
    regularisation = credibility_case["regularisation"]
    assert isinstance(measurement, DualPortFaradayMeasurement)
    assert isinstance(model, NonnegativeBilinearDensityModel)

    data_only = analyse_local_credibility(
        measurement,
        model,
        observed,
        unregularised,
        coefficient_upper=1.5,
        regularisation=None,
    )
    stabilised = analyse_local_credibility(
        measurement,
        model,
        observed,
        regularised,
        coefficient_upper=1.5,
        regularisation=regularisation,
    )
    assert data_only.data_singular_values.ndim == 1
    assert data_only.data_singular_values[0] >= data_only.data_singular_values[-1]
    assert np.all(data_only.relative_data_singular_values <= 1.0)
    assert np.all(data_only.relative_data_singular_values >= 0.0)
    assert np.all(data_only.generalised_data_mode_fractions <= 1.0 + 1e-12)
    assert np.all(data_only.generalised_data_mode_fractions >= -1e-12)
    assert data_only.effective_prior_degrees_of_freedom == pytest.approx(0.0, abs=1e-7)
    assert data_only.locally_unconstrained_degrees_of_freedom == pytest.approx(
        0.0,
        abs=1e-7,
    )
    assert data_only.combined_constrained_rank == np.count_nonzero(
        data_only.free_coefficient_mask
    )
    assert stabilised.effective_prior_degrees_of_freedom > 0.0
    assert stabilised.effective_data_degrees_of_freedom < np.count_nonzero(
        stabilised.free_coefficient_mask
    )
    assert stabilised.density_standard_uncertainty_m2.shape == model.y_grid_m.shape
    assert np.all(stabilised.density_standard_uncertainty_m2 >= 0.0)
    assert len(stabilised.residual_channels) == 2
    assert all(
        channel.standardised_residual_map.shape == measurement.grid.camera_shape
        for channel in stabilised.residual_channels
    )


def test_low_order_feature_summary_does_not_require_a_shape_model(
    credibility_case: dict[str, object],
) -> None:
    measurement = credibility_case["measurement"]
    regularised = credibility_case["regularised"]
    assert isinstance(measurement, DualPortFaradayMeasurement)
    summary = summarise_density_features(
        regularised.column_density_m2,
        measurement,
    )
    assert summary.integrated_column_density > 0.0
    assert summary.peak_column_density_m2 > 0.0
    assert abs(summary.centroid_y_um) < 10.0
    assert abs(summary.centroid_z_um) < 5.0
    assert summary.rms_y_um > summary.rms_z_um


def test_parametric_bootstrap_reports_conditional_intervals(
    credibility_case: dict[str, object],
) -> None:
    measurement = credibility_case["measurement"]
    model = credibility_case["model"]
    regularised = credibility_case["regularised"]
    regularisation = credibility_case["regularisation"]
    assert isinstance(measurement, DualPortFaradayMeasurement)
    assert isinstance(model, NonnegativeBilinearDensityModel)
    assert isinstance(regularised, DensityFitResult)
    observable_support = ObservableIntegrationSupport(
        y_grid_m=model.y_grid_m,
        z_grid_m=model.z_grid_m,
        support_mask=model.support_mask,
    )
    bootstrap = parametric_bootstrap_reconstruction(
        measurement,
        model,
        regularised,
        coefficient_upper=1.5,
        regularisation=regularisation,
        draws=3,
        seed=20260721,
        confidence_level=0.68,
        options=DensityFitOptions(irls_iterations=1, max_nfev=60),
        observable_integration_support=observable_support,
    )
    assert bootstrap.requested_draws == 3
    assert bootstrap.successful_draws == 3
    assert bootstrap.coefficient_samples.shape == (3, model.parameter_count)
    assert bootstrap.density_mean_m2.shape == model.y_grid_m.shape
    assert np.all(bootstrap.density_standard_uncertainty_m2 >= 0.0)
    interval = bootstrap.feature_intervals["integrated_column_density"]
    assert interval.lower <= interval.upper
    assert np.isfinite(interval.bootstrap_mean)
    assert interval.bootstrap_standard_deviation >= 0.0
    assert interval.bootstrap_bias == pytest.approx(
        interval.bootstrap_mean - interval.estimate
    )
    assert "fixed forward operator, basis, support and regulariser" in bootstrap.assumptions
    assert bootstrap.observables is not None
    observable_bootstrap = bootstrap.observables
    assert observable_bootstrap.parameter_names == (
        "integrated_response",
        "centroid_y_m",
        "centroid_z_m",
        "major_rms_width_m",
    )
    assert observable_bootstrap.samples.shape == (bootstrap.successful_draws, 4)
    assert np.all(observable_bootstrap.supported_mask)
    assert observable_bootstrap.joint_supported_draw_count == bootstrap.successful_draws
    assert observable_bootstrap.joint_covariance is not None
    assert observable_bootstrap.joint_covariance.shape == (4, 4)
    assert observable_bootstrap.point_estimate.integration_support.is_identical_to(
        observable_support
    )
    alpha = 0.5 * (1.0 - bootstrap.confidence_level)
    for index, name in enumerate(observable_bootstrap.parameter_names):
        observable_interval = observable_bootstrap.intervals[name]
        assert observable_interval.status == "complete"
        assert observable_interval.supported_draws == bootstrap.successful_draws
        expected_lower, expected_upper = np.quantile(
            observable_bootstrap.samples[:, index],
            [alpha, 1.0 - alpha],
        )
        assert observable_interval.lower == pytest.approx(expected_lower)
        assert observable_interval.upper == pytest.approx(expected_upper)


def test_parametric_bootstrap_preserves_rows_with_undefined_moments(
    credibility_case: dict[str, object],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    measurement = credibility_case["measurement"]
    model = credibility_case["model"]
    regularised = credibility_case["regularised"]
    regularisation = credibility_case["regularisation"]
    assert isinstance(measurement, DualPortFaradayMeasurement)
    assert isinstance(model, NonnegativeBilinearDensityModel)
    assert isinstance(regularised, DensityFitResult)
    observable_support = ObservableIntegrationSupport(
        y_grid_m=model.y_grid_m,
        z_grid_m=model.z_grid_m,
        support_mask=model.support_mask,
    )
    zero_density = np.zeros_like(regularised.column_density_m2)
    zero_result = DensityFitResult(
        coefficients=np.zeros_like(regularised.coefficients),
        column_density_m2=zero_density,
        predicted_channels=measurement.expected_channels_from_density(zero_density),
        diagnostics=regularised.diagnostics,
    )
    draw_results = iter((regularised, zero_result, regularised))
    monkeypatch.setattr(
        credibility_module,
        "fit_nonnegative_basis_density",
        lambda *args, **kwargs: next(draw_results),
    )

    bootstrap = credibility_module.parametric_bootstrap_reconstruction(
        measurement,
        model,
        regularised,
        coefficient_upper=1.5,
        regularisation=regularisation,
        draws=3,
        seed=20260722,
        confidence_level=0.68,
        options=DensityFitOptions(irls_iterations=1, max_nfev=10),
        observable_integration_support=observable_support,
        retain_latent_artifacts=False,
    )
    assert bootstrap.coefficient_samples is None
    assert bootstrap.density_mean_m2 is None
    assert bootstrap.density_standard_uncertainty_m2 is None
    assert bootstrap.density_interval_lower_m2 is None
    assert bootstrap.density_interval_upper_m2 is None
    assert bootstrap.feature_intervals == {}
    assert bootstrap.observables is not None
    observable_bootstrap = bootstrap.observables
    assert np.all(np.isfinite(observable_bootstrap.samples[:, 0]))
    assert np.isnan(observable_bootstrap.samples[1, 1:]).all()
    assert np.array_equal(
        observable_bootstrap.supported_mask[1],
        np.asarray([True, False, False, False]),
    )
    assert observable_bootstrap.joint_supported_draw_count == 2
    assert observable_bootstrap.joint_covariance is None
    assert observable_bootstrap.intervals["integrated_response"].status == "complete"
    for name in ("centroid_y_m", "centroid_z_m", "major_rms_width_m"):
        interval = observable_bootstrap.intervals[name]
        assert interval.status == "partial"
        assert interval.supported_draws == 2
        assert interval.successful_draws == 3
        assert interval.lower is None
        assert interval.upper is None


def test_failed_bootstrap_refit_makes_every_formal_interval_partial(
    credibility_case: dict[str, object],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    measurement = credibility_case["measurement"]
    model = credibility_case["model"]
    regularised = credibility_case["regularised"]
    regularisation = credibility_case["regularisation"]
    assert isinstance(measurement, DualPortFaradayMeasurement)
    assert isinstance(model, NonnegativeBilinearDensityModel)
    assert isinstance(regularised, DensityFitResult)
    observable_support = ObservableIntegrationSupport(
        y_grid_m=model.y_grid_m,
        z_grid_m=model.z_grid_m,
        support_mask=model.support_mask,
    )
    failed_result = DensityFitResult(
        coefficients=regularised.coefficients,
        column_density_m2=regularised.column_density_m2,
        predicted_channels=regularised.predicted_channels,
        diagnostics=replace(
            regularised.diagnostics,
            success=False,
            message="deliberate bootstrap test failure",
        ),
    )
    draw_results = iter((regularised, failed_result, regularised))
    monkeypatch.setattr(
        credibility_module,
        "fit_nonnegative_basis_density",
        lambda *args, **kwargs: next(draw_results),
    )

    bootstrap = credibility_module.parametric_bootstrap_reconstruction(
        measurement,
        model,
        regularised,
        coefficient_upper=1.5,
        regularisation=regularisation,
        draws=3,
        seed=20260722,
        confidence_level=0.68,
        options=DensityFitOptions(irls_iterations=1, max_nfev=10),
        observable_integration_support=observable_support,
        retain_latent_artifacts=False,
    )
    assert bootstrap.requested_draws == 3
    assert bootstrap.successful_draws == 2
    assert bootstrap.observables is not None
    observable_bootstrap = bootstrap.observables
    assert observable_bootstrap.requested_draws == 3
    assert observable_bootstrap.successful_draws == 2
    assert observable_bootstrap.samples.shape == (2, 4)
    assert observable_bootstrap.joint_covariance is None
    for interval in observable_bootstrap.intervals.values():
        assert interval.status == "partial"
        assert interval.requested_draws == 3
        assert interval.successful_draws == 2
        assert interval.lower is None
        assert interval.upper is None


def test_stability_summary_exposes_choice_dependent_spread(
    credibility_case: dict[str, object],
) -> None:
    measurement = credibility_case["measurement"]
    assert isinstance(measurement, DualPortFaradayMeasurement)
    stability = summarise_reconstruction_stability(
        {
            "data_only": credibility_case["unregularised"],
            "curvature_20_um2": credibility_case["regularised"],
        },
        measurement,
    )
    assert stability.variant_labels == ("data_only", "curvature_20_um2")
    assert stability.density_standard_deviation_m2.shape == measurement.grid.y_grid_m.shape
    assert np.any(stability.density_standard_deviation_m2 > 0.0)
    lower, upper = stability.feature_ranges["peak_column_density_m2"]
    assert lower <= upper
