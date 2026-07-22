from __future__ import annotations

import numpy as np
import pytest

from non_destructive_image.reconstruction import (
    DarkFieldFaradayMeasurement,
    DensityFitOptions,
    DetectorContract,
    DualPortFaradayMeasurement,
    FaradayResponseContract,
    NonnegativeBilinearDensityModel,
    ReconstructionGrid,
    ReconstructionCandidate,
    SmoothTFBounds,
    build_reference_morphology_suite,
    density_recovery_metrics,
    dark_field_sqrt_moment_initialisation,
    fit_nonnegative_basis_density,
    assess_reconstruction_candidate,
    build_camera_aligned_reduced_grid,
    build_calibration_held_out_morphology_suite,
    build_curvature_regularisation,
    build_uniform_reconstruction_grid,
    evaluate_frozen_candidate_on_held_out,
    faraday_grid_agreement,
    generate_noisy_observation_ensemble,
    linearised_nonnegative_initialisation,
    make_linear_candidate_initialiser,
    run_dark_field_morphology_benchmark,
    run_linear_readout_morphology_benchmark,
    select_and_freeze_candidate,
)
from non_destructive_image.camera import bin_to_camera_pixels


def test_camera_aligned_reduced_grid_preserves_canonical_faraday_readouts() -> None:
    canonical = build_uniform_reconstruction_grid(
        ngrid=1024,
        field_of_view_m=100e-6,
        bin_size=15,
        numerical_aperture=0.080,
        wavelength_m=401e-9,
    )
    reduced, provenance = build_camera_aligned_reduced_grid(
        canonical_ngrid=1024,
        canonical_field_of_view_m=100e-6,
        canonical_bin_size=15,
        reduced_bin_size=3,
        numerical_aperture=0.080,
        wavelength_m=401e-9,
    )
    assert canonical.camera_shape == reduced.camera_shape == (68, 68)
    assert provenance.reduced_ngrid == 204
    assert provenance.maximum_camera_coordinate_mismatch_m < 1e-15

    detector = DetectorContract(
        photoelectrons_per_i0_pixel=1.0,
        read_noise_electrons_per_pixel_per_readout=0.0,
    )
    response = FaradayResponseContract(
        phase_per_column_density_rad_m2=3.77498270624925e-16,
        kappa_f=1.0,
    )
    canonical_truth = build_reference_morphology_suite(
        canonical.y_grid_m,
        canonical.z_grid_m,
        peak_column_density_m2=5.3759624525784675e14,
        radius_y_um=25.0,
        radius_z_um=6.0,
    )[3]
    reduced_truth = build_reference_morphology_suite(
        reduced.y_grid_m,
        reduced.z_grid_m,
        peak_column_density_m2=5.3759624525784675e14,
        radius_y_um=25.0,
        radius_z_um=6.0,
    )[3]
    canonical_dual = DualPortFaradayMeasurement(
        grid=canonical,
        detector=detector,
        response=response,
    )
    reduced_dual = DualPortFaradayMeasurement(
        grid=reduced,
        detector=detector,
        response=response,
    )
    canonical_dark = DarkFieldFaradayMeasurement(
        grid=canonical,
        detector=detector,
        response=response,
    )
    reduced_dark = DarkFieldFaradayMeasurement(
        grid=reduced,
        detector=detector,
        response=response,
    )
    canonical_h, canonical_v = canonical_dual.expected_channels_from_density(
        canonical_truth.column_density_m2
    )
    reduced_h, reduced_v = reduced_dual.expected_channels_from_density(
        reduced_truth.column_density_m2
    )
    agreement = faraday_grid_agreement(
        canonical_h=canonical_h,
        canonical_v=canonical_v,
        reduced_h=reduced_h,
        reduced_v=reduced_v,
        canonical_dark=canonical_dark.expected_channels_from_density(
            canonical_truth.column_density_m2
        )[0],
        reduced_dark=reduced_dark.expected_channels_from_density(
            reduced_truth.column_density_m2
        )[0],
    )
    assert agreement.dual_port_signal_relative_l2_error < 0.007
    assert agreement.dual_port_atom_dependent_channels_relative_l2_error < 0.007
    assert agreement.dark_field_relative_l2_error < 0.007
    assert abs(agreement.dual_port_peak_relative_error) < 0.004
    assert abs(agreement.dark_field_peak_relative_error) < 0.007


@pytest.fixture(scope="module")
def reconstruction_context() -> dict[str, object]:
    size = 64
    field_of_view_m = 80e-6
    spacing_m = field_of_view_m / size
    axis_m = (np.arange(size) - size / 2) * spacing_m
    y_grid_m, z_grid_m = np.meshgrid(axis_m, axis_m)
    frequency = np.fft.fftfreq(size, d=spacing_m)
    fy, fz = np.meshgrid(frequency, frequency)
    pupil = ((fy**2 + fz**2) <= (2.2e5) ** 2).astype(float)
    bin_size = 4
    camera_y_um = bin_to_camera_pixels(y_grid_m * 1e6, bin_size)[0, :]
    camera_z_um = bin_to_camera_pixels(z_grid_m * 1e6, bin_size)[:, 0]
    camera_y, camera_z = np.meshgrid(camera_y_um, camera_z_um)
    roi_mask = (np.abs(camera_y) <= 27.0) & (np.abs(camera_z) <= 16.0)
    grid = ReconstructionGrid.from_arrays(
        y_grid_m=y_grid_m,
        z_grid_m=z_grid_m,
        pupil=pupil,
        bin_size=bin_size,
        roi_mask=roi_mask,
    )
    detector = DetectorContract(
        photoelectrons_per_i0_pixel=900.0,
        read_noise_electrons_per_pixel_per_readout=3.0,
    )
    response = FaradayResponseContract(
        phase_per_column_density_rad_m2=4.0e-16,
        kappa_f=1.0,
    )
    return {
        "grid": grid,
        "dual": DualPortFaradayMeasurement(grid=grid, detector=detector, response=response),
        "dark": DarkFieldFaradayMeasurement(grid=grid, detector=detector, response=response),
    }


@pytest.fixture(scope="module")
def bilinear_model(reconstruction_context: dict[str, object]) -> NonnegativeBilinearDensityModel:
    grid = reconstruction_context["grid"]
    assert isinstance(grid, ReconstructionGrid)
    return NonnegativeBilinearDensityModel.from_grid(
        y_grid_m=grid.y_grid_m,
        z_grid_m=grid.z_grid_m,
        knot_y_um=[-18.0, -9.0, 0.0, 9.0, 18.0],
        knot_z_um=[-6.0, 0.0, 6.0],
        coefficient_scale_m2=5.0e14,
    )


@pytest.fixture(scope="module")
def double_lobe_coefficients() -> np.ndarray:
    return np.asarray(
        [
            [0.02, 0.08, 0.12, 0.05, 0.01],
            [0.05, 0.90, 0.18, 0.72, 0.04],
            [0.01, 0.06, 0.14, 0.09, 0.02],
        ],
        dtype=float,
    ).ravel()


def test_bilinear_density_is_nonnegative_and_finite_support(
    bilinear_model: NonnegativeBilinearDensityModel,
    double_lobe_coefficients: np.ndarray,
) -> None:
    density = bilinear_model.column_density(double_lobe_coefficients)
    assert np.all(density >= 0)
    outside = (
        (bilinear_model.y_grid_m * 1e6 < bilinear_model.knot_y_um[0])
        | (bilinear_model.y_grid_m * 1e6 > bilinear_model.knot_y_um[-1])
        | (bilinear_model.z_grid_m * 1e6 < bilinear_model.knot_z_um[0])
        | (bilinear_model.z_grid_m * 1e6 > bilinear_model.knot_z_um[-1])
    )
    assert np.count_nonzero(density[outside]) == 0
    central_row = density[density.shape[0] // 2]
    y_axis_um = bilinear_model.y_grid_m[density.shape[0] // 2] * 1e6
    left = int(np.argmin(np.abs(y_axis_um + 9.0)))
    centre = int(np.argmin(np.abs(y_axis_um)))
    right = int(np.argmin(np.abs(y_axis_um - 9.0)))
    assert central_row[left] > 3.0 * central_row[centre]
    assert central_row[right] > 2.0 * central_row[centre]


@pytest.mark.parametrize("measurement_key", ["dual", "dark"])
def test_bilinear_object_jacobian_matches_finite_difference(
    reconstruction_context: dict[str, object],
    bilinear_model: NonnegativeBilinearDensityModel,
    double_lobe_coefficients: np.ndarray,
    measurement_key: str,
) -> None:
    measurement = reconstruction_context[measurement_key]
    assert isinstance(
        measurement,
        (DualPortFaradayMeasurement, DarkFieldFaradayMeasurement),
    )
    _, analytic = measurement.expected_vector_and_jacobian_model(
        bilinear_model,
        double_lobe_coefficients,
    )
    for index in (1, 7, 13):
        step = 2e-6
        plus = double_lobe_coefficients.copy()
        minus = double_lobe_coefficients.copy()
        plus[index] += step
        minus[index] -= step
        numerical = (
            measurement.expected_vector_from_density(bilinear_model.column_density(plus))
            - measurement.expected_vector_from_density(bilinear_model.column_density(minus))
        ) / (2.0 * step)
        relative_error = np.linalg.norm(numerical - analytic[:, index]) / max(
            np.linalg.norm(numerical),
            np.finfo(float).eps,
        )
        assert relative_error < 4e-6


def test_dual_port_fit_recovers_resolved_non_tf_double_lobe(
    reconstruction_context: dict[str, object],
    bilinear_model: NonnegativeBilinearDensityModel,
    double_lobe_coefficients: np.ndarray,
) -> None:
    measurement = reconstruction_context["dual"]
    assert isinstance(measurement, DualPortFaradayMeasurement)
    truth_density = bilinear_model.column_density(double_lobe_coefficients)
    result = fit_nonnegative_basis_density(
        measurement,
        bilinear_model,
        measurement.expected_channels_from_density(truth_density),
        initial_coefficients=np.full(bilinear_model.parameter_count, 0.20),
        coefficient_upper=1.5,
        regularisation=None,
        options=DensityFitOptions(irls_iterations=2, max_nfev=120),
    )
    assert result.diagnostics.success
    assert result.diagnostics.data_jacobian_rank == bilinear_model.parameter_count
    relative_error = np.linalg.norm(result.column_density_m2 - truth_density) / np.linalg.norm(
        truth_density
    )
    assert relative_error < 2e-5


def test_dark_field_fit_recovers_resolved_non_tf_double_lobe_from_nearby_start(
    reconstruction_context: dict[str, object],
    bilinear_model: NonnegativeBilinearDensityModel,
    double_lobe_coefficients: np.ndarray,
) -> None:
    measurement = reconstruction_context["dark"]
    assert isinstance(measurement, DarkFieldFaradayMeasurement)
    truth_density = bilinear_model.column_density(double_lobe_coefficients)
    result = fit_nonnegative_basis_density(
        measurement,
        bilinear_model,
        measurement.expected_channels_from_density(truth_density),
        initial_coefficients=0.92 * double_lobe_coefficients + 0.005,
        coefficient_upper=1.5,
        regularisation=None,
        options=DensityFitOptions(irls_iterations=2, max_nfev=140),
    )
    assert result.diagnostics.success
    assert result.diagnostics.data_jacobian_rank == bilinear_model.parameter_count
    relative_error = np.linalg.norm(result.column_density_m2 - truth_density) / np.linalg.norm(
        truth_density
    )
    assert relative_error < 3e-5


def test_curvature_regularisation_is_integrated_into_density_fit(
    reconstruction_context: dict[str, object],
    bilinear_model: NonnegativeBilinearDensityModel,
    double_lobe_coefficients: np.ndarray,
) -> None:
    measurement = reconstruction_context["dual"]
    assert isinstance(measurement, DualPortFaradayMeasurement)
    truth_density = bilinear_model.column_density(double_lobe_coefficients)
    observed = measurement.expected_channels_from_density(truth_density)
    initial = np.full(bilinear_model.parameter_count, 0.20)
    unregularised = fit_nonnegative_basis_density(
        measurement,
        bilinear_model,
        observed,
        initial_coefficients=initial,
        coefficient_upper=1.5,
        regularisation=None,
        options=DensityFitOptions(irls_iterations=2, max_nfev=120),
    )
    regularisation = build_curvature_regularisation(
        bilinear_model.knot_y_um,
        bilinear_model.knot_z_um,
        density_scale_m2=5.0e14,
        weight_um2=30.0,
    )
    regularised = fit_nonnegative_basis_density(
        measurement,
        bilinear_model,
        observed,
        initial_coefficients=initial,
        coefficient_upper=1.5,
        regularisation=regularisation,
        options=DensityFitOptions(irls_iterations=2, max_nfev=120),
    )
    comparison_operator = build_curvature_regularisation(
        bilinear_model.knot_y_um,
        bilinear_model.knot_z_um,
        density_scale_m2=5.0e14,
        weight_um2=1.0,
    )
    unregularised_curvature = np.linalg.norm(
        comparison_operator.residual_from_coefficients(
            unregularised.coefficients,
            coefficient_scale_m2=bilinear_model.coefficient_scale_m2,
        )
    )
    regularised_curvature = np.linalg.norm(
        comparison_operator.residual_from_coefficients(
            regularised.coefficients,
            coefficient_scale_m2=bilinear_model.coefficient_scale_m2,
        )
    )
    assert regularised.diagnostics.curvature_weight_um2 == 30.0
    assert regularised.diagnostics.regularisation_boundary_policy == (
        "fixed_zero_ghost_knots"
    )
    assert regularised_curvature < unregularised_curvature
    assert regularised.diagnostics.weighted_chi_square >= (
        unregularised.diagnostics.weighted_chi_square
    )


def test_aperture_null_is_not_reported_as_recovered_information(
    reconstruction_context: dict[str, object],
    bilinear_model: NonnegativeBilinearDensityModel,
    double_lobe_coefficients: np.ndarray,
) -> None:
    measurement = reconstruction_context["dual"]
    assert isinstance(measurement, DualPortFaradayMeasurement)
    null_grid = ReconstructionGrid.from_arrays(
        y_grid_m=measurement.grid.y_grid_m,
        z_grid_m=measurement.grid.z_grid_m,
        pupil=np.zeros_like(measurement.grid.pupil),
        bin_size=measurement.grid.bin_size,
        roi_mask=measurement.grid.roi_mask,
    )
    null_measurement = DualPortFaradayMeasurement(
        grid=null_grid,
        detector=measurement.detector,
        response=measurement.response,
    )
    _, jacobian = null_measurement.expected_vector_and_jacobian_model(
        bilinear_model,
        double_lobe_coefficients,
    )
    assert np.linalg.matrix_rank(jacobian) == 0
    assert np.count_nonzero(jacobian) == 0


def test_reference_morphologies_are_independent_nonnegative_truth_maps(
    reconstruction_context: dict[str, object],
) -> None:
    grid = reconstruction_context["grid"]
    assert isinstance(grid, ReconstructionGrid)
    suite = build_reference_morphology_suite(
        grid.y_grid_m,
        grid.z_grid_m,
        peak_column_density_m2=5.0e14,
        radius_y_um=18.0,
        radius_z_um=5.5,
    )
    assert {case.feature_class for case in suite} == {
        "smooth",
        "asymmetry",
        "modulation",
        "fragmentation",
        "local_defect",
    }
    assert all(case.column_density_m2.shape == grid.y_grid_m.shape for case in suite)
    assert all(np.all(case.column_density_m2 >= 0) for case in suite)
    assert all(np.max(case.column_density_m2) == pytest.approx(5.0e14) for case in suite)


def test_calibration_and_held_out_suites_cover_same_classes_without_map_reuse(
    reconstruction_context: dict[str, object],
) -> None:
    grid = reconstruction_context["grid"]
    assert isinstance(grid, ReconstructionGrid)
    split = build_calibration_held_out_morphology_suite(
        grid.y_grid_m,
        grid.z_grid_m,
        peak_column_density_m2=5.0e14,
        radius_y_um=18.0,
        radius_z_um=5.5,
    )
    calibration_classes = {case.feature_class for case in split.calibration}
    held_out_classes = {case.feature_class for case in split.held_out}
    assert calibration_classes == held_out_classes == {
        "smooth",
        "asymmetry",
        "modulation",
        "fragmentation",
        "local_defect",
    }
    assert {case.name for case in split.calibration}.isdisjoint(
        case.name for case in split.held_out
    )
    assert all(np.max(case.column_density_m2) == pytest.approx(5.0e14) for case in split.calibration)
    assert all(np.max(case.column_density_m2) == pytest.approx(5.0e14) for case in split.held_out)


def test_identical_maps_have_zero_full_and_supported_band_error(
    reconstruction_context: dict[str, object],
) -> None:
    grid = reconstruction_context["grid"]
    assert isinstance(grid, ReconstructionGrid)
    truth = build_reference_morphology_suite(
        grid.y_grid_m,
        grid.z_grid_m,
        peak_column_density_m2=5.0e14,
    )[0].column_density_m2
    metrics = density_recovery_metrics(truth, truth, grid)
    assert metrics.full_map_relative_l2_error == pytest.approx(0.0, abs=1e-15)
    assert metrics.supported_band_relative_l2_error == pytest.approx(0.0, abs=1e-15)
    assert metrics.integrated_density_relative_error == pytest.approx(0.0, abs=1e-15)


def test_linearised_initialiser_rejects_quadratic_dark_field_origin(
    reconstruction_context: dict[str, object],
    bilinear_model: NonnegativeBilinearDensityModel,
    double_lobe_coefficients: np.ndarray,
) -> None:
    measurement = reconstruction_context["dark"]
    assert isinstance(measurement, DarkFieldFaradayMeasurement)
    observed = measurement.expected_channels_from_density(
        bilinear_model.column_density(double_lobe_coefficients)
    )
    with pytest.raises(ValueError, match="no first-order density information"):
        linearised_nonnegative_initialisation(
            measurement,
            bilinear_model,
            observed,
            coefficient_upper=1.5,
        )


def test_dark_field_blind_seed_supports_smooth_model_mismatch_fit(
    reconstruction_context: dict[str, object],
    bilinear_model: NonnegativeBilinearDensityModel,
) -> None:
    measurement = reconstruction_context["dark"]
    grid = reconstruction_context["grid"]
    assert isinstance(measurement, DarkFieldFaradayMeasurement)
    assert isinstance(grid, ReconstructionGrid)
    morphology = build_reference_morphology_suite(
        grid.y_grid_m,
        grid.z_grid_m,
        peak_column_density_m2=5.0e14,
        radius_y_um=18.0,
        radius_z_um=5.5,
    )[0]
    bounds = SmoothTFBounds.from_mapping(
        {
            "column_density_peak_m2": [1.0e13, 5.0e15],
            "y0_um": [-10.0, 10.0],
            "z0_um": [-8.0, 8.0],
            "radius_y_um": [5.0, 50.0],
            "radius_z_um": [0.5, 15.0],
        }
    )
    direct_initialisation = dark_field_sqrt_moment_initialisation(
        measurement,
        bilinear_model,
        measurement.expected_channels_from_density(morphology.column_density_m2)[0],
        smooth_bounds=bounds,
        coefficient_upper=1.5,
        projection_ridge_strength=1e-8,
    )
    results = run_dark_field_morphology_benchmark(
        measurement,
        bilinear_model,
        [morphology],
        smooth_bounds=bounds,
        coefficient_upper=1.5,
        regularisation=None,
        projection_ridge_strength=1e-8,
        fit_options=DensityFitOptions(irls_iterations=2, max_nfev=180),
    )
    benchmark = results[0]
    np.testing.assert_allclose(
        benchmark.initialisation.coefficients,
        direct_initialisation.coefficients,
    )
    assert benchmark.initialisation.method == "dark_field_sqrt_moment_envelope"
    assert benchmark.initialisation.data_jacobian_rank == bilinear_model.parameter_count
    assert benchmark.fit.diagnostics.success
    assert benchmark.metrics.supported_band_relative_l2_error < 0.30


def test_blind_dual_port_benchmark_exposes_morphology_mismatch(
    reconstruction_context: dict[str, object],
    bilinear_model: NonnegativeBilinearDensityModel,
) -> None:
    measurement = reconstruction_context["dual"]
    grid = reconstruction_context["grid"]
    assert isinstance(measurement, DualPortFaradayMeasurement)
    assert isinstance(grid, ReconstructionGrid)
    suite = build_reference_morphology_suite(
        grid.y_grid_m,
        grid.z_grid_m,
        peak_column_density_m2=5.0e14,
        radius_y_um=18.0,
        radius_z_um=5.5,
    )
    results = run_linear_readout_morphology_benchmark(
        measurement,
        bilinear_model,
        suite,
        coefficient_upper=1.5,
        regularisation=None,
        initialisation_ridge_strength=1e-8,
        fit_options=DensityFitOptions(irls_iterations=2, max_nfev=150),
    )
    by_name = {result.morphology.name: result for result in results}
    smooth_error = by_name["smooth_reference"].metrics.supported_band_relative_l2_error
    fragmented_error = by_name[
        "fragmented_three_peak_chain"
    ].metrics.supported_band_relative_l2_error
    assert len(results) == 5
    assert all(result.fit.diagnostics.success for result in results)
    assert all(
        result.initialisation.data_jacobian_rank == bilinear_model.parameter_count
        for result in results
    )
    assert smooth_error < 0.25
    assert fragmented_error > smooth_error + 0.20


def test_noise_ensemble_is_reproducible_and_shared_across_candidates(
    reconstruction_context: dict[str, object],
) -> None:
    measurement = reconstruction_context["dual"]
    grid = reconstruction_context["grid"]
    assert isinstance(measurement, DualPortFaradayMeasurement)
    assert isinstance(grid, ReconstructionGrid)
    morphologies = build_reference_morphology_suite(
        grid.y_grid_m,
        grid.z_grid_m,
        peak_column_density_m2=5.0e14,
        radius_y_um=18.0,
        radius_z_um=5.5,
    )[:2]
    first = generate_noisy_observation_ensemble(
        measurement,
        morphologies,
        realizations_per_morphology=2,
        base_seed=20260720,
    )
    second = generate_noisy_observation_ensemble(
        measurement,
        morphologies,
        realizations_per_morphology=2,
        base_seed=20260720,
    )
    assert len(first) == len(second) == 4
    for left, right in zip(first, second, strict=True):
        assert left.seed == right.seed
        assert left.morphology.name == right.morphology.name
        for left_channel, right_channel in zip(left.channels, right.channels, strict=True):
            np.testing.assert_array_equal(left_channel, right_channel)
    replayed = measurement.simulate_channels_from_density(
        first[0].morphology.column_density_m2,
        np.random.default_rng(first[0].seed),
    )
    for stored_channel, replayed_channel in zip(first[0].channels, replayed, strict=True):
        np.testing.assert_array_equal(stored_channel, replayed_channel)
    assert not np.array_equal(first[0].channels[0], first[1].channels[0])


def test_calibration_selection_is_frozen_before_disjoint_held_out_assessment(
    reconstruction_context: dict[str, object],
    bilinear_model: NonnegativeBilinearDensityModel,
) -> None:
    measurement = reconstruction_context["dual"]
    grid = reconstruction_context["grid"]
    assert isinstance(measurement, DualPortFaradayMeasurement)
    assert isinstance(grid, ReconstructionGrid)
    suite = build_reference_morphology_suite(
        grid.y_grid_m,
        grid.z_grid_m,
        peak_column_density_m2=5.0e14,
        radius_y_um=18.0,
        radius_z_um=5.5,
    )
    fine_model = NonnegativeBilinearDensityModel.from_grid(
        y_grid_m=grid.y_grid_m,
        z_grid_m=grid.z_grid_m,
        knot_y_um=np.linspace(-18.0, 18.0, 7),
        knot_z_um=np.linspace(-6.0, 6.0, 5),
        coefficient_scale_m2=5.0e14,
    )
    candidates = (
        ReconstructionCandidate(
            label="coarse_unregularised",
            model=bilinear_model,
            coefficient_upper=1.5,
            regularisation=None,
            fit_options=DensityFitOptions(irls_iterations=2, max_nfev=100),
        ),
        ReconstructionCandidate(
            label="fine_weakly_regularised",
            model=fine_model,
            coefficient_upper=1.5,
            regularisation=build_curvature_regularisation(
                fine_model.knot_y_um,
                fine_model.knot_z_um,
                density_scale_m2=5.0e14,
                weight_um2=1.0,
            ),
            fit_options=DensityFitOptions(irls_iterations=2, max_nfev=120),
        ),
    )
    calibration_observations = generate_noisy_observation_ensemble(
        measurement,
        suite[:2],
        realizations_per_morphology=2,
        base_seed=731,
    )
    initialiser = make_linear_candidate_initialiser(measurement, ridge_strength=1e-8)
    assessments = tuple(
        assess_reconstruction_candidate(
            measurement,
            candidate,
            calibration_observations,
            initialise=initialiser,
        )
        for candidate in candidates
    )
    different_noise = generate_noisy_observation_ensemble(
        measurement,
        suite[:2],
        realizations_per_morphology=2,
        base_seed=732,
    )
    mismatched_assessment = assess_reconstruction_candidate(
        measurement,
        candidates[1],
        different_noise,
        initialise=initialiser,
    )
    with pytest.raises(ValueError, match="same noisy observations"):
        select_and_freeze_candidate(
            [assessments[0], mismatched_assessment],
            minimum_success_fraction=0.75,
        )
    choice = select_and_freeze_candidate(
        assessments,
        minimum_success_fraction=0.75,
        relative_error_tolerance=0.03,
    )
    held_out_observations = generate_noisy_observation_ensemble(
        measurement,
        suite[2:],
        realizations_per_morphology=2,
        base_seed=991,
    )
    held_out = evaluate_frozen_candidate_on_held_out(
        measurement,
        choice,
        held_out_observations,
        initialise=initialiser,
    )
    assert choice.candidate.label in {candidate.label for candidate in candidates}
    assert set(choice.calibration_morphology_names) == {
        "asymmetric_single_cloud",
        "smooth_reference",
    }
    assert held_out.assessment.summary.trial_count == 6
    assert set(held_out.assessment.summary.morphology_names).isdisjoint(
        choice.calibration_morphology_names
    )
    with pytest.raises(ValueError, match="overlap the calibration set"):
        evaluate_frozen_candidate_on_held_out(
            measurement,
            choice,
            calibration_observations,
            initialise=initialiser,
        )
