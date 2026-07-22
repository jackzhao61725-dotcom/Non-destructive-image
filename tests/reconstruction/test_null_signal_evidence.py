"""Focused tests for bounded zero-density evidence claims."""

from __future__ import annotations

import numpy as np
import pytest

from non_destructive_image.reconstruction.contracts import (
    DetectorContract,
    FaradayResponseContract,
    ReconstructionGrid,
)
from non_destructive_image.reconstruction.density_fit import (
    DensityFitDiagnostics,
    DensityFitResult,
)
from non_destructive_image.reconstruction.evidence import (
    NullReferenceDistribution,
    ResponseScaleDeclaration,
    compare_to_zero_density,
)
from non_destructive_image.reconstruction.measurements import (
    DualPortFaradayMeasurement,
)
from non_destructive_image.reconstruction.noise import gaussian_quasi_deviance


@pytest.fixture
def measurement() -> DualPortFaradayMeasurement:
    coordinate = np.linspace(-4e-6, 4e-6, 8)
    y_grid, z_grid = np.meshgrid(coordinate, coordinate)
    grid = ReconstructionGrid.from_arrays(
        y_grid_m=y_grid,
        z_grid_m=z_grid,
        pupil=np.ones((8, 8), dtype=complex),
        bin_size=2,
        roi_mask=np.ones((4, 4), dtype=bool),
    )
    return DualPortFaradayMeasurement(
        grid=grid,
        detector=DetectorContract(1000.0, 3.0),
        response=FaradayResponseContract(1.0e-14, 1.0),
    )


def _density(measurement: DualPortFaradayMeasurement) -> np.ndarray:
    radius2 = measurement.grid.y_grid_m**2 + measurement.grid.z_grid_m**2
    return 5.0e12 * np.exp(-radius2 / (2.0 * (2.0e-6) ** 2))


def _fit_result(
    measurement: DualPortFaradayMeasurement,
    density: np.ndarray,
    observed_channels: tuple[np.ndarray, ...],
    *,
    regularised: bool = False,
) -> DensityFitResult:
    predicted = measurement.expected_channels_from_density(density)
    observed = measurement.flatten_observed(*observed_channels)
    predicted_vector = measurement.flatten_observed(*predicted)
    quasi_deviance = gaussian_quasi_deviance(
        observed,
        predicted_vector,
        read_noise_electrons=measurement.read_noise_electrons,
    )
    residual = observed - predicted_vector
    return DensityFitResult(
        coefficients=np.asarray([1.0]),
        column_density_m2=np.asarray(density),
        predicted_channels=tuple(np.asarray(channel) for channel in predicted),
        diagnostics=DensityFitDiagnostics(
            success=True,
            message="constructed test result",
            quasi_deviance=quasi_deviance,
            weighted_chi_square=float(np.dot(residual, residual)),
            reduced_chi_square=0.0,
            degrees_of_freedom=max(observed.size - 1, 1),
            nfev=1,
            irls_iterations=1,
            data_jacobian_rank=1,
            data_jacobian_condition=1.0,
            curvature_weight_um2=1.0 if regularised else 0.0,
            curvature_seminorm_per_um2=0.0,
            regularisation_objective=0.0,
            regularisation_density_scale_m2=None,
            regularisation_boundary_policy=None,
            regularisation_axis_weights=None,
            active_lower_coefficients=0,
            active_upper_coefficients=0,
            whitened_residual_vector=np.zeros_like(observed),
        ),
    )


def _compare(
    measurement: DualPortFaradayMeasurement,
    observed_channels: tuple[np.ndarray, ...],
    fit: DensityFitResult,
    **kwargs: object,
):
    return compare_to_zero_density(
        measurement,
        observed_channels,
        fit,
        pipeline_fingerprint="pipeline-v1",
        condition_fingerprint="condition-A",
        target_origin=kwargs.pop("target_origin", "synthetic_development"),
        **kwargs,
    )


def test_delta_retains_zero_positive_and_negative_signs(
    measurement: DualPortFaradayMeasurement,
) -> None:
    zero = np.zeros_like(measurement.grid.y_grid_m)
    injected = _density(measurement)
    zero_channels = measurement.expected_channels_from_density(zero)
    injected_channels = measurement.expected_channels_from_density(injected)

    zero_evidence = _compare(
        measurement,
        zero_channels,
        _fit_result(measurement, zero, zero_channels),
    )
    injected_evidence = _compare(
        measurement,
        injected_channels,
        _fit_result(measurement, injected, injected_channels, regularised=True),
    )
    wrong_fit_evidence = _compare(
        measurement,
        zero_channels,
        _fit_result(measurement, injected, zero_channels),
    )

    assert zero_evidence.delta_quasi_deviance == pytest.approx(0.0)
    assert injected_evidence.delta_quasi_deviance > 0.0
    assert injected_evidence.alternative_is_regularised
    assert wrong_fit_evidence.delta_quasi_deviance < 0.0


def test_synthetic_blanks_provide_rank_only_and_never_a_p_value(
    measurement: DualPortFaradayMeasurement,
) -> None:
    density = _density(measurement)
    observed = measurement.expected_channels_from_density(density)
    fit = _fit_result(measurement, density, observed)
    reference = NullReferenceDistribution(
        delta_quasi_deviance=np.asarray([0.1, 0.2, 0.3]),
        origin="synthetic_development",
        acquisition_ids=("s0", "s1", "s2"),
        attempted_count=3,
        failed_acquisition_ids=(),
        pipeline_fingerprint="pipeline-v1",
        condition_fingerprint="condition-A",
        independent_of_target=True,
        pipeline_frozen_before_target=True,
    )
    evidence = _compare(measurement, observed, fit, reference=reference)

    assert evidence.evidence_level == "development_rank_only"
    assert evidence.reference_origin == "synthetic_development"
    assert evidence.exceedance_count is not None
    assert evidence.empirical_upper_tail_probability is None
    assert evidence.tail_probability_resolution is None
    assert evidence.crosses_predeclared_level is None
    assert any("do not define an experimental threshold or p value" in item for item in evidence.assumptions)
    with pytest.raises(ValueError, match="predeclared_alpha requires"):
        _compare(
            measurement,
            observed,
            fit,
            reference=reference,
            predeclared_alpha=0.05,
        )


def test_matched_experimental_blank_rank_uses_plus_one_formula(
    measurement: DualPortFaradayMeasurement,
) -> None:
    density = _density(measurement)
    observed = measurement.expected_channels_from_density(density)
    fit = _fit_result(measurement, density, observed)
    baseline = _compare(measurement, observed, fit)
    blank_values = baseline.delta_quasi_deviance - np.asarray([1.0, 2.0, 3.0, 4.0])
    reference = NullReferenceDistribution(
        delta_quasi_deviance=blank_values,
        origin="independent_experimental_blank",
        acquisition_ids=("b0", "b1", "b2", "b3"),
        attempted_count=4,
        failed_acquisition_ids=(),
        pipeline_fingerprint="pipeline-v1",
        condition_fingerprint="condition-A",
        independent_of_target=True,
        pipeline_frozen_before_target=True,
    )
    evidence = _compare(
        measurement,
        observed,
        fit,
        target_origin="experimental_observation",
        reference=reference,
        predeclared_alpha=0.25,
    )

    assert evidence.evidence_level == "matched_blank_empirical"
    assert evidence.exceedance_count == 0
    assert evidence.empirical_upper_tail_probability == pytest.approx(1.0 / 5.0)
    assert evidence.tail_probability_resolution == pytest.approx(1.0 / 5.0)
    assert evidence.crosses_predeclared_level is True


@pytest.mark.parametrize(
    ("change", "expected_fragment"),
    [
        ({"pipeline_fingerprint": "other-pipeline"}, "different pipeline"),
        ({"condition_fingerprint": "other-condition"}, "different acquisition-condition"),
        ({"independent_of_target": False}, "not eligible"),
        ({"pipeline_frozen_before_target": False}, "not eligible"),
        (
            {
                "delta_quasi_deviance": np.asarray([0.1, 0.2]),
                "acquisition_ids": ("b0", "b1"),
                "attempted_count": 3,
                "failed_acquisition_ids": ("failed-b2",),
            },
            "not eligible",
        ),
    ],
)
def test_ineligible_experimental_blanks_downgrade_without_probability(
    measurement: DualPortFaradayMeasurement,
    change: dict[str, object],
    expected_fragment: str,
) -> None:
    density = _density(measurement)
    observed = measurement.expected_channels_from_density(density)
    fit = _fit_result(measurement, density, observed)
    values: dict[str, object] = {
        "delta_quasi_deviance": np.asarray([0.1, 0.2, 0.3]),
        "origin": "independent_experimental_blank",
        "acquisition_ids": ("b0", "b1", "b2"),
        "attempted_count": 3,
        "failed_acquisition_ids": (),
        "pipeline_fingerprint": "pipeline-v1",
        "condition_fingerprint": "condition-A",
        "independent_of_target": True,
        "pipeline_frozen_before_target": True,
    }
    values.update(change)
    reference = NullReferenceDistribution(**values)  # type: ignore[arg-type]
    evidence = _compare(
        measurement,
        observed,
        fit,
        target_origin="experimental_observation",
        reference=reference,
    )

    assert evidence.evidence_level == "model_only"
    assert evidence.empirical_upper_tail_probability is None
    assert evidence.crosses_predeclared_level is None
    assert any(expected_fragment in item for item in evidence.assumptions)


def test_response_scale_declaration_must_match_fitted_measurement(
    measurement: DualPortFaradayMeasurement,
) -> None:
    ResponseScaleDeclaration(
        kappa_f=1.0,
        status="illustrative_uncalibrated",
        source="declared simulation scale",
    ).assert_matches_measurement(measurement)

    with pytest.raises(ValueError, match="does not match"):
        ResponseScaleDeclaration(
            kappa_f=0.5,
            status="provisional",
            source="provisional response",
        ).assert_matches_measurement(measurement)

    with pytest.raises(ValueError, match="requires calibration_id"):
        ResponseScaleDeclaration(
            kappa_f=1.0,
            status="calibrated",
            source="paired measurement",
        )
