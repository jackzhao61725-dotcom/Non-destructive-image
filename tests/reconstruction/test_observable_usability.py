from __future__ import annotations

import numpy as np
import pytest

from non_destructive_image.reconstruction import (
    CENTROID_POSITION_ERROR_LIMIT_UM,
    INTEGRATED_RESPONSE_RELATIVE_ERROR_LIMIT,
    MAJOR_RMS_WIDTH_RELATIVE_ERROR_LIMIT,
    ObservableIntegrationSupport,
    ObservableUsabilityEvaluation,
    evaluate_observable_usability,
    extract_density_observables,
)


def _support(*, spacing_um: float = 0.05) -> ObservableIntegrationSupport:
    axis_m = np.arange(-4.0, 4.0 + 0.5 * spacing_um, spacing_um) * 1e-6
    y_grid_m, z_grid_m = np.meshgrid(axis_m, axis_m)
    return ObservableIntegrationSupport(y_grid_m=y_grid_m, z_grid_m=z_grid_m)


def _two_point_cloud(
    support: ObservableIntegrationSupport,
    *,
    half_separation_um: float = 2.0,
    centre_y_um: float = 0.0,
    amplitude: float = 1.0,
) -> np.ndarray:
    density = np.zeros(support.shape, dtype=float)
    y_axis_um = support.y_grid_m[0, :] * 1e6
    z_axis_um = support.z_grid_m[:, 0] * 1e6
    centre_z_index = int(np.argmin(np.abs(z_axis_um)))
    for y_um in (
        centre_y_um - half_separation_um,
        centre_y_um + half_separation_um,
    ):
        y_index = int(np.argmin(np.abs(y_axis_um - y_um)))
        assert y_axis_um[y_index] == pytest.approx(y_um, abs=1e-12)
        density[centre_z_index, y_index] = amplitude
    return density


def test_frozen_usability_limits_are_the_approved_values() -> None:
    assert INTEGRATED_RESPONSE_RELATIVE_ERROR_LIMIT == 0.10
    assert CENTROID_POSITION_ERROR_LIMIT_UM == 0.650
    assert MAJOR_RMS_WIDTH_RELATIVE_ERROR_LIMIT == 0.10


def test_identical_observables_have_zero_coefficients_and_pass_separately() -> None:
    support = _support()
    truth = extract_density_observables(_two_point_cloud(support), support)

    evaluation = evaluate_observable_usability(truth, truth)

    assert isinstance(evaluation, ObservableUsabilityEvaluation)
    assert evaluation.c_a == pytest.approx(0.0)
    assert evaluation.c_r == pytest.approx(0.0)
    assert evaluation.c_w == pytest.approx(0.0)
    assert evaluation.integrated_response_passed
    assert evaluation.centroid_position_passed
    assert evaluation.major_rms_width_passed
    assert not hasattr(evaluation, "overall_passed")


def test_each_coefficient_uses_only_its_declared_physical_error() -> None:
    support = _support()
    truth = extract_density_observables(_two_point_cloud(support), support)

    amplitude_changed = extract_density_observables(
        _two_point_cloud(support, amplitude=1.05),
        support,
    )
    amplitude_evaluation = evaluate_observable_usability(truth, amplitude_changed)
    assert amplitude_evaluation.c_a == pytest.approx(0.5)
    assert amplitude_evaluation.c_r == pytest.approx(0.0, abs=1e-12)
    assert amplitude_evaluation.c_w == pytest.approx(0.0)

    translated = extract_density_observables(
        _two_point_cloud(support, centre_y_um=0.65),
        support,
    )
    translated_evaluation = evaluate_observable_usability(truth, translated)
    assert translated_evaluation.c_a == pytest.approx(0.0)
    assert translated_evaluation.c_r == pytest.approx(1.0)
    assert translated_evaluation.c_w == pytest.approx(0.0)

    widened = extract_density_observables(
        _two_point_cloud(support, half_separation_um=2.2),
        support,
    )
    widened_evaluation = evaluate_observable_usability(truth, widened)
    assert widened_evaluation.c_a == pytest.approx(0.0)
    assert widened_evaluation.c_r == pytest.approx(0.0, abs=1e-12)
    assert widened_evaluation.c_w == pytest.approx(1.0)


def test_pass_decisions_are_independent_and_use_coefficient_at_most_one() -> None:
    support = _support()
    truth = extract_density_observables(_two_point_cloud(support), support)
    reconstructed = extract_density_observables(
        _two_point_cloud(
            support,
            half_separation_um=2.3,
            centre_y_um=0.7,
            amplitude=1.05,
        ),
        support,
    )

    evaluation = evaluate_observable_usability(truth, reconstructed)

    assert evaluation.c_a == pytest.approx(0.5)
    assert evaluation.c_r == pytest.approx(0.7 / 0.650)
    assert evaluation.c_w == pytest.approx(1.5)
    assert evaluation.integrated_response_passed
    assert not evaluation.centroid_position_passed
    assert not evaluation.major_rms_width_passed


def test_undefined_moments_fail_without_hiding_integral_error() -> None:
    support = _support()
    truth = extract_density_observables(_two_point_cloud(support), support)
    blank = extract_density_observables(np.zeros(support.shape), support)

    evaluation = evaluate_observable_usability(truth, blank)

    assert evaluation.c_a == pytest.approx(10.0)
    assert evaluation.c_r is None
    assert evaluation.c_w is None
    assert not evaluation.integrated_response_passed
    assert not evaluation.centroid_position_passed
    assert not evaluation.major_rms_width_passed


def test_evaluation_requires_density_summaries_on_an_identical_support() -> None:
    support = _support()
    truth = extract_density_observables(_two_point_cloud(support), support)
    shifted_support = ObservableIntegrationSupport(
        y_grid_m=support.y_grid_m + 0.1e-6,
        z_grid_m=support.z_grid_m,
        cell_area_m2=support.cell_area_m2,
    )
    shifted = extract_density_observables(
        _two_point_cloud(shifted_support),
        shifted_support,
    )

    with pytest.raises(ValueError, match="identical coordinate grids"):
        evaluate_observable_usability(truth, shifted)
    with pytest.raises(TypeError, match="density observable summaries"):
        evaluate_observable_usability(truth, object())  # type: ignore[arg-type]


def test_blank_truth_makes_all_coefficients_unavailable() -> None:
    support = _support()
    blank = extract_density_observables(np.zeros(support.shape), support)

    evaluation = evaluate_observable_usability(blank, blank)

    assert evaluation.c_a is None
    assert evaluation.c_r is None
    assert evaluation.c_w is None
    assert not evaluation.integrated_response_passed
    assert not evaluation.centroid_position_passed
    assert not evaluation.major_rms_width_passed


def test_evaluation_dataclass_rejects_inconsistent_manual_results() -> None:
    with pytest.raises(ValueError, match="pass flags"):
        ObservableUsabilityEvaluation(
            c_a=0.5,
            c_r=None,
            c_w=2.0,
            integrated_response_passed=False,
            centroid_position_passed=False,
            major_rms_width_passed=False,
        )
    with pytest.raises(ValueError, match="finite and non-negative"):
        ObservableUsabilityEvaluation(
            c_a=np.nan,
            c_r=None,
            c_w=None,
            integrated_response_passed=False,
            centroid_position_passed=False,
            major_rms_width_passed=False,
        )
