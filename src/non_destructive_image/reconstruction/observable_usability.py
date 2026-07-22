"""Frozen observable-specific usability scores for synthetic reconstruction.

The scores in this module answer a deliberately narrow question: whether a
truth-known synthetic reconstruction recovers each declared physical
observable within its fixed tolerance.  They do not alter the inverse solver,
and they are reported separately rather than combined into an overall image
quality score.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .observables import DensityObservableSummary


INTEGRATED_RESPONSE_RELATIVE_ERROR_LIMIT = 0.10
CENTROID_POSITION_ERROR_LIMIT_UM = 0.650
MAJOR_RMS_WIDTH_RELATIVE_ERROR_LIMIT = 0.10


def _score_passed(score: float | None) -> bool:
    return score is not None and score <= 1.0


@dataclass(frozen=True)
class ObservableUsabilityEvaluation:
    """Separate usability coefficients and pass decisions for three observables.

    ``c_a`` is the absolute relative supported-integral error divided by
    ``0.10``.  ``c_r`` is the Euclidean centroid error in micrometres divided
    by ``0.650``.  ``c_w`` is the absolute relative major-rms-width error
    divided by ``0.10``.  A coefficient at or below one passes.

    A coefficient is ``None`` and its pass flag is false when that observable
    is undefined in either input summary.  There is intentionally no overall
    score or overall pass decision.
    """

    c_a: float | None
    c_r: float | None
    c_w: float | None
    integrated_response_passed: bool
    centroid_position_passed: bool
    major_rms_width_passed: bool

    def __post_init__(self) -> None:
        for name in ("c_a", "c_r", "c_w"):
            value = getattr(self, name)
            if value is not None and (not np.isfinite(value) or value < 0.0):
                raise ValueError(f"{name} must be finite and non-negative when defined")
        expected_passes = (
            _score_passed(self.c_a),
            _score_passed(self.c_r),
            _score_passed(self.c_w),
        )
        declared_passes = (
            self.integrated_response_passed,
            self.centroid_position_passed,
            self.major_rms_width_passed,
        )
        if declared_passes != expected_passes:
            raise ValueError(
                "observable pass flags must match their usability coefficients"
            )


def evaluate_observable_usability(
    truth: DensityObservableSummary,
    reconstructed: DensityObservableSummary,
) -> ObservableUsabilityEvaluation:
    """Evaluate the three frozen usability coefficients on one fixed support.

    ``truth`` must be the unfiltered physical truth evaluated on exactly the
    same coordinate grid, cell areas and support mask as ``reconstructed``.
    Undefined truth or reconstructed moments produce an unavailable score and
    a failed decision for that observable, while a blank reconstruction still
    receives its well-defined supported-integral score.
    """

    if not isinstance(truth, DensityObservableSummary) or not isinstance(
        reconstructed,
        DensityObservableSummary,
    ):
        raise TypeError("usability evaluation requires density observable summaries")
    if not truth.integration_support.is_identical_to(
        reconstructed.integration_support
    ):
        raise ValueError(
            "usability evaluation requires identical coordinate grids, cell areas, "
            "and integration support masks"
        )

    if truth.integrated_response > 0.0:
        c_a = abs(
            reconstructed.integrated_response / truth.integrated_response - 1.0
        ) / INTEGRATED_RESPONSE_RELATIVE_ERROR_LIMIT
    else:
        c_a = None

    if truth.centroid_m is not None and reconstructed.centroid_m is not None:
        centroid_error_um = float(
            np.linalg.norm(reconstructed.centroid_m - truth.centroid_m) * 1e6
        )
        c_r = centroid_error_um / CENTROID_POSITION_ERROR_LIMIT_UM
    else:
        c_r = None

    if (
        truth.major_rms_width_m is not None
        and truth.major_rms_width_m > 0.0
        and reconstructed.major_rms_width_m is not None
    ):
        c_w = abs(
            reconstructed.major_rms_width_m / truth.major_rms_width_m - 1.0
        ) / MAJOR_RMS_WIDTH_RELATIVE_ERROR_LIMIT
    else:
        c_w = None

    return ObservableUsabilityEvaluation(
        c_a=c_a,
        c_r=c_r,
        c_w=c_w,
        integrated_response_passed=_score_passed(c_a),
        centroid_position_passed=_score_passed(c_r),
        major_rms_width_passed=_score_passed(c_w),
    )


__all__ = [
    "CENTROID_POSITION_ERROR_LIMIT_UM",
    "INTEGRATED_RESPONSE_RELATIVE_ERROR_LIMIT",
    "MAJOR_RMS_WIDTH_RELATIVE_ERROR_LIMIT",
    "ObservableUsabilityEvaluation",
    "evaluate_observable_usability",
]
