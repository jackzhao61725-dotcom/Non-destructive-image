"""Evidence contracts for deciding what a reconstruction can support.

The comparison implemented here is deliberately narrower than a detection
claim.  It asks whether one fitted forward prediction describes the calibrated
camera counts better than the *declared zero-density forward model*, using the
same Gaussian quasi-deviance as the existing inverse fit.  That statistic is
not a likelihood-ratio p value: the alternative is regularised, coefficients
are bounded and the Gaussian treatment is only an approximation to the mixed
Poisson--Gaussian camera likelihood.

An empirical upper-tail probability is produced only when an independent set
of matched experimental blank acquisitions has been processed through the
complete, pre-frozen reconstruction pipeline.  Synthetic blanks remain a
development diagnostic and never create an experimental acceptance threshold.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .density_fit import DensityFitResult
from .measurements import DifferentiableDensityMeasurement
from .noise import gaussian_quasi_deviance


ReferenceOrigin = Literal[
    "synthetic_development",
    "independent_experimental_blank",
]
EvidenceLevel = Literal[
    "model_only",
    "development_rank_only",
    "matched_blank_empirical",
]
TargetOrigin = Literal[
    "synthetic_development",
    "experimental_observation",
]
ResponseScaleStatus = Literal[
    "calibrated",
    "provisional",
    "illustrative_uncalibrated",
]


def _nonempty_text(value: str, label: str) -> str:
    text = str(value).strip()
    if not text:
        raise ValueError(f"{label} cannot be empty")
    return text


def _frozen_vector(values: ArrayLike, label: str) -> NDArray[np.floating]:
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or array.size == 0:
        raise ValueError(f"{label} must be a non-empty one-dimensional array")
    if np.any(~np.isfinite(array)):
        raise ValueError(f"{label} must contain only finite values")
    copied = array.copy()
    return np.frombuffer(copied.tobytes(order="C"), dtype=float)


@dataclass(frozen=True)
class NullReferenceDistribution:
    """Completed blank-pipeline statistics used to contextualise one target.

    Each entry in ``delta_quasi_deviance`` must come from a complete fit of one
    blank acquisition through the same alternative reconstruction used for the
    target.  ``acquisition_ids`` contains the successful blank fits;
    ``failed_acquisition_ids`` records every attempted blank that did not yield
    a statistic.  Failed blanks make the set ineligible for calibrated rank
    evidence in this implementation, because silently dropping them can bias
    the empirical tail.
    """

    delta_quasi_deviance: NDArray[np.floating] = field(repr=False, compare=False)
    origin: ReferenceOrigin
    acquisition_ids: tuple[str, ...]
    attempted_count: int
    failed_acquisition_ids: tuple[str, ...]
    pipeline_fingerprint: str
    condition_fingerprint: str
    independent_of_target: bool
    pipeline_frozen_before_target: bool

    def __post_init__(self) -> None:
        if self.origin not in {
            "synthetic_development",
            "independent_experimental_blank",
        }:
            raise ValueError("unknown null-reference origin")
        values = _frozen_vector(
            self.delta_quasi_deviance,
            "blank delta-quasi-deviance distribution",
        )
        acquisition_ids = tuple(
            _nonempty_text(value, "blank acquisition_id")
            for value in self.acquisition_ids
        )
        failed_ids = tuple(
            _nonempty_text(value, "failed blank acquisition_id")
            for value in self.failed_acquisition_ids
        )
        if len(acquisition_ids) != values.size:
            raise ValueError(
                "blank acquisition_ids must correspond one-to-one with statistics"
            )
        if len(set(acquisition_ids + failed_ids)) != len(acquisition_ids) + len(
            failed_ids
        ):
            raise ValueError("blank acquisition IDs must be unique")
        if not isinstance(self.attempted_count, int) or isinstance(
            self.attempted_count,
            bool,
        ):
            raise TypeError("attempted_count must be an integer")
        if self.attempted_count != len(acquisition_ids) + len(failed_ids):
            raise ValueError(
                "attempted_count must equal successful plus failed blank acquisitions"
            )
        object.__setattr__(self, "delta_quasi_deviance", values)
        object.__setattr__(self, "acquisition_ids", acquisition_ids)
        object.__setattr__(self, "failed_acquisition_ids", failed_ids)
        object.__setattr__(
            self,
            "pipeline_fingerprint",
            _nonempty_text(self.pipeline_fingerprint, "pipeline_fingerprint"),
        )
        object.__setattr__(
            self,
            "condition_fingerprint",
            _nonempty_text(self.condition_fingerprint, "condition_fingerprint"),
        )


@dataclass(frozen=True)
class NullSignalEvidence:
    """Zero-density comparison with an explicitly bounded evidence level."""

    statistic_name: Literal["delta_gaussian_quasi_deviance"]
    null_model: Literal["declared_zero_density_forward_model"]
    target_origin: TargetOrigin
    pipeline_fingerprint: str
    condition_fingerprint: str
    null_quasi_deviance: float
    fitted_quasi_deviance: float
    delta_quasi_deviance: float
    alternative_is_regularised: bool
    evidence_level: EvidenceLevel
    reference_origin: ReferenceOrigin | None
    reference_count: int
    failed_reference_count: int
    exceedance_count: int | None
    empirical_upper_tail_probability: float | None
    tail_probability_resolution: float | None
    predeclared_alpha: float | None
    crosses_predeclared_level: bool | None
    assumptions: tuple[str, ...]

    def __post_init__(self) -> None:
        for name in (
            "null_quasi_deviance",
            "fitted_quasi_deviance",
            "delta_quasi_deviance",
        ):
            if not np.isfinite(float(getattr(self, name))):
                raise ValueError(f"{name} must be finite")
        if self.reference_count < 0 or self.failed_reference_count < 0:
            raise ValueError("reference counts cannot be negative")
        object.__setattr__(
            self,
            "pipeline_fingerprint",
            _nonempty_text(self.pipeline_fingerprint, "pipeline_fingerprint"),
        )
        object.__setattr__(
            self,
            "condition_fingerprint",
            _nonempty_text(self.condition_fingerprint, "condition_fingerprint"),
        )
        object.__setattr__(self, "assumptions", tuple(self.assumptions))


@dataclass(frozen=True)
class ResponseScaleDeclaration:
    """Provenance and calibration status of the Faraday response scale."""

    kappa_f: float
    status: ResponseScaleStatus
    source: str
    calibration_id: str | None = None
    standard_uncertainty: float | None = None

    def __post_init__(self) -> None:
        kappa = float(self.kappa_f)
        if not np.isfinite(kappa):
            raise ValueError("kappa_f must be finite")
        if self.status not in {
            "calibrated",
            "provisional",
            "illustrative_uncalibrated",
        }:
            raise ValueError("unknown response-scale status")
        source = _nonempty_text(self.source, "response-scale source")
        calibration_id = (
            None
            if self.calibration_id is None
            else _nonempty_text(self.calibration_id, "calibration_id")
        )
        if self.status == "calibrated" and calibration_id is None:
            raise ValueError("a calibrated response scale requires calibration_id")
        if self.status != "calibrated" and calibration_id is not None:
            raise ValueError(
                "only a calibrated response scale may carry calibration_id"
            )
        uncertainty = (
            None
            if self.standard_uncertainty is None
            else float(self.standard_uncertainty)
        )
        if uncertainty is not None and (
            not np.isfinite(uncertainty) or uncertainty < 0.0
        ):
            raise ValueError("response-scale standard_uncertainty cannot be negative")
        object.__setattr__(self, "kappa_f", kappa)
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "calibration_id", calibration_id)
        object.__setattr__(self, "standard_uncertainty", uncertainty)

    def assert_matches_measurement(
        self,
        measurement: DifferentiableDensityMeasurement,
    ) -> None:
        """Reject a declaration that does not describe the fitted operator."""

        response = getattr(measurement, "response", None)
        if response is None or not hasattr(response, "kappa_f"):
            raise TypeError("measurement does not expose a Faraday response contract")
        if not np.isclose(
            self.kappa_f,
            float(response.kappa_f),
            rtol=1e-12,
            atol=1e-15,
        ):
            raise ValueError(
                "response-scale kappa_f does not match the measurement operator"
            )


def compare_to_zero_density(
    measurement: DifferentiableDensityMeasurement,
    observed_channels: tuple[ArrayLike, ...],
    fit_result: DensityFitResult,
    *,
    pipeline_fingerprint: str,
    condition_fingerprint: str,
    target_origin: TargetOrigin,
    reference: NullReferenceDistribution | None = None,
    predeclared_alpha: float | None = None,
) -> NullSignalEvidence:
    """Compare a complete fitted prediction with the declared zero-density model.

    The fitted statistic contains only the camera-data quasi-deviance; the
    regularisation objective is intentionally excluded.  ``delta`` is retained
    with its sign and is never clipped.  A positive value favours the fitted
    prediction under this score, but does not by itself establish a physical
    signal.
    """

    pipeline = _nonempty_text(pipeline_fingerprint, "pipeline_fingerprint")
    condition = _nonempty_text(condition_fingerprint, "condition_fingerprint")
    if target_origin not in {"synthetic_development", "experimental_observation"}:
        raise ValueError("unknown target_origin")
    if predeclared_alpha is not None:
        alpha = float(predeclared_alpha)
        if not np.isfinite(alpha) or not 0.0 < alpha < 1.0:
            raise ValueError("predeclared_alpha must lie strictly between zero and one")
    else:
        alpha = None

    observed = measurement.flatten_observed(*observed_channels)
    zero_density = np.zeros_like(measurement.grid.y_grid_m, dtype=float)
    zero_channels = measurement.expected_channels_from_density(zero_density)
    null_prediction = measurement.flatten_observed(*zero_channels)
    fitted_channels = measurement.expected_channels_from_density(
        fit_result.column_density_m2
    )
    if len(fitted_channels) != len(fit_result.predicted_channels) or any(
        not np.allclose(recomputed, stored, rtol=1e-12, atol=1e-9)
        for recomputed, stored in zip(
            fitted_channels,
            fit_result.predicted_channels,
            strict=True,
        )
    ):
        raise ValueError(
            "fit stored prediction does not match its density under the measurement operator"
        )
    fitted_prediction = measurement.flatten_observed(*fitted_channels)
    if observed.shape != fitted_prediction.shape:
        raise ValueError("fit prediction does not match the observed count vector")
    null_quasi_deviance = gaussian_quasi_deviance(
        observed,
        null_prediction,
        read_noise_electrons=measurement.read_noise_electrons,
    )
    fitted_quasi_deviance = gaussian_quasi_deviance(
        observed,
        fitted_prediction,
        read_noise_electrons=measurement.read_noise_electrons,
    )
    delta = null_quasi_deviance - fitted_quasi_deviance
    if not np.isclose(
        fitted_quasi_deviance,
        fit_result.diagnostics.quasi_deviance,
        rtol=1e-10,
        atol=1e-8,
    ):
        raise ValueError(
            "fit diagnostics quasi-deviance does not match the raw-channel prediction"
        )

    level: EvidenceLevel = "model_only"
    reference_origin: ReferenceOrigin | None = None
    reference_count = 0
    failed_reference_count = 0
    exceedance_count: int | None = None
    empirical_tail: float | None = None
    tail_resolution: float | None = None
    crosses: bool | None = None
    assumptions = [
        "The null is zero column density under the declared forward operator, not a measured instrumental background model.",
        "The statistic uses the Gaussian approximation to Poisson photoelectron and Gaussian read noise.",
        "The alternative fit may be bounded and regularised, so Wilks or chi-square likelihood-ratio interpretations do not apply.",
        "Unmodelled analyser leakage, offsets, gain imbalance and optical mismatch are not represented by this comparison.",
    ]
    if reference is not None:
        reference_origin = reference.origin
        reference_count = int(reference.delta_quasi_deviance.size)
        failed_reference_count = len(reference.failed_acquisition_ids)
        if reference.pipeline_fingerprint != pipeline:
            assumptions.append(
                "The supplied blank reference used a different pipeline fingerprint and was not used for empirical calibration."
            )
        elif reference.condition_fingerprint != condition:
            assumptions.append(
                "The supplied blank reference used a different acquisition-condition fingerprint and was not used for empirical calibration."
            )
        else:
            exceedance_count = int(
                np.count_nonzero(reference.delta_quasi_deviance >= delta)
            )
            if reference.origin == "synthetic_development":
                level = "development_rank_only"
                assumptions.append(
                    "Synthetic blank ranks are development diagnostics only and do not define an experimental threshold or p value."
                )
                if failed_reference_count:
                    assumptions.append(
                        "Some synthetic blank reconstructions failed, so the reported development rank is incomplete."
                    )
            elif (
                target_origin == "experimental_observation"
                and reference.independent_of_target
                and reference.pipeline_frozen_before_target
                and failed_reference_count == 0
            ):
                level = "matched_blank_empirical"
                denominator = reference_count + 1
                empirical_tail = (1.0 + exceedance_count) / denominator
                tail_resolution = 1.0 / denominator
                crosses = None if alpha is None else bool(empirical_tail <= alpha)
                assumptions.append(
                    "The empirical upper tail is calibrated by independent matched experimental blanks processed through the complete pre-frozen pipeline."
                )
            else:
                assumptions.append(
                    "The experimental blank set was not eligible for calibrated rank evidence because independence, pipeline freezing, target origin or complete blank-fit success was not satisfied."
                )
    if alpha is not None and level != "matched_blank_empirical":
        raise ValueError(
            "predeclared_alpha requires eligible matched experimental blank evidence"
        )

    return NullSignalEvidence(
        statistic_name="delta_gaussian_quasi_deviance",
        null_model="declared_zero_density_forward_model",
        target_origin=target_origin,
        pipeline_fingerprint=pipeline,
        condition_fingerprint=condition,
        null_quasi_deviance=float(null_quasi_deviance),
        fitted_quasi_deviance=float(fitted_quasi_deviance),
        delta_quasi_deviance=float(delta),
        alternative_is_regularised=bool(
            fit_result.diagnostics.curvature_weight_um2 > 0.0
        ),
        evidence_level=level,
        reference_origin=reference_origin,
        reference_count=reference_count,
        failed_reference_count=failed_reference_count,
        exceedance_count=exceedance_count,
        empirical_upper_tail_probability=empirical_tail,
        tail_probability_resolution=tail_resolution,
        predeclared_alpha=alpha,
        crosses_predeclared_level=crosses,
        assumptions=tuple(assumptions),
    )
