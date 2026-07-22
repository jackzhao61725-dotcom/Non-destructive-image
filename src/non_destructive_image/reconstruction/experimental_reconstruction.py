"""End-to-end reconstruction of one calibrated non-destructive exposure.

This module joins the strict experimental camera adapter to the existing
shape-flexible inverse and credibility layers.  The result is an evidence
bundle, not a denoised picture: its formal output is a fixed-support physical
observable estimate plus input hashes, fit diagnostics, the zero-density
comparison, local data/prior support and the explicit status of the Faraday
response scale.  It never accepts or stores a truth-density map.  Latent
density artifacts are retained only by explicit diagnostic opt-in.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Literal, Mapping

import numpy as np
from numpy.typing import ArrayLike

from .credibility import (
    LocalCredibilityDiagnostics,
    ParametricBootstrapResult,
    ReconstructionStabilitySummary,
    analyse_local_credibility,
    parametric_bootstrap_reconstruction,
    summarise_reconstruction_stability,
)
from .density_fit import (
    DensityFitDiagnostics,
    DensityFitOptions,
    DensityFitResult,
    fit_nonnegative_basis_density,
)
from .evidence import (
    NullReferenceDistribution,
    NullSignalEvidence,
    ResponseScaleDeclaration,
    compare_to_zero_density,
)
from .experimental_adapter import CameraCalibration, bind_non_destructive_exposure
from .measurements import DifferentiableDensityMeasurement
from .object_models import NonnegativeBilinearDensityModel
from .observables import (
    DensityObservableSummary,
    ObservableIntegrationSupport,
    extract_density_observables,
)
from .observations import NonDestructiveSequence
from .regularisation import CurvatureRegularisation


AbsoluteDensityStatus = Literal[
    "calibrated",
    "conditional_on_assumed_kappa_f",
]
BootstrapStatus = Literal["not_requested", "complete", "partial", "failed"]


def _nonempty_text(value: str, label: str) -> str:
    text = str(value).strip()
    if not text:
        raise ValueError(f"{label} cannot be empty")
    return text


def _array_sha256(values: ArrayLike) -> str:
    array = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
    digest = hashlib.sha256()
    digest.update(str(array.shape).encode("ascii"))
    digest.update(array.dtype.str.encode("ascii"))
    digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def _observable_support_sha256(support: ObservableIntegrationSupport) -> str:
    digest = hashlib.sha256()
    digest.update(b"observable-integration-support-v1")
    for label, values in (
        ("y_grid_m", support.y_grid_m),
        ("z_grid_m", support.z_grid_m),
        ("cell_area_m2", support.cell_area_m2),
        ("support_mask", support.support_mask),
    ):
        array = np.ascontiguousarray(values)
        digest.update(label.encode("ascii"))
        digest.update(str(array.shape).encode("ascii"))
        digest.update(array.dtype.str.encode("ascii"))
        digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def _validate_observable_integration_support(
    measurement: DifferentiableDensityMeasurement,
    model: NonnegativeBilinearDensityModel,
    support: ObservableIntegrationSupport,
) -> None:
    if not isinstance(support, ObservableIntegrationSupport):
        raise TypeError(
            "observable_integration_support must be ObservableIntegrationSupport"
        )
    grid = measurement.grid
    if not np.array_equal(model.y_grid_m, grid.y_grid_m) or not np.array_equal(
        model.z_grid_m,
        grid.z_grid_m,
    ):
        raise ValueError("density model coordinates must match the measurement grid")
    if not np.array_equal(support.y_grid_m, grid.y_grid_m) or not np.array_equal(
        support.z_grid_m,
        grid.z_grid_m,
    ):
        raise ValueError(
            "observable integration coordinates must exactly match the measurement grid"
        )
    if np.any(support.support_mask & ~np.asarray(model.support_mask, dtype=bool)):
        raise ValueError(
            "observable integration mask cannot extend beyond latent model support"
        )


@dataclass(frozen=True)
class BootstrapRequest:
    """Explicit request for conditional detector-noise resampling."""

    draws: int
    seed: int
    confidence_level: float = 0.68
    fit_options: DensityFitOptions | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.draws, int) or isinstance(self.draws, bool):
            raise TypeError("bootstrap draws must be an integer")
        if self.draws <= 0:
            raise ValueError("bootstrap draws must be positive")
        if not isinstance(self.seed, int) or isinstance(self.seed, bool):
            raise TypeError("bootstrap seed must be an integer")
        confidence = float(self.confidence_level)
        if not np.isfinite(confidence) or not 0.0 < confidence < 1.0:
            raise ValueError("bootstrap confidence_level must lie between zero and one")
        object.__setattr__(self, "confidence_level", confidence)


@dataclass(frozen=True)
class BootstrapAssessment:
    """Outcome of an optional conditional detector-noise bootstrap."""

    status: BootstrapStatus
    result: ParametricBootstrapResult | None
    message: str | None

    def __post_init__(self) -> None:
        carries_result = self.status in ("complete", "partial")
        if carries_result != (self.result is not None):
            raise ValueError(
                "only complete or partial bootstrap assessments may carry a result"
            )
        if self.status in ("partial", "failed") and not self.message:
            raise ValueError("partial or failed bootstrap assessment requires a message")
        if self.result is not None:
            fully_converged = (
                self.result.successful_draws == self.result.requested_draws
            )
            if self.status == "complete" and not fully_converged:
                raise ValueError("complete bootstrap requires every requested draw")
            if self.status == "partial" and fully_converged:
                raise ValueError("partial bootstrap requires at least one failed draw")


@dataclass(frozen=True)
class ExperimentalReconstructionResult:
    """One exposure plus the evidence needed to interpret its reconstruction."""

    context_id: str
    sequence_id: str
    exposure_index: int
    timestamp_s: float
    readout: str
    channel_names: tuple[str, ...]
    input_channel_sha256: tuple[tuple[str, str], ...]
    pipeline_fingerprint: str
    condition_fingerprint: str
    observable_support_sha256: str
    response_scale: ResponseScaleDeclaration
    absolute_density_status: AbsoluteDensityStatus
    fit_diagnostics: DensityFitDiagnostics
    observables: DensityObservableSummary
    latent_fit: DensityFitResult | None
    local_credibility: LocalCredibilityDiagnostics
    null_evidence: NullSignalEvidence
    bootstrap: BootstrapAssessment
    stability: ReconstructionStabilitySummary | None
    interpretation_limits: tuple[str, ...]

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "context_id",
            _nonempty_text(self.context_id, "context_id"),
        )
        object.__setattr__(
            self,
            "sequence_id",
            _nonempty_text(self.sequence_id, "sequence_id"),
        )
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
        if len(self.observable_support_sha256) != 64 or any(
            character not in "0123456789abcdef"
            for character in self.observable_support_sha256
        ):
            raise ValueError(
                "observable support hash must be a lowercase SHA-256 digest"
            )
        if tuple(name for name, _ in self.input_channel_sha256) != self.channel_names:
            raise ValueError("input hashes must follow the declared channel order")
        for _, digest in self.input_channel_sha256:
            if len(digest) != 64 or any(
                character not in "0123456789abcdef" for character in digest
            ):
                raise ValueError(
                    "input channel hashes must be lowercase SHA-256 digests"
                )
        object.__setattr__(
            self,
            "interpretation_limits",
            tuple(self.interpretation_limits),
        )


def _exposure_by_index(sequence: NonDestructiveSequence, exposure_index: int):
    matches = [
        exposure
        for exposure in sequence.exposures
        if exposure.exposure_index == exposure_index
    ]
    if len(matches) != 1:
        raise ValueError(
            f"sequence must contain exactly one exposure with index {exposure_index}"
        )
    return matches[0]


def reconstruct_experimental_exposure(
    sequence: NonDestructiveSequence,
    exposure_index: int,
    *,
    measurement: DifferentiableDensityMeasurement,
    camera_calibrations: Mapping[str, CameraCalibration],
    model: NonnegativeBilinearDensityModel,
    observable_integration_support: ObservableIntegrationSupport,
    initial_coefficients: ArrayLike,
    coefficient_upper: float | ArrayLike,
    regularisation: CurvatureRegularisation | None,
    response_scale: ResponseScaleDeclaration,
    pipeline_fingerprint: str,
    condition_fingerprint: str,
    fit_options: DensityFitOptions | None = None,
    null_reference: NullReferenceDistribution | None = None,
    predeclared_alpha: float | None = None,
    bootstrap_request: BootstrapRequest | None = None,
    stability_variants: Mapping[str, DensityFitResult] | None = None,
    minimum_integrated_response: float = 0.0,
    retain_latent_artifacts: bool = False,
) -> ExperimentalReconstructionResult:
    """Fit one sequence exposure and return a provenance-bearing evidence bundle.

    The caller-supplied pipeline fingerprint must cover the basis and support,
    regulariser, initialisation/multi-start rule, fit options, observable
    contract and code version.  The exact integration support is also retained
    and hashed independently in the result.
    The condition fingerprint must cover the readout, ROI, pulse, detector and
    electron calibration, and background treatment.  They are compared with
    any blank reference but are not invented from incomplete metadata here.
    """

    if not isinstance(sequence, NonDestructiveSequence):
        raise TypeError("sequence must be a NonDestructiveSequence")
    if not isinstance(retain_latent_artifacts, bool):
        raise TypeError("retain_latent_artifacts must be boolean")
    if stability_variants is not None and not retain_latent_artifacts:
        raise ValueError(
            "legacy stability density maps require retain_latent_artifacts=True"
        )
    _validate_observable_integration_support(
        measurement,
        model,
        observable_integration_support,
    )
    pipeline = _nonempty_text(pipeline_fingerprint, "pipeline_fingerprint")
    condition = _nonempty_text(condition_fingerprint, "condition_fingerprint")
    exposure = _exposure_by_index(sequence, exposure_index)
    bound = bind_non_destructive_exposure(
        exposure,
        measurement,  # type: ignore[arg-type]
        camera_calibrations,
    )
    observed_channels = tuple(
        bound.channel_electrons[name] for name in bound.channel_electrons
    )
    fit = fit_nonnegative_basis_density(
        measurement,
        model,
        observed_channels,
        initial_coefficients=initial_coefficients,
        coefficient_upper=coefficient_upper,
        regularisation=regularisation,
        options=fit_options,
    )
    if not fit.diagnostics.success:
        raise RuntimeError(
            f"experimental density fit failed: {fit.diagnostics.message}"
        )

    response_scale.assert_matches_measurement(measurement)
    observables = extract_density_observables(
        fit.column_density_m2,
        observable_integration_support,
        minimum_integrated_response=minimum_integrated_response,
    )
    local = analyse_local_credibility(
        measurement,
        model,
        observed_channels,
        fit,
        coefficient_upper=coefficient_upper,
        regularisation=regularisation,
        retain_density_artifacts=retain_latent_artifacts,
    )
    null_evidence = compare_to_zero_density(
        measurement,
        observed_channels,
        fit,
        pipeline_fingerprint=pipeline,
        condition_fingerprint=condition,
        target_origin="experimental_observation",
        reference=null_reference,
        predeclared_alpha=predeclared_alpha,
    )

    if bootstrap_request is None:
        bootstrap = BootstrapAssessment(
            status="not_requested",
            result=None,
            message=None,
        )
    else:
        try:
            bootstrap_result = parametric_bootstrap_reconstruction(
                measurement,
                model,
                fit,
                coefficient_upper=coefficient_upper,
                regularisation=regularisation,
                draws=bootstrap_request.draws,
                seed=bootstrap_request.seed,
                confidence_level=bootstrap_request.confidence_level,
                options=bootstrap_request.fit_options or fit_options,
                observable_integration_support=observable_integration_support,
                minimum_integrated_response=minimum_integrated_response,
                retain_latent_artifacts=retain_latent_artifacts,
            )
        except RuntimeError as error:
            bootstrap = BootstrapAssessment(
                status="failed",
                result=None,
                message=str(error),
            )
        else:
            bootstrap_status: BootstrapStatus = (
                "complete"
                if bootstrap_result.successful_draws
                == bootstrap_result.requested_draws
                else "partial"
            )
            bootstrap = BootstrapAssessment(
                status=bootstrap_status,
                result=bootstrap_result,
                message=(
                    f"Conditional detector-noise bootstrap converged for "
                    f"{bootstrap_result.successful_draws}/"
                    f"{bootstrap_result.requested_draws} requested draws, with "
                    "fixed forward model, latent basis/support, observable "
                    "integration support and regulariser."
                ),
            )

    stability = (
        None
        if stability_variants is None
        else summarise_reconstruction_stability(stability_variants, measurement)
    )
    absolute_status: AbsoluteDensityStatus = (
        "calibrated"
        if response_scale.status == "calibrated"
        else "conditional_on_assumed_kappa_f"
    )
    limits = list(bound.warnings)
    limits.extend(
        [
            (
                "No truth density, TF state or eGPE state is supplied to this "
                "reconstruction."
            ),
            (
                "The null comparison is not an analytic likelihood-ratio "
                "significance test."
            ),
            (
                "Local uncertainty and any bootstrap are conditional on the "
                "declared forward and inverse contracts."
            ),
            (
                "This exposure is assessed individually; no sequence-level "
                "detection probability is inferred."
            ),
        ]
    )
    if absolute_status != "calibrated":
        limits.append(
            "The integrated response remains conditional on the assumed kappa_F; "
            "centroid and major-axis rms width are invariant to a uniform "
            "response-scale change but remain model-dependent."
        )
    if bootstrap_request is None:
        limits.append(
            "No conditional interval was requested for the reported physical "
            "observables."
        )
    if not retain_latent_artifacts:
        limits.append(
            "The latent density fit and density-uncertainty image were not retained; "
            "the formal output is the fixed-support observable estimate."
        )

    channel_names = tuple(bound.channel_electrons)
    return ExperimentalReconstructionResult(
        context_id=sequence.context.context_id,
        sequence_id=sequence.sequence_id,
        exposure_index=bound.exposure_index,
        timestamp_s=bound.timestamp_s,
        readout=bound.readout,
        channel_names=channel_names,
        input_channel_sha256=tuple(
            (name, _array_sha256(bound.channel_electrons[name]))
            for name in channel_names
        ),
        pipeline_fingerprint=pipeline,
        condition_fingerprint=condition,
        observable_support_sha256=_observable_support_sha256(
            observable_integration_support
        ),
        response_scale=response_scale,
        absolute_density_status=absolute_status,
        fit_diagnostics=fit.diagnostics,
        observables=observables,
        latent_fit=fit if retain_latent_artifacts else None,
        local_credibility=local,
        null_evidence=null_evidence,
        bootstrap=bootstrap,
        stability=stability,
        interpretation_limits=tuple(dict.fromkeys(limits)),
    )
