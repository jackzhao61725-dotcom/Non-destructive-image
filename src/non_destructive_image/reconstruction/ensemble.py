"""Noise-ensemble selection and held-out assessment for reconstruction."""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .benchmark import DensityRecoveryMetrics, density_recovery_metrics
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
class ReconstructionCandidate:
    """One declared basis and regularisation choice."""

    label: str
    model: NonnegativeBilinearDensityModel = field(repr=False)
    coefficient_upper: float | NDArray[np.floating]
    regularisation: CurvatureRegularisation | None = field(repr=False)
    fit_options: DensityFitOptions = field(default_factory=DensityFitOptions)

    def __post_init__(self) -> None:
        if not self.label.strip():
            raise ValueError("reconstruction candidate label cannot be empty")
        if (
            self.regularisation is not None
            and self.regularisation.parameter_count != self.model.parameter_count
        ):
            raise ValueError("candidate regularisation does not match the density model")
        upper = np.asarray(self.coefficient_upper, dtype=float)
        if upper.ndim == 0:
            if not np.isfinite(upper) or float(upper) <= 0:
                raise ValueError("candidate coefficient upper bound must be positive")
        elif upper.shape != (self.model.parameter_count,):
            raise ValueError(
                f"candidate coefficient upper bound must have shape ({self.model.parameter_count},)"
            )
        elif np.any(~np.isfinite(upper)) or np.any(upper <= 0):
            raise ValueError("candidate coefficient upper bounds must be finite and positive")


@dataclass(frozen=True)
class SyntheticNoisyObservation:
    """Raw noisy channels for one morphology and deterministic realization."""

    morphology: SyntheticMorphology
    realization_index: int
    seed: int
    channels: tuple[NDArray[np.floating], ...]


@dataclass(frozen=True)
class ReconstructionTrial:
    """Outcome of one candidate applied to one fixed noisy observation."""

    candidate_label: str
    morphology_name: str
    realization_index: int
    observation_seed: int
    success: bool
    message: str
    metrics: DensityRecoveryMetrics | None
    data_jacobian_rank: int
    data_jacobian_condition: float
    parameter_count: int


@dataclass(frozen=True)
class CandidateEnsembleSummary:
    """Compact selection quantities aggregated over calibration trials."""

    candidate_label: str
    morphology_names: tuple[str, ...]
    observation_keys: tuple[tuple[str, int, int], ...]
    trial_count: int
    success_fraction: float
    full_rank_fraction: float
    median_supported_band_error: float
    upper_quartile_supported_band_error: float
    median_absolute_integrated_density_error: float
    median_data_jacobian_condition: float
    maximum_data_jacobian_condition: float
    parameter_count: int


@dataclass(frozen=True)
class CandidateEnsembleAssessment:
    """All trials and their summary for one reconstruction candidate."""

    candidate: ReconstructionCandidate
    trials: tuple[ReconstructionTrial, ...]
    summary: CandidateEnsembleSummary


@dataclass(frozen=True)
class FrozenReconstructionChoice:
    """Candidate selected using only the named calibration morphologies."""

    candidate: ReconstructionCandidate
    calibration_morphology_names: tuple[str, ...]
    minimum_success_fraction: float
    relative_error_tolerance: float


@dataclass(frozen=True)
class HeldOutAssessment:
    """Frozen candidate evaluated on disjoint synthetic morphologies."""

    choice: FrozenReconstructionChoice
    assessment: CandidateEnsembleAssessment


CandidateInitialiser = Callable[
    [ReconstructionCandidate, tuple[NDArray[np.floating], ...]],
    DensityInitialisation,
]
SuccessfulFitCallback = Callable[
    [ReconstructionTrial, SyntheticNoisyObservation, DensityFitResult],
    None,
]


def make_linear_candidate_initialiser(
    measurement: DifferentiableDensityMeasurement,
    *,
    ridge_strength: float = 0.0,
) -> CandidateInitialiser:
    """Return a candidate-aware zero-density linear initializer."""

    def initialise(
        candidate: ReconstructionCandidate,
        channels: tuple[NDArray[np.floating], ...],
    ) -> DensityInitialisation:
        return linearised_nonnegative_initialisation(
            measurement,
            candidate.model,
            channels,
            coefficient_upper=candidate.coefficient_upper,
            ridge_strength=ridge_strength,
        )

    return initialise


def make_dark_field_candidate_initialiser(
    measurement: DarkFieldFaradayMeasurement,
    *,
    smooth_bounds: SmoothTFBounds,
    projection_ridge_strength: float = 0.0,
) -> CandidateInitialiser:
    """Return a candidate-aware square-root dark-field initializer."""

    def initialise(
        candidate: ReconstructionCandidate,
        channels: tuple[NDArray[np.floating], ...],
    ) -> DensityInitialisation:
        if len(channels) != 1:
            raise ValueError("dark-field candidate initialisation expects one count channel")
        return dark_field_sqrt_moment_initialisation(
            measurement,
            candidate.model,
            channels[0],
            smooth_bounds=smooth_bounds,
            coefficient_upper=candidate.coefficient_upper,
            projection_ridge_strength=projection_ridge_strength,
        )

    return initialise


def generate_noisy_observation_ensemble(
    measurement: DifferentiableDensityMeasurement,
    morphologies: Iterable[SyntheticMorphology],
    *,
    realizations_per_morphology: int,
    base_seed: int,
) -> tuple[SyntheticNoisyObservation, ...]:
    """Generate one reusable raw-count ensemble for fair candidate comparison."""

    if realizations_per_morphology <= 0:
        raise ValueError("realizations_per_morphology must be positive")
    morphology_list = list(morphologies)
    names = [morphology.name for morphology in morphology_list]
    if len(names) != len(set(names)):
        raise ValueError("morphology names must be unique within an ensemble")
    seed_sequence = np.random.SeedSequence(base_seed)
    child_sequences = seed_sequence.spawn(len(morphology_list) * realizations_per_morphology)
    observations: list[SyntheticNoisyObservation] = []
    sequence_index = 0
    for morphology in morphology_list:
        if morphology.column_density_m2.shape != measurement.grid.y_grid_m.shape:
            raise ValueError(
                f"morphology {morphology.name!r} does not match the measurement grid"
            )
        for realization_index in range(realizations_per_morphology):
            child = child_sequences[sequence_index]
            sequence_index += 1
            seed = int(child.generate_state(1, dtype=np.uint32)[0])
            rng = np.random.default_rng(seed)
            observations.append(
                SyntheticNoisyObservation(
                    morphology=morphology,
                    realization_index=realization_index,
                    seed=seed,
                    channels=measurement.simulate_channels_from_density(
                        morphology.column_density_m2,
                        rng,
                    ),
                )
            )
    return tuple(observations)


def _summary(
    candidate: ReconstructionCandidate,
    trials: tuple[ReconstructionTrial, ...],
) -> CandidateEnsembleSummary:
    if not trials:
        raise ValueError("candidate assessment requires at least one trial")
    successful = [trial for trial in trials if trial.success and trial.metrics is not None]
    supported_errors = np.asarray(
        [trial.metrics.supported_band_relative_l2_error for trial in successful],
        dtype=float,
    )
    integrated_errors = np.asarray(
        [abs(trial.metrics.integrated_density_relative_error) for trial in successful],
        dtype=float,
    )
    conditions = np.asarray(
        [trial.data_jacobian_condition for trial in successful],
        dtype=float,
    )
    return CandidateEnsembleSummary(
        candidate_label=candidate.label,
        morphology_names=tuple(sorted({trial.morphology_name for trial in trials})),
        observation_keys=tuple(
            sorted(
                (
                    trial.morphology_name,
                    trial.realization_index,
                    trial.observation_seed,
                )
                for trial in trials
            )
        ),
        trial_count=len(trials),
        success_fraction=len(successful) / len(trials),
        full_rank_fraction=sum(
            trial.data_jacobian_rank == trial.parameter_count for trial in trials
        )
        / len(trials),
        median_supported_band_error=(
            float(np.median(supported_errors)) if supported_errors.size else float("inf")
        ),
        upper_quartile_supported_band_error=(
            float(np.quantile(supported_errors, 0.75))
            if supported_errors.size
            else float("inf")
        ),
        median_absolute_integrated_density_error=(
            float(np.median(integrated_errors)) if integrated_errors.size else float("inf")
        ),
        median_data_jacobian_condition=(
            float(np.median(conditions)) if conditions.size else float("inf")
        ),
        maximum_data_jacobian_condition=(
            float(np.max(conditions)) if conditions.size else float("inf")
        ),
        parameter_count=candidate.model.parameter_count,
    )


def assess_reconstruction_candidate(
    measurement: DifferentiableDensityMeasurement,
    candidate: ReconstructionCandidate,
    observations: Iterable[SyntheticNoisyObservation],
    *,
    initialise: CandidateInitialiser,
    on_successful_fit: SuccessfulFitCallback | None = None,
) -> CandidateEnsembleAssessment:
    """Apply one candidate to a previously generated shared ensemble."""

    trials: list[ReconstructionTrial] = []
    for observation in observations:
        try:
            initialisation = initialise(candidate, observation.channels)
            fit = fit_nonnegative_basis_density(
                measurement,
                candidate.model,
                observation.channels,
                initial_coefficients=initialisation.coefficients,
                coefficient_upper=candidate.coefficient_upper,
                regularisation=candidate.regularisation,
                options=candidate.fit_options,
            )
            success = bool(fit.diagnostics.success)
            metrics = (
                density_recovery_metrics(
                    observation.morphology.column_density_m2,
                    fit.column_density_m2,
                    measurement.grid,
                )
                if success
                else None
            )
            message = fit.diagnostics.message
            rank = fit.diagnostics.data_jacobian_rank
            condition = fit.diagnostics.data_jacobian_condition
        except (RuntimeError, ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
            success = False
            metrics = None
            message = f"{type(error).__name__}: {error}"
            rank = 0
            condition = float("inf")
        trial = ReconstructionTrial(
            candidate_label=candidate.label,
            morphology_name=observation.morphology.name,
            realization_index=observation.realization_index,
            observation_seed=observation.seed,
            success=success,
            message=message,
            metrics=metrics,
            data_jacobian_rank=rank,
            data_jacobian_condition=condition,
            parameter_count=candidate.model.parameter_count,
        )
        trials.append(trial)
        if success and on_successful_fit is not None:
            on_successful_fit(trial, observation, fit)
    trial_tuple = tuple(trials)
    return CandidateEnsembleAssessment(
        candidate=candidate,
        trials=trial_tuple,
        summary=_summary(candidate, trial_tuple),
    )


def select_and_freeze_candidate(
    assessments: Iterable[CandidateEnsembleAssessment],
    *,
    minimum_success_fraction: float = 0.95,
    relative_error_tolerance: float = 0.02,
) -> FrozenReconstructionChoice:
    """Select on calibration data, preferring simpler near-equivalent models.

    Candidates must satisfy both the success and full-rank requirements.  The
    smallest median supported-band error defines the best value.  Candidates
    within ``relative_error_tolerance`` of it are treated as equivalent, and
    the one with the fewest coefficients is frozen.  The final tie-break is the
    candidate label, making selection deterministic.
    """

    if not 0 < minimum_success_fraction <= 1:
        raise ValueError("minimum success fraction must lie in (0, 1]")
    if not np.isfinite(relative_error_tolerance) or relative_error_tolerance < 0:
        raise ValueError("relative error tolerance must be finite and non-negative")
    assessment_list = list(assessments)
    if not assessment_list:
        raise ValueError("at least one candidate assessment is required")
    labels = [assessment.candidate.label for assessment in assessment_list]
    if len(labels) != len(set(labels)):
        raise ValueError("candidate labels must be unique")
    observation_sets = {assessment.summary.observation_keys for assessment in assessment_list}
    if len(observation_sets) != 1:
        raise ValueError("all candidates must be assessed on the same noisy observations")
    eligible = [
        assessment
        for assessment in assessment_list
        if assessment.summary.success_fraction >= minimum_success_fraction
        and assessment.summary.full_rank_fraction >= minimum_success_fraction
        and np.isfinite(assessment.summary.median_supported_band_error)
    ]
    if not eligible:
        raise RuntimeError("no candidate satisfies the calibration success and rank criteria")
    best_error = min(
        assessment.summary.median_supported_band_error for assessment in eligible
    )
    near_best = [
        assessment
        for assessment in eligible
        if assessment.summary.median_supported_band_error
        <= best_error * (1.0 + relative_error_tolerance)
    ]
    selected = min(
        near_best,
        key=lambda assessment: (
            assessment.summary.parameter_count,
            assessment.summary.median_supported_band_error,
            assessment.candidate.label,
        ),
    )
    return FrozenReconstructionChoice(
        candidate=selected.candidate,
        calibration_morphology_names=selected.summary.morphology_names,
        minimum_success_fraction=float(minimum_success_fraction),
        relative_error_tolerance=float(relative_error_tolerance),
    )


def evaluate_frozen_candidate_on_held_out(
    measurement: DifferentiableDensityMeasurement,
    choice: FrozenReconstructionChoice,
    observations: Iterable[SyntheticNoisyObservation],
    *,
    initialise: CandidateInitialiser,
    on_successful_fit: SuccessfulFitCallback | None = None,
) -> HeldOutAssessment:
    """Evaluate a frozen choice and reject calibration/held-out overlap."""

    observation_tuple = tuple(observations)
    held_out_names = {observation.morphology.name for observation in observation_tuple}
    overlap = held_out_names.intersection(choice.calibration_morphology_names)
    if overlap:
        raise ValueError(
            "held-out morphologies overlap the calibration set: " + ", ".join(sorted(overlap))
        )
    assessment = assess_reconstruction_candidate(
        measurement,
        choice.candidate,
        observation_tuple,
        initialise=initialise,
        on_successful_fit=on_successful_fit,
    )
    return HeldOutAssessment(choice=choice, assessment=assessment)
