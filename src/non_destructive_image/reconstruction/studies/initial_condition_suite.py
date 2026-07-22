"""DPFI sweep of source-traceable analytic initial conditions.

This study applies the already-frozen dual-port reconstruction candidate to
independent synthetic initial states.  It does not tune or reselect the
inverse.  Physical recovery is judged separately for the supported integral,
centroid and major-axis rms width; no overall image-usability decision is
defined.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
from time import perf_counter
from typing import Any, Callable

import numpy as np
from numpy.typing import NDArray

from ..density_fit import DensityFitResult
from ..ensemble import (
    ReconstructionCandidate,
    SyntheticNoisyObservation,
    assess_reconstruction_candidate,
    make_linear_candidate_initialiser,
)
from ..initial_conditions import InitialConditionTruth, build_initial_condition_truth
from ..observable_usability import (
    CENTROID_POSITION_ERROR_LIMIT_UM,
    INTEGRATED_RESPONSE_RELATIVE_ERROR_LIMIT,
    MAJOR_RMS_WIDTH_RELATIVE_ERROR_LIMIT,
    ObservableUsabilityEvaluation,
    evaluate_observable_usability,
)
from ..observables import (
    DensityObservableSummary,
    ObservableIntegrationSupport,
    extract_density_observables,
)
from .io import file_sha256, load_json, write_json, write_rows
from .morphology import (
    MorphologyStudyContext,
    build_morphology_study_context,
    make_study_measurement,
)
from .observable_benchmark import (
    SourceBenchmarkArtifacts,
    build_observable_integration_support,
    validate_source_benchmark_artifacts,
)
from .provenance import capture_reconstruction_provenance


ProgressCallback = Callable[[str], None]

TRIALS_FILENAME = "condition_trials.csv"
CONDITION_SUMMARY_FILENAME = "condition_summary.csv"
MAPS_FILENAME = "condition_maps.npz"
STUDY_CONFIG_FILENAME = "study_config.json"
SOURCE_CONFIG_FILENAME = "source_benchmark_config.json"
SOURCE_MANIFEST_FILENAME = "source_benchmark_manifest.json"
PROVENANCE_FILENAME = "generation_provenance.json"
SUMMARY_FILENAME = "suite_summary.json"
MANIFEST_FILENAME = "artifact_manifest.json"
MANIFEST_SCHEMA_VERSION = 1

RAW_ARTIFACT_ROLES: dict[str, str] = {
    STUDY_CONFIG_FILENAME: "configuration_snapshot",
    SOURCE_CONFIG_FILENAME: "source_configuration_snapshot",
    SOURCE_MANIFEST_FILENAME: "source_manifest_snapshot",
    PROVENANCE_FILENAME: "generation_provenance",
    TRIALS_FILENAME: "condition_trials",
    CONDITION_SUMMARY_FILENAME: "condition_summary",
    MAPS_FILENAME: "truth_reconstruction_raw_channels_and_coefficients",
}


@dataclass(frozen=True)
class InitialConditionSuiteContext:
    """Validated source lineage, grids, candidate and analytic truth maps."""

    config: dict[str, Any] = field(repr=False)
    source: SourceBenchmarkArtifacts = field(repr=False)
    reconstruction: MorphologyStudyContext = field(repr=False)
    integration_support: ObservableIntegrationSupport = field(repr=False)
    candidate: ReconstructionCandidate = field(repr=False)
    truths: tuple[InitialConditionTruth, ...] = field(repr=False)
    fluences_mw_us: tuple[float, ...]
    realizations_per_condition: int


@dataclass(frozen=True)
class InitialConditionTrial:
    """One fixed noisy DPFI observation and frozen-candidate fit."""

    condition_index: int
    fluence_index: int
    realization_index: int
    seed: int
    candidate_label: str
    truth: InitialConditionTruth = field(repr=False)
    channels_e: tuple[NDArray[np.floating], NDArray[np.floating]] = field(
        repr=False
    )
    fit_success: bool
    fit_message: str
    data_jacobian_rank: int
    data_jacobian_condition: float
    reconstructed_column_density_m2: NDArray[np.floating] | None = field(
        repr=False
    )
    fitted_coefficients: NDArray[np.floating] | None = field(repr=False)
    truth_observables: DensityObservableSummary = field(repr=False)
    reconstructed_observables: DensityObservableSummary | None = field(repr=False)
    usability: ObservableUsabilityEvaluation | None

    @property
    def condition_id(self) -> str:
        return self.truth.morphology.name


@dataclass(frozen=True)
class InitialConditionSuiteRun:
    """Complete in-memory DPFI suite before artifact serialization."""

    context: InitialConditionSuiteContext = field(repr=False)
    trials: tuple[InitialConditionTrial, ...] = field(repr=False)
    elapsed_seconds: float


def _progress(callback: ProgressCallback | None, message: str) -> None:
    if callback is not None:
        callback(message)


def _configured_output_directory(
    config: dict[str, Any], repository_root: Path
) -> Path:
    configured = Path(str(config["output_directory"]))
    output = (
        configured.resolve()
        if configured.is_absolute()
        else (repository_root / configured).resolve()
    )
    try:
        output.relative_to(repository_root.resolve())
    except ValueError as error:
        raise ValueError("suite output directory must remain inside the repository") from error
    return output


def _selected_dual_port_candidate(
    reconstruction: MorphologyStudyContext,
    source: SourceBenchmarkArtifacts,
) -> ReconstructionCandidate:
    selected = source.summary["selected_candidates"]["dual_port"]
    label = str(selected["label"])
    matches = [candidate for candidate in reconstruction.candidates if candidate.label == label]
    if len(matches) != 1:
        raise ValueError(f"expected one frozen DPFI candidate labelled {label!r}")
    candidate = matches[0]
    if candidate.model.parameter_count != int(selected["parameter_count"]):
        raise ValueError("frozen DPFI candidate parameter count disagrees with source")
    return candidate


def _validate_usability_contract(config: dict[str, Any]) -> None:
    configured = config.get("observable_usability")
    if not isinstance(configured, dict):
        raise ValueError("suite config must declare observable_usability")
    expected = {
        "integrated_response_relative_error_tolerance": (
            INTEGRATED_RESPONSE_RELATIVE_ERROR_LIMIT
        ),
        "centroid_position_error_um_tolerance": CENTROID_POSITION_ERROR_LIMIT_UM,
        "major_rms_width_relative_error_tolerance": (
            MAJOR_RMS_WIDTH_RELATIVE_ERROR_LIMIT
        ),
    }
    for key, value in expected.items():
        if float(configured.get(key, float("nan"))) != value:
            raise ValueError(
                f"configured {key} must match the frozen approved value {value:g}"
            )


def _fluence_seed_key(fluence_mw_us: float) -> str:
    return f"{fluence_mw_us:g}"


def derive_trial_seed(
    base_seed: int,
    condition_id: str,
    realization_index: int,
) -> int:
    """Derive a stable uint32 seed independent of condition ordering."""

    if not condition_id:
        raise ValueError("condition_id must be non-empty")
    if realization_index < 0:
        raise ValueError("realization_index must be non-negative")
    payload = json.dumps(
        {
            "base_seed": int(base_seed),
            "condition_id": condition_id,
            "realization_index": int(realization_index),
            "schema": "dpfi_initial_condition_suite_seed_v1",
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:4], "little")


def validate_initial_condition_suite(
    config: dict[str, Any],
    repository_root: Path,
) -> InitialConditionSuiteContext:
    """Validate source lineage and build every truth map without fitting."""

    repository_root = repository_root.resolve()
    _validate_usability_contract(config)
    source = validate_source_benchmark_artifacts(config, repository_root)
    reconstruction = build_morphology_study_context(source.config)
    integration_support = build_observable_integration_support(
        config,
        source,
        reconstruction,
    )
    candidate = _selected_dual_port_candidate(reconstruction, source)

    raw_conditions = config.get("initial_conditions")
    if not isinstance(raw_conditions, list) or not raw_conditions:
        raise ValueError("suite config must contain initial_conditions")
    identifiers = [str(condition.get("id", "")) for condition in raw_conditions]
    if any(not identifier for identifier in identifiers):
        raise ValueError("every initial condition must have a non-empty id")
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("initial-condition ids must be unique")
    atomic_constants = config.get("atomic_constants")
    if not isinstance(atomic_constants, dict):
        raise ValueError("suite config must declare atomic_constants")
    truths = tuple(
        build_initial_condition_truth(
            reconstruction.grid.y_grid_m,
            reconstruction.grid.z_grid_m,
            condition,
            atomic_constants,
        )
        for condition in raw_conditions
    )
    if tuple(truth.morphology.name for truth in truths) != tuple(identifiers):
        raise ValueError("generated truth names do not match configured condition ids")
    for truth in truths:
        if truth.morphology.column_density_m2.shape != integration_support.shape:
            raise ValueError(
                f"truth map {truth.morphology.name!r} does not match reconstruction grid"
            )
        summary = extract_density_observables(
            truth.morphology.column_density_m2,
            integration_support,
        )
        if summary.integrated_response <= 0.0:
            raise ValueError(
                f"truth map {truth.morphology.name!r} has zero response on fixed support"
            )

    ensemble = config.get("ensemble")
    if not isinstance(ensemble, dict):
        raise ValueError("suite config must declare ensemble")
    fluences = tuple(float(value) for value in ensemble.get("fluence_mw_us", ()))
    if not fluences or len(fluences) != len(set(fluences)):
        raise ValueError("suite fluences must be non-empty and unique")
    if any(not np.isfinite(value) or value <= 0.0 for value in fluences):
        raise ValueError("suite fluences must be finite and positive")
    realization_count = int(ensemble.get("realizations_per_condition", 0))
    if realization_count != 1:
        raise ValueError("version-1 suite requires exactly one realization per condition")
    seeds = ensemble.get("base_seed_by_fluence")
    if not isinstance(seeds, dict):
        raise ValueError("suite config must declare base_seed_by_fluence")
    expected_seed_keys = {_fluence_seed_key(value) for value in fluences}
    if set(seeds) != expected_seed_keys:
        raise ValueError("base_seed_by_fluence keys must match the configured fluences")

    output = _configured_output_directory(config, repository_root)
    if output == source.directory or source.directory in output.parents:
        raise ValueError("suite outputs must not be written into the sealed source family")
    return InitialConditionSuiteContext(
        config=config,
        source=source,
        reconstruction=reconstruction,
        integration_support=integration_support,
        candidate=candidate,
        truths=truths,
        fluences_mw_us=fluences,
        realizations_per_condition=realization_count,
    )


def _condition_semantics(truth: InitialConditionTruth) -> str:
    generator = truth.generator_name
    if generator in {"bimodal_ideal_gas", "tf_roton_modulation"}:
        return "support_restricted_total_column_density"
    return "support_restricted_column_density"


def _simulate_observations(
    context: InitialConditionSuiteContext,
    fluence_index: int,
) -> tuple[SyntheticNoisyObservation, ...]:
    fluence = context.fluences_mw_us[fluence_index]
    measurement = make_study_measurement(
        context.reconstruction,
        "dual_port",
        fluence_mw_us=fluence,
    )
    base_seed = int(
        context.config["ensemble"]["base_seed_by_fluence"][
            _fluence_seed_key(fluence)
        ]
    )
    observations: list[SyntheticNoisyObservation] = []
    for truth in context.truths:
        for realization_index in range(context.realizations_per_condition):
            seed = derive_trial_seed(
                base_seed,
                truth.morphology.name,
                realization_index,
            )
            rng = np.random.default_rng(seed)
            channels = measurement.simulate_channels_from_density(
                truth.morphology.column_density_m2,
                rng,
            )
            if len(channels) != 2:
                raise RuntimeError("DPFI simulation did not produce two raw channels")
            observations.append(
                SyntheticNoisyObservation(
                    morphology=truth.morphology,
                    realization_index=realization_index,
                    seed=seed,
                    channels=channels,
                )
            )
    seeds = [observation.seed for observation in observations]
    if len(seeds) != len(set(seeds)):
        raise RuntimeError("derived DPFI trial seeds are not unique within a fluence")
    return tuple(observations)


def run_initial_condition_suite(
    context: InitialConditionSuiteContext,
    *,
    progress: ProgressCallback | None = None,
) -> InitialConditionSuiteRun:
    """Run the frozen DPFI candidate on all configured initial conditions."""

    started = perf_counter()
    truth_by_id = {truth.morphology.name: truth for truth in context.truths}
    truth_observables = {
        name: extract_density_observables(
            truth.morphology.column_density_m2,
            context.integration_support,
        )
        for name, truth in truth_by_id.items()
    }
    condition_index = {
        truth.morphology.name: index for index, truth in enumerate(context.truths)
    }
    trials: list[InitialConditionTrial] = []
    for fluence_index, fluence in enumerate(context.fluences_mw_us):
        _progress(
            progress,
            f"DPFI: fitting {len(context.truths)} conditions at F={fluence:g} mW us",
        )
        measurement = make_study_measurement(
            context.reconstruction,
            "dual_port",
            fluence_mw_us=fluence,
        )
        observations = _simulate_observations(context, fluence_index)
        captures: dict[tuple[str, int, int], DensityFitResult] = {}

        def capture(
            trial: Any,
            observation: SyntheticNoisyObservation,
            fit: DensityFitResult,
        ) -> None:
            key = (
                observation.morphology.name,
                observation.realization_index,
                observation.seed,
            )
            if key in captures:
                raise RuntimeError(f"duplicate successful fit capture: {key}")
            captures[key] = fit
            _progress(progress, f"  {observation.morphology.name}: success")

        assessment = assess_reconstruction_candidate(
            measurement,
            context.candidate,
            observations,
            initialise=make_linear_candidate_initialiser(
                measurement,
                ridge_strength=float(
                    context.reconstruction.config["fit"][
                        "dual_port_initialisation_ridge_strength"
                    ]
                ),
            ),
            on_successful_fit=capture,
        )
        observation_by_key = {
            (
                observation.morphology.name,
                observation.realization_index,
                observation.seed,
            ): observation
            for observation in observations
        }
        for reconstruction_trial in assessment.trials:
            key = (
                reconstruction_trial.morphology_name,
                reconstruction_trial.realization_index,
                reconstruction_trial.observation_seed,
            )
            observation = observation_by_key[key]
            truth = truth_by_id[reconstruction_trial.morphology_name]
            fit = captures.get(key)
            reconstructed_summary: DensityObservableSummary | None = None
            usability: ObservableUsabilityEvaluation | None = None
            reconstructed_map: NDArray[np.floating] | None = None
            coefficients: NDArray[np.floating] | None = None
            if reconstruction_trial.success:
                if fit is None:
                    raise RuntimeError(f"successful fit was not captured: {key}")
                reconstructed_map = np.array(
                    fit.column_density_m2,
                    dtype=float,
                    copy=True,
                )
                coefficients = np.array(fit.coefficients, dtype=float, copy=True)
                reconstructed_summary = extract_density_observables(
                    reconstructed_map,
                    context.integration_support,
                )
                usability = evaluate_observable_usability(
                    truth_observables[reconstruction_trial.morphology_name],
                    reconstructed_summary,
                )
            trials.append(
                InitialConditionTrial(
                    condition_index=condition_index[reconstruction_trial.morphology_name],
                    fluence_index=fluence_index,
                    realization_index=reconstruction_trial.realization_index,
                    seed=reconstruction_trial.observation_seed,
                    candidate_label=context.candidate.label,
                    truth=truth,
                    channels_e=(
                        np.array(observation.channels[0], dtype=float, copy=True),
                        np.array(observation.channels[1], dtype=float, copy=True),
                    ),
                    fit_success=reconstruction_trial.success,
                    fit_message=reconstruction_trial.message,
                    data_jacobian_rank=reconstruction_trial.data_jacobian_rank,
                    data_jacobian_condition=reconstruction_trial.data_jacobian_condition,
                    reconstructed_column_density_m2=reconstructed_map,
                    fitted_coefficients=coefficients,
                    truth_observables=truth_observables[
                        reconstruction_trial.morphology_name
                    ],
                    reconstructed_observables=reconstructed_summary,
                    usability=usability,
                )
            )
            if not reconstruction_trial.success:
                _progress(
                    progress,
                    f"  {reconstruction_trial.morphology_name}: failed",
                )
    expected = (
        len(context.truths)
        * len(context.fluences_mw_us)
        * context.realizations_per_condition
    )
    if len(trials) != expected:
        raise RuntimeError(f"suite produced {len(trials)} trials, expected {expected}")
    keys = {
        (trial.fluence_index, trial.condition_index, trial.realization_index)
        for trial in trials
    }
    if len(keys) != expected:
        raise RuntimeError("suite produced duplicate trial coordinates")
    return InitialConditionSuiteRun(
        context=context,
        trials=tuple(trials),
        elapsed_seconds=perf_counter() - started,
    )


def _optional(value: float | None, *, scale: float = 1.0) -> float | str:
    return "" if value is None else float(value) * scale


def _relative_error(recovered: float | None, truth: float | None) -> float | str:
    if recovered is None or truth is None or truth == 0.0:
        return ""
    return float(recovered / truth - 1.0)


def trial_row(trial: InitialConditionTrial, fluence_mw_us: float) -> dict[str, Any]:
    """Serialize one trial using only the three declared physical targets."""

    truth = trial.truth_observables
    reconstructed = trial.reconstructed_observables
    usability = trial.usability
    source_metadata = trial.truth.source_metadata
    reconstructed_integral = (
        None if reconstructed is None else reconstructed.integrated_response
    )
    reconstructed_centroid_y = (
        None if reconstructed is None else reconstructed.centroid_y_m
    )
    reconstructed_centroid_z = (
        None if reconstructed is None else reconstructed.centroid_z_m
    )
    reconstructed_width = (
        None if reconstructed is None else reconstructed.major_rms_width_m
    )
    centroid_error = ""
    if truth.centroid_m is not None and reconstructed is not None and reconstructed.centroid_m is not None:
        centroid_error = float(
            np.linalg.norm(reconstructed.centroid_m - truth.centroid_m) * 1e6
        )
    return {
        "condition_id": trial.condition_id,
        "family": str(source_metadata["family"]),
        "feature_class": trial.truth.feature_class,
        "generator": trial.truth.generator_name,
        "observable_semantics": _condition_semantics(trial.truth),
        "fluence_mw_us": fluence_mw_us,
        "realization_index": trial.realization_index,
        "seed": trial.seed,
        "candidate": trial.candidate_label,
        "fit_success": trial.fit_success,
        "fit_message": trial.fit_message,
        "data_jacobian_rank": trial.data_jacobian_rank,
        "data_jacobian_condition": trial.data_jacobian_condition,
        "truth_integrated_response": truth.integrated_response,
        "reconstructed_integrated_response": _optional(reconstructed_integral),
        "integrated_response_relative_error": _relative_error(
            reconstructed_integral,
            truth.integrated_response,
        ),
        "c_A": "" if usability is None or usability.c_a is None else usability.c_a,
        "integrated_response_passed": (
            False if usability is None else usability.integrated_response_passed
        ),
        "truth_centroid_y_um": _optional(truth.centroid_y_m, scale=1e6),
        "truth_centroid_z_um": _optional(truth.centroid_z_m, scale=1e6),
        "reconstructed_centroid_y_um": _optional(
            reconstructed_centroid_y,
            scale=1e6,
        ),
        "reconstructed_centroid_z_um": _optional(
            reconstructed_centroid_z,
            scale=1e6,
        ),
        "centroid_position_error_um": centroid_error,
        "c_r": "" if usability is None or usability.c_r is None else usability.c_r,
        "centroid_position_passed": (
            False if usability is None else usability.centroid_position_passed
        ),
        "truth_major_rms_width_um": _optional(
            truth.major_rms_width_m,
            scale=1e6,
        ),
        "reconstructed_major_rms_width_um": _optional(
            reconstructed_width,
            scale=1e6,
        ),
        "major_rms_width_relative_error": _relative_error(
            reconstructed_width,
            truth.major_rms_width_m,
        ),
        "c_w": "" if usability is None or usability.c_w is None else usability.c_w,
        "major_rms_width_passed": (
            False if usability is None else usability.major_rms_width_passed
        ),
    }


def aggregate_condition_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Return per-condition worst coefficients and independent pass decisions."""

    identifiers = list(dict.fromkeys(str(row["condition_id"]) for row in rows))
    output: list[dict[str, Any]] = []
    specifications = (
        ("c_A", "integrated_response_passed"),
        ("c_r", "centroid_position_passed"),
        ("c_w", "major_rms_width_passed"),
    )
    for identifier in identifiers:
        group = [row for row in rows if row["condition_id"] == identifier]
        first = group[0]
        aggregate: dict[str, Any] = {
            "condition_id": identifier,
            "family": first["family"],
            "feature_class": first["feature_class"],
            "generator": first["generator"],
            "observable_semantics": first["observable_semantics"],
            "trial_count": len(group),
            "successful_fit_count": sum(bool(row["fit_success"]) for row in group),
        }
        for score_key, pass_key in specifications:
            scores = [
                float(row[score_key]) for row in group if row.get(score_key, "") != ""
            ]
            aggregate[f"supported_{score_key}_trials"] = len(scores)
            aggregate[f"maximum_{score_key}"] = max(scores) if scores else ""
            aggregate[f"{pass_key}_all_trials"] = bool(
                len(scores) == len(group) and all(bool(row[pass_key]) for row in group)
            )
        output.append(aggregate)
    return output


def _write_maps(path: Path, run: InitialConditionSuiteRun) -> None:
    context = run.context
    fluences = np.asarray(context.fluences_mw_us, dtype=float)
    condition_ids = np.asarray(
        [truth.morphology.name for truth in context.truths]
    )
    realizations = np.arange(context.realizations_per_condition, dtype=int)
    leading_shape = (fluences.size, condition_ids.size, realizations.size)
    grid_shape = context.integration_support.shape
    camera_shape = run.trials[0].channels_e[0].shape
    parameter_count = context.candidate.model.parameter_count
    reconstructed = np.full((*leading_shape, *grid_shape), np.nan, dtype=float)
    coefficients = np.full((*leading_shape, parameter_count), np.nan, dtype=float)
    h_counts = np.full((*leading_shape, *camera_shape), np.nan, dtype=float)
    v_counts = np.full((*leading_shape, *camera_shape), np.nan, dtype=float)
    seeds = np.full(leading_shape, -1, dtype=np.int64)
    fit_success = np.zeros(leading_shape, dtype=bool)
    filled = np.zeros(leading_shape, dtype=bool)
    for trial in run.trials:
        index = (trial.fluence_index, trial.condition_index, trial.realization_index)
        if filled[index]:
            raise ValueError(f"duplicate map tensor coordinate: {index}")
        h_counts[index] = trial.channels_e[0]
        v_counts[index] = trial.channels_e[1]
        seeds[index] = trial.seed
        fit_success[index] = trial.fit_success
        if trial.reconstructed_column_density_m2 is not None:
            reconstructed[index] = trial.reconstructed_column_density_m2
        if trial.fitted_coefficients is not None:
            coefficients[index] = trial.fitted_coefficients
        filled[index] = True
    if not np.all(filled):
        raise ValueError("map tensor contains unfilled coordinates")
    np.savez_compressed(
        path,
        y_grid_m=context.integration_support.y_grid_m,
        z_grid_m=context.integration_support.z_grid_m,
        cell_area_m2=context.integration_support.cell_area_m2,
        support_mask=context.integration_support.support_mask,
        fluence_mw_us=fluences,
        condition_ids=condition_ids,
        realization_indices=realizations,
        candidate_label=np.asarray(context.candidate.label),
        seeds=seeds,
        fit_success=fit_success,
        truth_column_density_m2=np.stack(
            [truth.morphology.column_density_m2 for truth in context.truths]
        ),
        reconstructed_column_density_m2=reconstructed,
        fitted_coefficients=coefficients,
        dual_port_h_counts_e=h_counts,
        dual_port_v_counts_e=v_counts,
        source_run_id=np.asarray(context.source.run_id),
        tensor_axis_order=np.asarray(
            "fluence,condition,realization,z_index,y_index"
        ),
    )


def _raw_hashes(directory: Path) -> dict[str, str]:
    return {
        filename: file_sha256(directory / filename)
        for filename in RAW_ARTIFACT_ROLES
    }


def _deterministic_run_id(
    config_sha256: str,
    source_run_id: str,
    raw_hashes: dict[str, str],
) -> str:
    payload = json.dumps(
        {
            "config_sha256": config_sha256,
            "source_run_id": source_run_id,
            "raw_artifact_sha256": raw_hashes,
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _summary_payload(
    run: InitialConditionSuiteRun,
    rows: list[dict[str, Any]],
) -> dict[str, Any]:
    context = run.context
    return {
        "label": context.config["label"],
        "purpose": context.config["purpose"],
        "source_benchmark": {
            "label": context.source.manifest["label"],
            "run_id": context.source.run_id,
            "directory": context.config["source_benchmark"]["directory"],
            "manifest_sha256": file_sha256(context.source.manifest_path),
            "config_sha256": file_sha256(context.source.config_path),
        },
        "suite": {
            "readout": "dual_port",
            "trial_count": len(rows),
            "condition_ids": [truth.morphology.name for truth in context.truths],
            "fluence_mw_us": list(context.fluences_mw_us),
            "realization_indices": list(
                range(context.realizations_per_condition)
            ),
            "successful_fit_count": sum(bool(row["fit_success"]) for row in rows),
            "selected_candidate": context.candidate.label,
            "parameter_count": context.candidate.model.parameter_count,
            "elapsed_seconds": run.elapsed_seconds,
        },
        "observable_usability": {
            "scores_are_reported_independently": True,
            "overall_image_usability_score": None,
            "integrated_response_relative_error_tolerance": (
                INTEGRATED_RESPONSE_RELATIVE_ERROR_LIMIT
            ),
            "centroid_position_error_um_tolerance": (
                CENTROID_POSITION_ERROR_LIMIT_UM
            ),
            "major_rms_width_relative_error_tolerance": (
                MAJOR_RMS_WIDTH_RELATIVE_ERROR_LIMIT
            ),
            "truth_contract": context.config["observable_usability"][
                "truth_contract"
            ],
            "single_realization_status": context.config["ensemble"]["status"],
        },
        "integration_support": {
            "kind": context.config["integration_support"]["kind"],
            "grid_shape": list(context.integration_support.shape),
            "supported_cell_count": int(
                np.count_nonzero(context.integration_support.support_mask)
            ),
            "physical_area_m2": context.integration_support.physical_area_m2,
            "bimodal_and_thermal_semantics": (
                "support-restricted total column-density observables, not N0, "
                "Ntot, temperature or global thermal width"
            ),
        },
        "maps_schema": {
            "artifact": MAPS_FILENAME,
            "truth_shape": [len(context.truths), *context.integration_support.shape],
            "reconstruction_axis_order": [
                "fluence",
                "condition",
                "realization",
                "z_index",
                "y_index",
            ],
        },
        "tables": {
            "trials": TRIALS_FILENAME,
            "condition_summary": CONDITION_SUMMARY_FILENAME,
        },
        "claims_boundary": context.config["claims_boundary"],
    }


def write_initial_condition_suite(
    run: InitialConditionSuiteRun,
    config_path: Path,
    repository_root: Path,
    *,
    generation_provenance: dict[str, Any],
) -> dict[str, Path]:
    """Write one new sealed result family through a scratch staging directory."""

    repository_root = repository_root.resolve()
    output_directory = _configured_output_directory(run.context.config, repository_root)
    if output_directory.exists():
        raise FileExistsError(
            f"refusing to overwrite existing suite output: {output_directory}"
        )
    scratch_root = repository_root / ".scratch"
    scratch_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="dpfi_initial_condition_suite_",
        dir=scratch_root,
    ) as temporary:
        staging = Path(temporary) / output_directory.name
        staging.mkdir()
        shutil.copyfile(config_path, staging / STUDY_CONFIG_FILENAME)
        shutil.copyfile(run.context.source.config_path, staging / SOURCE_CONFIG_FILENAME)
        shutil.copyfile(
            run.context.source.manifest_path,
            staging / SOURCE_MANIFEST_FILENAME,
        )
        write_json(staging / PROVENANCE_FILENAME, generation_provenance)
        rows = [
            trial_row(trial, run.context.fluences_mw_us[trial.fluence_index])
            for trial in run.trials
        ]
        write_rows(staging / TRIALS_FILENAME, rows)
        write_rows(staging / CONDITION_SUMMARY_FILENAME, aggregate_condition_rows(rows))
        _write_maps(staging / MAPS_FILENAME, run)
        raw_hashes = _raw_hashes(staging)
        config_sha256 = file_sha256(staging / STUDY_CONFIG_FILENAME)
        run_id = _deterministic_run_id(
            config_sha256,
            run.context.source.run_id,
            raw_hashes,
        )
        summary = _summary_payload(run, rows)
        summary["run_id"] = run_id
        write_json(staging / SUMMARY_FILENAME, summary)
        artifacts = {
            filename: {"role": RAW_ARTIFACT_ROLES[filename], "sha256": digest}
            for filename, digest in raw_hashes.items()
        }
        artifacts[SUMMARY_FILENAME] = {
            "role": "run_summary",
            "sha256": file_sha256(staging / SUMMARY_FILENAME),
        }
        write_json(
            staging / MANIFEST_FILENAME,
            {
                "schema_version": MANIFEST_SCHEMA_VERSION,
                "label": run.context.config["label"],
                "run_id": run_id,
                "source_run_id": run.context.source.run_id,
                "config": {
                    "artifact": STUDY_CONFIG_FILENAME,
                    "sha256": config_sha256,
                },
                "artifacts": artifacts,
            },
        )
        output_directory.parent.mkdir(parents=True, exist_ok=True)
        staging.replace(output_directory)
    return {
        "trials": output_directory / TRIALS_FILENAME,
        "condition_summary": output_directory / CONDITION_SUMMARY_FILENAME,
        "maps": output_directory / MAPS_FILENAME,
        "summary": output_directory / SUMMARY_FILENAME,
        "manifest": output_directory / MANIFEST_FILENAME,
        "provenance": output_directory / PROVENANCE_FILENAME,
    }


def run_and_write_initial_condition_suite(
    config: dict[str, Any],
    config_path: Path,
    repository_root: Path,
    *,
    progress: ProgressCallback | None = None,
) -> dict[str, Path]:
    """Validate, run and seal the DPFI initial-condition sweep."""

    repository_root = repository_root.resolve()
    config_path = config_path.resolve()
    config_sha256 = file_sha256(config_path)
    context = validate_initial_condition_suite(config, repository_root)
    entry_points = (Path("scripts/run_dpfi_initial_condition_suite.py"),)
    before = capture_reconstruction_provenance(
        repository_root,
        entry_points=entry_points,
    )
    run = run_initial_condition_suite(context, progress=progress)
    if file_sha256(config_path) != config_sha256:
        raise RuntimeError("suite config changed during numerical generation")
    after = capture_reconstruction_provenance(
        repository_root,
        entry_points=entry_points,
    )
    for field_name in ("source_files_sha256", "runtime_versions"):
        if before[field_name] != after[field_name]:
            raise RuntimeError(f"suite {field_name} changed during numerical generation")
    provenance = {
        "study_config_source": str(config_path),
        "study_config_source_sha256": config_sha256,
        "source_benchmark": {
            "run_id": context.source.run_id,
            "manifest": str(context.source.manifest_path),
            "manifest_sha256": file_sha256(context.source.manifest_path),
            "config_sha256": file_sha256(context.source.config_path),
        },
        "generation_state_before_run": before,
        "generation_state_after_run": after,
        "trial_seed_contract": (
            "SHA-256(base seed by fluence, condition id, realization index); "
            "stable under condition reordering"
        ),
        "observable_contract": {
            "same_fixed_support_truth_and_reconstruction": True,
            "truth_is_not_pupil_filtered": True,
            "scores_are_independent": True,
            "overall_image_usability_score": None,
        },
    }
    return write_initial_condition_suite(
        run,
        config_path,
        repository_root,
        generation_provenance=provenance,
    )


__all__ = [
    "CONDITION_SUMMARY_FILENAME",
    "InitialConditionSuiteContext",
    "InitialConditionSuiteRun",
    "InitialConditionTrial",
    "MANIFEST_FILENAME",
    "MAPS_FILENAME",
    "SUMMARY_FILENAME",
    "TRIALS_FILENAME",
    "aggregate_condition_rows",
    "derive_trial_seed",
    "run_and_write_initial_condition_suite",
    "run_initial_condition_suite",
    "trial_row",
    "validate_initial_condition_suite",
    "write_initial_condition_suite",
]
