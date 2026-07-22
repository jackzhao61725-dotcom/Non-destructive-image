"""Held-out replay and sealed artifacts for physical density observables.

The source morphology benchmark intentionally retained only one representative
map.  This module validates that sealed benchmark, replays only its already
selected 60 held-out fits, and writes a separate result family.  It never
modifies or reseals the source benchmark.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
import shutil
from time import perf_counter
from typing import Any, Literal

import numpy as np
from numpy.typing import NDArray

from ..benchmark import DensityRecoveryMetrics, density_recovery_metrics
from ..density_fit import DensityFitResult, fit_nonnegative_basis_density
from ..density_initialise import (
    dark_field_sqrt_moment_initialisation,
    linearised_nonnegative_initialisation,
)
from ..ensemble import ReconstructionCandidate, generate_noisy_observation_ensemble
from ..measurements import DarkFieldFaradayMeasurement, DifferentiableDensityMeasurement
from ..observables import (
    DensityObservableSummary,
    ObservableIntegrationSupport,
    extract_density_observables,
)
from .io import file_sha256, load_json, load_rows, write_json, write_rows
from .morphology import (
    MorphologyStudyContext,
    build_morphology_study_context,
    make_study_measurement,
)
from .provenance import capture_reconstruction_provenance


ReadoutName = Literal["dual_port", "dark_field"]
READOUT_NAMES: tuple[ReadoutName, ...] = ("dual_port", "dark_field")
SOURCE_MANIFEST_SCHEMA_VERSION = 1
OBSERVABLE_MANIFEST_SCHEMA_VERSION = 1
OBSERVABLE_MANIFEST_FILENAME = "artifact_manifest.json"
OBSERVABLE_SUMMARY_FILENAME = "observable_benchmark_summary.json"
OBSERVABLE_MAPS_FILENAME = "held_out_observable_maps.npz"
OBSERVABLE_TRIALS_FILENAME = "observable_trials.csv"
OBSERVABLE_AGGREGATE_FILENAME = "observable_recovery_summary.csv"
LEGACY_VERIFICATION_FILENAME = "legacy_metric_verification.csv"
STUDY_CONFIG_SNAPSHOT_FILENAME = "study_config.json"
SOURCE_CONFIG_SNAPSHOT_FILENAME = "source_benchmark_config.json"
SOURCE_MANIFEST_SNAPSHOT_FILENAME = "source_benchmark_manifest.json"
GENERATION_PROVENANCE_FILENAME = "generation_provenance.json"

_SOURCE_RAW_ROLES = frozenset(
    {
        "calibration_candidates",
        "held_out_trials",
        "held_out_summary",
        "grid_convergence",
        "representative_arrays",
    }
)
OBSERVABLE_RAW_ARTIFACT_ROLES: dict[str, str] = {
    STUDY_CONFIG_SNAPSHOT_FILENAME: "configuration_snapshot",
    SOURCE_CONFIG_SNAPSHOT_FILENAME: "source_configuration_snapshot",
    SOURCE_MANIFEST_SNAPSHOT_FILENAME: "source_manifest_snapshot",
    GENERATION_PROVENANCE_FILENAME: "generation_provenance",
    OBSERVABLE_TRIALS_FILENAME: "observable_trials",
    OBSERVABLE_AGGREGATE_FILENAME: "observable_aggregate",
    LEGACY_VERIFICATION_FILENAME: "legacy_metric_verification",
    OBSERVABLE_MAPS_FILENAME: "held_out_maps_and_coefficients",
}

_LEGACY_FLOAT_METRICS = (
    "data_jacobian_condition",
    "full_map_relative_l2_error",
    "supported_band_relative_l2_error",
    "integrated_density_relative_error",
    "centroid_y_error_um",
    "centroid_z_error_um",
    "rms_y_relative_error",
    "rms_z_relative_error",
)
_LEGACY_INTEGER_METRICS = ("data_jacobian_rank", "parameter_count")


@dataclass(frozen=True)
class SourceBenchmarkArtifacts:
    """Validated sealed source benchmark used by the held-out replay."""

    directory: Path
    manifest_path: Path
    config_path: Path
    held_out_trials_path: Path
    manifest: dict[str, Any] = field(repr=False)
    summary: dict[str, Any] = field(repr=False)
    metadata: dict[str, Any] = field(repr=False)
    config: dict[str, Any] = field(repr=False)
    held_out_rows: tuple[dict[str, str], ...] = field(repr=False)

    @property
    def run_id(self) -> str:
        return str(self.manifest["run_id"])


@dataclass(frozen=True)
class ObservableReplayTrial:
    """One verified frozen fit and its truth/reconstruction observables."""

    readout: ReadoutName
    fluence_mw_us: float
    morphology_name: str
    realization_index: int
    seed: int
    candidate_label: str
    truth_column_density_m2: NDArray[np.floating] = field(repr=False)
    reconstructed_column_density_m2: NDArray[np.floating] = field(repr=False)
    fitted_coefficients: NDArray[np.floating] = field(repr=False)
    truth_observables: DensityObservableSummary = field(repr=False)
    reconstructed_observables: DensityObservableSummary = field(repr=False)
    legacy_metrics: DensityRecoveryMetrics
    legacy_verification_rows: tuple[dict[str, Any], ...] = field(repr=False)


@dataclass(frozen=True)
class ObservableBenchmarkRun:
    """Complete verified replay before serialization."""

    study_config: dict[str, Any] = field(repr=False)
    source: SourceBenchmarkArtifacts = field(repr=False)
    context: MorphologyStudyContext = field(repr=False)
    integration_support: ObservableIntegrationSupport = field(repr=False)
    trials: tuple[ObservableReplayTrial, ...] = field(repr=False)
    elapsed_seconds: float


def _safe_artifact_path(directory: Path, filename: str) -> Path:
    path = Path(filename)
    if path.name != filename or path.is_absolute():
        raise ValueError(f"artifact name is not a safe basename: {filename!r}")
    return directory / path


def _deterministic_source_run_id(
    config_sha256: str,
    raw_hashes: dict[str, str],
) -> str:
    payload = json.dumps(
        {
            "config_sha256": config_sha256,
            "raw_artifact_sha256": raw_hashes,
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _source_artifact_records(
    source_directory: Path,
    manifest: dict[str, Any],
) -> tuple[dict[str, dict[str, str]], dict[str, Path]]:
    if manifest.get("schema_version") != SOURCE_MANIFEST_SCHEMA_VERSION:
        raise ValueError("unsupported source artifact-manifest schema version")
    records = manifest.get("artifacts")
    if not isinstance(records, dict):
        raise ValueError("source artifact manifest has no artifact records")
    paths: dict[str, Path] = {}
    normalised: dict[str, dict[str, str]] = {}
    for filename, raw_record in records.items():
        if not isinstance(raw_record, dict):
            raise ValueError(f"invalid source artifact record for {filename!r}")
        digest = str(raw_record.get("sha256", ""))
        role = str(raw_record.get("role", ""))
        if len(digest) != 64 or any(character not in "0123456789abcdef" for character in digest):
            raise ValueError(f"invalid source artifact hash for {filename!r}")
        path = _safe_artifact_path(source_directory, str(filename))
        if not path.is_file():
            raise FileNotFoundError(path)
        actual = file_sha256(path)
        if actual != digest:
            raise ValueError(
                f"source artifact hash mismatch for {filename}: expected {digest}, got {actual}"
            )
        paths[str(filename)] = path
        normalised[str(filename)] = {"role": role, "sha256": digest}
    return normalised, paths


def _parse_source_trial_success(value: str) -> bool:
    normalised = str(value).strip().lower()
    if normalised in {"true", "1", "yes"}:
        return True
    if normalised in {"false", "0", "no"}:
        return False
    raise ValueError(f"invalid source trial success value: {value!r}")


def _historical_windows_text_sha256(path: Path) -> str:
    """Reproduce pre-LF-contract hashes without weakening content checks."""

    content = path.read_bytes()
    normalised = content.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(normalised.replace(b"\n", b"\r\n")).hexdigest()


def _validate_critical_source_hashes(
    study_config: dict[str, Any],
    metadata: dict[str, Any],
    repository_root: Path,
) -> None:
    recorded = metadata.get("source_files_sha256")
    if not isinstance(recorded, dict):
        raise ValueError("source metadata has no source_files_sha256 record")
    critical = study_config["source_benchmark"].get("critical_source_files")
    if not isinstance(critical, list) or not critical:
        raise ValueError("observable study must declare non-empty critical_source_files")
    for relative in critical:
        key = Path(str(relative)).as_posix()
        expected = recorded.get(key)
        if not isinstance(expected, str):
            raise ValueError(f"source metadata has no critical hash for {key}")
        path = repository_root / Path(key)
        if not path.is_file():
            raise FileNotFoundError(path)
        actual = file_sha256(path)
        historical_windows_hash = _historical_windows_text_sha256(path)
        if actual != expected and historical_windows_hash != expected:
            raise ValueError(
                f"critical replay source hash mismatch for {key}: expected {expected}, got {actual}"
            )


def _validate_source_trial_contract(
    source_config: dict[str, Any],
    source_summary: dict[str, Any],
    rows: tuple[dict[str, str], ...],
    *,
    expected_count: int,
    require_success: bool,
) -> None:
    fluences = tuple(float(value) for value in source_config["ensemble"]["held_out_fluence_mw_us"])
    morphologies = tuple(
        str(value)
        for value in source_config["synthetic_morphologies"]["held_out_names"]
    )
    realization_count = int(source_config["ensemble"]["held_out_realizations_per_morphology"])
    expected_keys = {
        (readout, fluence, morphology, realization)
        for readout in READOUT_NAMES
        for fluence in fluences
        for morphology in morphologies
        for realization in range(realization_count)
    }
    if expected_count != len(expected_keys):
        raise ValueError(
            "configured expected source trial count disagrees with the frozen ensemble axes"
        )
    if len(rows) != expected_count:
        raise ValueError(
            f"source held-out table contains {len(rows)} trials, expected {expected_count}"
        )
    actual_keys: set[tuple[str, float, str, int]] = set()
    selected = source_summary.get("selected_candidates")
    held_summaries = source_summary.get("held_out_summaries")
    if not isinstance(selected, dict) or not isinstance(held_summaries, dict):
        raise ValueError("source summary lacks selected candidates or held-out summaries")
    summary_observations: dict[tuple[str, float], set[tuple[str, int, int]]] = {}
    for readout in READOUT_NAMES:
        readout_summary = held_summaries.get(readout)
        if not isinstance(readout_summary, dict):
            raise ValueError(f"source summary lacks held-out record for {readout}")
        for fluence in fluences:
            record = readout_summary.get(f"{fluence:g}")
            if not isinstance(record, dict):
                raise ValueError(
                    f"source summary lacks held-out record for {readout} at F={fluence:g}"
                )
            summary_observations[(readout, fluence)] = {
                (str(item[0]), int(item[1]), int(item[2]))
                for item in record.get("observation_keys", [])
            }
    row_observations: dict[tuple[str, float], set[tuple[str, int, int]]] = {
        key: set() for key in summary_observations
    }
    for row in rows:
        readout = str(row["readout"])
        fluence = float(row["fluence_mw_us"])
        morphology = str(row["morphology"])
        realization = int(row["realization_index"])
        seed = int(row["seed"])
        key = (readout, fluence, morphology, realization)
        if key in actual_keys:
            raise ValueError(f"duplicate source held-out trial key: {key}")
        actual_keys.add(key)
        if require_success and not _parse_source_trial_success(row["success"]):
            raise ValueError(f"source held-out trial was not successful: {key}")
        candidate = selected.get(readout)
        if not isinstance(candidate, dict) or row["candidate"] != str(candidate.get("label")):
            raise ValueError(f"source candidate mismatch for held-out trial: {key}")
        row_observations.setdefault((readout, fluence), set()).add(
            (morphology, realization, seed)
        )
    if actual_keys != expected_keys:
        missing = sorted(expected_keys.difference(actual_keys))
        unexpected = sorted(actual_keys.difference(expected_keys))
        raise ValueError(
            "source held-out trial axes disagree with frozen config; "
            f"missing={missing}, unexpected={unexpected}"
        )
    for key, expected in summary_observations.items():
        if row_observations.get(key) != expected:
            raise ValueError(
                f"source held-out CSV seeds disagree with benchmark summary for {key}"
            )


def validate_source_benchmark_artifacts(
    study_config: dict[str, Any],
    repository_root: Path,
) -> SourceBenchmarkArtifacts:
    """Validate hashes, run identity, source code and all 60 trial keys."""

    repository_root = repository_root.resolve()
    source_contract = study_config["source_benchmark"]
    configured = Path(str(source_contract["directory"]))
    source_directory = (
        configured.resolve()
        if configured.is_absolute()
        else (repository_root / configured).resolve()
    )
    manifest_path = source_directory / "artifact_manifest.json"
    metadata_path = source_directory / "metadata.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)
    if not metadata_path.is_file():
        raise FileNotFoundError(metadata_path)
    manifest = load_json(manifest_path)
    records, paths = _source_artifact_records(source_directory, manifest)
    expected_label = str(source_contract["expected_label"])
    expected_run_id = str(source_contract["expected_run_id"])
    if manifest.get("label") != expected_label:
        raise ValueError("source artifact label does not match observable study config")
    if manifest.get("run_id") != expected_run_id:
        raise ValueError("source artifact run_id does not match observable study config")
    config_record = manifest.get("config")
    if not isinstance(config_record, dict):
        raise ValueError("source artifact manifest has no config record")
    config_filename = str(config_record.get("artifact", ""))
    config_path = paths.get(config_filename)
    if config_path is None:
        raise ValueError("source config artifact is absent from source manifest records")
    config_sha256 = str(config_record.get("sha256", ""))
    if config_sha256 != records[config_filename]["sha256"]:
        raise ValueError("source config record and artifact hash disagree")
    raw_hashes = {
        filename: record["sha256"]
        for filename, record in records.items()
        if record["role"] in _SOURCE_RAW_ROLES
    }
    present_raw_roles = {
        record["role"]
        for record in records.values()
        if record["role"] in _SOURCE_RAW_ROLES
    }
    if present_raw_roles != _SOURCE_RAW_ROLES:
        raise ValueError("source manifest does not contain every required raw artifact role")
    recomputed_run_id = _deterministic_source_run_id(config_sha256, raw_hashes)
    if recomputed_run_id != expected_run_id:
        raise ValueError("source artifact run_id does not match its sealed raw hashes")
    summary_paths = [
        paths[name]
        for name, record in records.items()
        if record["role"] == "run_summary"
    ]
    held_paths = [
        paths[name]
        for name, record in records.items()
        if record["role"] == "held_out_trials"
    ]
    if len(summary_paths) != 1 or len(held_paths) != 1:
        raise ValueError("source manifest must name one summary and one held-out table")
    summary = load_json(summary_paths[0])
    metadata = load_json(metadata_path)
    source_config = load_json(config_path)
    if summary.get("run_id") != expected_run_id or metadata.get("run_id") != expected_run_id:
        raise ValueError("source summary/metadata run_id disagrees with source manifest")
    if summary.get("label") != expected_label or source_config.get("label") != expected_label:
        raise ValueError("source summary/config label disagrees with source manifest")
    if metadata.get("config_sha256") != config_sha256:
        raise ValueError("source metadata config hash disagrees with source manifest")
    _validate_critical_source_hashes(study_config, metadata, repository_root)
    held_rows = tuple(load_rows(held_paths[0]))
    verification = study_config["replay_verification"]
    _validate_source_trial_contract(
        source_config,
        summary,
        held_rows,
        expected_count=int(source_contract["expected_held_out_trial_count"]),
        require_success=bool(verification["require_all_source_trials_successful"]),
    )
    return SourceBenchmarkArtifacts(
        directory=source_directory,
        manifest_path=manifest_path,
        config_path=config_path,
        held_out_trials_path=held_paths[0],
        manifest=manifest,
        summary=summary,
        metadata=metadata,
        config=source_config,
        held_out_rows=held_rows,
    )


def _selected_candidate(
    context: MorphologyStudyContext,
    source_summary: dict[str, Any],
    readout: ReadoutName,
) -> ReconstructionCandidate:
    selected = source_summary["selected_candidates"][readout]
    label = str(selected["label"])
    matches = [candidate for candidate in context.candidates if candidate.label == label]
    if len(matches) != 1:
        raise ValueError(f"expected exactly one frozen candidate labelled {label!r}")
    candidate = matches[0]
    if candidate.model.parameter_count != int(selected["parameter_count"]):
        raise ValueError(f"frozen candidate parameter count disagrees for {readout}")
    return candidate


def build_observable_integration_support(
    study_config: dict[str, Any],
    source: SourceBenchmarkArtifacts,
    context: MorphologyStudyContext,
) -> ObservableIntegrationSupport:
    """Build the declared object-grid support and match it to both frozen models."""

    support_config = study_config["integration_support"]
    if support_config.get("kind") != "source_density_basis_rectangular_support":
        raise ValueError("unsupported observable integration-support contract")
    density_config = source.config["density_basis"]
    half_y_um = float(density_config["support_half_width_y_um"])
    half_z_um = float(density_config["support_half_width_z_um"])
    support_mask = (
        (np.abs(context.grid.y_grid_m * 1e6) <= half_y_um)
        & (np.abs(context.grid.z_grid_m * 1e6) <= half_z_um)
    )
    support = ObservableIntegrationSupport(
        y_grid_m=context.grid.y_grid_m,
        z_grid_m=context.grid.z_grid_m,
        support_mask=support_mask,
    )
    for readout in READOUT_NAMES:
        candidate = _selected_candidate(context, source.summary, readout)
        if not np.array_equal(candidate.model.support_mask, support.support_mask):
            raise ValueError(
                f"observable integration support does not match frozen {readout} model support"
            )
    return support


def _fit_candidate(
    context: MorphologyStudyContext,
    readout: ReadoutName,
    measurement: DifferentiableDensityMeasurement,
    candidate: ReconstructionCandidate,
    channels: tuple[NDArray[np.floating], ...],
) -> DensityFitResult:
    fit_config = context.config["fit"]
    if readout == "dual_port":
        initial = linearised_nonnegative_initialisation(
            measurement,
            candidate.model,
            channels,
            coefficient_upper=candidate.coefficient_upper,
            ridge_strength=float(
                fit_config["dual_port_initialisation_ridge_strength"]
            ),
        ).coefficients
    else:
        if not isinstance(measurement, DarkFieldFaradayMeasurement):
            raise TypeError("dark-field replay requires a dark-field measurement")
        if len(channels) != 1:
            raise ValueError("dark-field replay requires exactly one raw channel")
        initial = dark_field_sqrt_moment_initialisation(
            measurement,
            candidate.model,
            channels[0],
            smooth_bounds=context.smooth_bounds,
            coefficient_upper=candidate.coefficient_upper,
            projection_ridge_strength=float(
                fit_config["dark_field_projection_ridge_strength"]
            ),
        ).coefficients
    return fit_nonnegative_basis_density(
        measurement,
        candidate.model,
        channels,
        initial_coefficients=initial,
        coefficient_upper=candidate.coefficient_upper,
        regularisation=candidate.regularisation,
        options=candidate.fit_options,
    )


def _legacy_metric_values(
    fit: DensityFitResult,
    metrics: DensityRecoveryMetrics,
) -> dict[str, float | int]:
    return {
        "data_jacobian_rank": fit.diagnostics.data_jacobian_rank,
        "data_jacobian_condition": fit.diagnostics.data_jacobian_condition,
        "parameter_count": int(fit.coefficients.size),
        "full_map_relative_l2_error": metrics.full_map_relative_l2_error,
        "supported_band_relative_l2_error": metrics.supported_band_relative_l2_error,
        "integrated_density_relative_error": metrics.integrated_density_relative_error,
        "centroid_y_error_um": metrics.centroid_y_error_um,
        "centroid_z_error_um": metrics.centroid_z_error_um,
        "rms_y_relative_error": metrics.rms_y_relative_error,
        "rms_z_relative_error": metrics.rms_z_relative_error,
    }


def verify_legacy_trial_metrics(
    source_row: dict[str, str],
    fit: DensityFitResult,
    metrics: DensityRecoveryMetrics,
    *,
    relative_tolerance: float,
    absolute_tolerance: float,
) -> tuple[dict[str, Any], ...]:
    """Verify every stored numerical diagnostic before accepting a replayed map."""

    replayed = _legacy_metric_values(fit, metrics)
    prefix = {
        "readout": source_row["readout"],
        "fluence_mw_us": float(source_row["fluence_mw_us"]),
        "morphology": source_row["morphology"],
        "realization_index": int(source_row["realization_index"]),
        "seed": int(source_row["seed"]),
    }
    rows: list[dict[str, Any]] = []
    failures: list[str] = []
    for metric in (*_LEGACY_INTEGER_METRICS, *_LEGACY_FLOAT_METRICS):
        if metric in _LEGACY_INTEGER_METRICS:
            source_value: float | int = int(source_row[metric])
            replayed_value: float | int = int(replayed[metric])
            passed = source_value == replayed_value
        else:
            source_value = float(source_row[metric])
            replayed_value = float(replayed[metric])
            passed = bool(
                np.isclose(
                    replayed_value,
                    source_value,
                    rtol=relative_tolerance,
                    atol=absolute_tolerance,
                )
            )
        absolute_difference = abs(float(replayed_value) - float(source_value))
        denominator = max(abs(float(source_value)), np.finfo(float).eps)
        relative_difference = absolute_difference / denominator
        rows.append(
            {
                **prefix,
                "metric": metric,
                "source_value": source_value,
                "replayed_value": replayed_value,
                "absolute_difference": absolute_difference,
                "relative_difference": relative_difference,
                "within_tolerance": passed,
            }
        )
        if not passed:
            failures.append(metric)
    if failures:
        key = (
            prefix["readout"],
            prefix["fluence_mw_us"],
            prefix["morphology"],
            prefix["realization_index"],
        )
        raise ValueError(
            f"replayed legacy metrics disagree with sealed source row {key}: "
            + ", ".join(failures)
        )
    return tuple(rows)


def _source_rows_by_key(
    source: SourceBenchmarkArtifacts,
) -> dict[tuple[str, float, str, int], dict[str, str]]:
    return {
        (
            row["readout"],
            float(row["fluence_mw_us"]),
            row["morphology"],
            int(row["realization_index"]),
        ): row
        for row in source.held_out_rows
    }


def run_observable_benchmark_study(
    study_config: dict[str, Any],
    source: SourceBenchmarkArtifacts,
    *,
    progress: Any | None = None,
) -> ObservableBenchmarkRun:
    """Replay only the frozen held-out fits and extract physical observables."""

    started = perf_counter()
    context = build_morphology_study_context(source.config)
    integration_support = build_observable_integration_support(
        study_config,
        source,
        context,
    )
    angle_threshold = float(
        study_config["integration_support"]["angle_anisotropy_threshold"]
    )
    minimum_response = float(
        study_config["integration_support"]["minimum_integrated_response"]
    )
    verification = study_config["replay_verification"]
    relative_tolerance = float(verification["relative_tolerance"])
    absolute_tolerance = float(verification["absolute_tolerance"])
    source_rows = _source_rows_by_key(source)
    ensemble = source.config["ensemble"]
    fluences = tuple(float(value) for value in ensemble["held_out_fluence_mw_us"])
    realization_count = int(ensemble["held_out_realizations_per_morphology"])
    trials: list[ObservableReplayTrial] = []
    truth_summaries = {
        morphology.name: extract_density_observables(
            morphology.column_density_m2,
            integration_support,
            angle_anisotropy_threshold=angle_threshold,
            minimum_integrated_response=minimum_response,
        )
        for morphology in context.morphology_split.held_out
    }
    for mode_index, readout in enumerate(READOUT_NAMES):
        candidate = _selected_candidate(context, source.summary, readout)
        for fluence_index, fluence in enumerate(fluences):
            if progress is not None:
                progress(f"{readout}: replaying frozen held-out fits at F={fluence:g} mW us")
            measurement = make_study_measurement(
                context,
                readout,
                fluence_mw_us=fluence,
            )
            observations = generate_noisy_observation_ensemble(
                measurement,
                context.morphology_split.held_out,
                realizations_per_morphology=realization_count,
                base_seed=(
                    int(ensemble["held_out_seed"])
                    + 10000 * mode_index
                    + 100 * fluence_index
                ),
            )
            for observation in observations:
                key = (
                    readout,
                    fluence,
                    observation.morphology.name,
                    observation.realization_index,
                )
                source_row = source_rows[key]
                if observation.seed != int(source_row["seed"]):
                    raise ValueError(f"regenerated observation seed disagrees for {key}")
                if source_row["candidate"] != candidate.label:
                    raise ValueError(f"frozen candidate label disagrees for {key}")
                fit = _fit_candidate(
                    context,
                    readout,
                    measurement,
                    candidate,
                    observation.channels,
                )
                if not fit.diagnostics.success:
                    raise RuntimeError(f"frozen held-out replay fit failed for {key}")
                metrics = density_recovery_metrics(
                    observation.morphology.column_density_m2,
                    fit.column_density_m2,
                    context.grid,
                )
                verification_rows = verify_legacy_trial_metrics(
                    source_row,
                    fit,
                    metrics,
                    relative_tolerance=relative_tolerance,
                    absolute_tolerance=absolute_tolerance,
                )
                reconstructed_summary = extract_density_observables(
                    fit.column_density_m2,
                    integration_support,
                    angle_anisotropy_threshold=angle_threshold,
                    minimum_integrated_response=minimum_response,
                )
                trials.append(
                    ObservableReplayTrial(
                        readout=readout,
                        fluence_mw_us=fluence,
                        morphology_name=observation.morphology.name,
                        realization_index=observation.realization_index,
                        seed=observation.seed,
                        candidate_label=candidate.label,
                        truth_column_density_m2=np.array(
                            observation.morphology.column_density_m2,
                            dtype=float,
                            copy=True,
                        ),
                        reconstructed_column_density_m2=np.array(
                            fit.column_density_m2,
                            dtype=float,
                            copy=True,
                        ),
                        fitted_coefficients=np.array(
                            fit.coefficients,
                            dtype=float,
                            copy=True,
                        ),
                        truth_observables=truth_summaries[observation.morphology.name],
                        reconstructed_observables=reconstructed_summary,
                        legacy_metrics=metrics,
                        legacy_verification_rows=verification_rows,
                    )
                )
    expected_count = int(
        study_config["source_benchmark"]["expected_held_out_trial_count"]
    )
    if len(trials) != expected_count:
        raise RuntimeError(
            f"observable replay produced {len(trials)} trials, expected {expected_count}"
        )
    return ObservableBenchmarkRun(
        study_config=study_config,
        source=source,
        context=context,
        integration_support=integration_support,
        trials=tuple(trials),
        elapsed_seconds=perf_counter() - started,
    )


def _optional_scaled(value: float | None, scale: float = 1.0) -> float | str:
    return "" if value is None else float(value) * scale


def _observable_fields(
    prefix: str,
    summary: DensityObservableSummary,
) -> dict[str, Any]:
    covariance = summary.covariance_m2
    flags = summary.support_flags
    moments_supported = flags.moments_numerically_supported
    positive_response = flags.positive_integrated_response
    return {
        f"{prefix}_integrated_response": summary.integrated_response,
        f"{prefix}_centroid_y_um": _optional_scaled(summary.centroid_y_m, 1e6),
        f"{prefix}_centroid_z_um": _optional_scaled(summary.centroid_z_m, 1e6),
        f"{prefix}_covariance_yy_um2": (
            "" if covariance is None else float(covariance[0, 0]) * 1e12
        ),
        f"{prefix}_covariance_yz_um2": (
            "" if covariance is None else float(covariance[0, 1]) * 1e12
        ),
        f"{prefix}_covariance_zz_um2": (
            "" if covariance is None else float(covariance[1, 1]) * 1e12
        ),
        f"{prefix}_major_rms_width_um": _optional_scaled(
            summary.major_rms_width_m,
            1e6,
        ),
        f"{prefix}_minor_rms_width_um": _optional_scaled(
            summary.minor_rms_width_m,
            1e6,
        ),
        f"{prefix}_aspect_ratio": _optional_scaled(summary.aspect_ratio),
        f"{prefix}_principal_axis_angle_rad": _optional_scaled(
            summary.principal_axis_angle_rad
        ),
        f"{prefix}_fractional_anisotropy": _optional_scaled(
            summary.fractional_anisotropy
        ),
        f"{prefix}_positive_integrated_response": bool(positive_response),
        f"{prefix}_moments_numerically_supported": bool(moments_supported),
        f"{prefix}_centroid_supported": flags.centroid_supported,
        f"{prefix}_covariance_supported": flags.covariance_supported,
        f"{prefix}_widths_supported": flags.widths_supported,
        f"{prefix}_aspect_ratio_supported": flags.aspect_ratio_supported,
        f"{prefix}_principal_axis_angle_supported": (
            flags.principal_axis_angle_supported
        ),
        f"{prefix}_support_reasons": "|".join(flags.reasons),
    }


def _axial_angle_difference_rad(reconstructed: float, truth: float) -> float:
    return float((reconstructed - truth + 0.5 * np.pi) % np.pi - 0.5 * np.pi)


def _trial_row(trial: ObservableReplayTrial) -> dict[str, Any]:
    truth = trial.truth_observables
    reconstructed = trial.reconstructed_observables
    truth_covariance = truth.covariance_m2
    reconstructed_covariance = reconstructed.covariance_m2
    centroid_error_y_um: float | str = ""
    centroid_error_z_um: float | str = ""
    centroid_position_error_um: float | str = ""
    if truth.centroid_m is not None and reconstructed.centroid_m is not None:
        centroid_delta_um = (reconstructed.centroid_m - truth.centroid_m) * 1e6
        centroid_error_y_um = float(centroid_delta_um[0])
        centroid_error_z_um = float(centroid_delta_um[1])
        centroid_position_error_um = float(np.linalg.norm(centroid_delta_um))
    covariance_error = ""
    covariance_component_errors: tuple[float | str, float | str, float | str] = (
        "",
        "",
        "",
    )
    if truth_covariance is not None and reconstructed_covariance is not None:
        covariance_delta = reconstructed_covariance - truth_covariance
        denominator = max(
            float(np.linalg.norm(truth_covariance)),
            np.finfo(float).eps,
        )
        covariance_error = float(np.linalg.norm(covariance_delta) / denominator)
        covariance_component_errors = (
            float(covariance_delta[0, 0]) * 1e12,
            float(covariance_delta[0, 1]) * 1e12,
            float(covariance_delta[1, 1]) * 1e12,
        )

    def relative_error(recovered: float | None, expected: float | None) -> float | str:
        if recovered is None or expected is None or expected == 0.0:
            return ""
        return float(recovered / expected - 1.0)

    angle_error: float | str = ""
    if (
        truth.principal_axis_angle_rad is not None
        and reconstructed.principal_axis_angle_rad is not None
    ):
        angle_error = _axial_angle_difference_rad(
            reconstructed.principal_axis_angle_rad,
            truth.principal_axis_angle_rad,
        )
    return {
        "readout": trial.readout,
        "fluence_mw_us": trial.fluence_mw_us,
        "candidate": trial.candidate_label,
        "morphology": trial.morphology_name,
        "realization_index": trial.realization_index,
        "seed": trial.seed,
        **_observable_fields("truth", truth),
        **_observable_fields("reconstructed", reconstructed),
        "integrated_response_relative_error": (
            reconstructed.integrated_response / truth.integrated_response - 1.0
        ),
        "centroid_y_error_um": centroid_error_y_um,
        "centroid_z_error_um": centroid_error_z_um,
        "centroid_position_error_um": centroid_position_error_um,
        "covariance_yy_error_um2": covariance_component_errors[0],
        "covariance_yz_error_um2": covariance_component_errors[1],
        "covariance_zz_error_um2": covariance_component_errors[2],
        "covariance_tensor_relative_frobenius_error": covariance_error,
        "major_rms_width_relative_error": relative_error(
            reconstructed.major_rms_width_m,
            truth.major_rms_width_m,
        ),
        "minor_rms_width_relative_error": relative_error(
            reconstructed.minor_rms_width_m,
            truth.minor_rms_width_m,
        ),
        "aspect_ratio_relative_error": relative_error(
            reconstructed.aspect_ratio,
            truth.aspect_ratio,
        ),
        "principal_axis_axial_error_rad": angle_error,
        "legacy_metrics_verified": all(
            bool(row["within_tolerance"])
            for row in trial.legacy_verification_rows
        ),
    }


def _finite_values(rows: list[dict[str, Any]], key: str) -> NDArray[np.floating]:
    values = [
        float(row[key])
        for row in rows
        if row.get(key, "") != "" and np.isfinite(float(row[key]))
    ]
    return np.asarray(values, dtype=float)


def _distribution_fields(
    rows: list[dict[str, Any]],
    source_key: str,
    output_stem: str,
    *,
    absolute: bool = False,
    scale: float = 1.0,
) -> dict[str, Any]:
    values = _finite_values(rows, source_key)
    if absolute:
        values = np.abs(values)
    values = values * scale
    if not values.size:
        return {
            f"supported_{output_stem}_trials": 0,
            f"minimum_{output_stem}": "",
            f"median_{output_stem}": "",
            f"maximum_{output_stem}": "",
        }
    return {
        f"supported_{output_stem}_trials": int(values.size),
        f"minimum_{output_stem}": float(np.min(values)),
        f"median_{output_stem}": float(np.median(values)),
        f"maximum_{output_stem}": float(np.max(values)),
    }


def aggregate_observable_trial_rows(
    rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Aggregate the finite development ensemble by readout and fluence."""

    output: list[dict[str, Any]] = []
    readouts = [name for name in READOUT_NAMES if any(row["readout"] == name for row in rows)]
    for readout in readouts:
        fluences = sorted(
            {float(row["fluence_mw_us"]) for row in rows if row["readout"] == readout}
        )
        for fluence in fluences:
            group = [
                row
                for row in rows
                if row["readout"] == readout
                and float(row["fluence_mw_us"]) == fluence
            ]
            output.append(
                {
                    "readout": readout,
                    "fluence_mw_us": fluence,
                    "trial_count": len(group),
                    "legacy_verified_trials": sum(
                        bool(row["legacy_metrics_verified"]) for row in group
                    ),
                    "reconstructed_moment_supported_trials": sum(
                        bool(row["reconstructed_moments_numerically_supported"])
                        for row in group
                    ),
                    "principal_axis_angle_supported_pairs": sum(
                        bool(row["truth_principal_axis_angle_supported"])
                        and bool(row["reconstructed_principal_axis_angle_supported"])
                        for row in group
                    ),
                    **_distribution_fields(
                        group,
                        "integrated_response_relative_error",
                        "absolute_integrated_response_relative_error",
                        absolute=True,
                    ),
                    **_distribution_fields(
                        group,
                        "centroid_position_error_um",
                        "centroid_position_error_um",
                    ),
                    **_distribution_fields(
                        group,
                        "covariance_tensor_relative_frobenius_error",
                        "covariance_tensor_relative_frobenius_error",
                    ),
                    **_distribution_fields(
                        group,
                        "major_rms_width_relative_error",
                        "absolute_major_rms_width_relative_error",
                        absolute=True,
                    ),
                    **_distribution_fields(
                        group,
                        "minor_rms_width_relative_error",
                        "absolute_minor_rms_width_relative_error",
                        absolute=True,
                    ),
                    **_distribution_fields(
                        group,
                        "aspect_ratio_relative_error",
                        "absolute_aspect_ratio_relative_error",
                        absolute=True,
                    ),
                    **_distribution_fields(
                        group,
                        "principal_axis_axial_error_rad",
                        "absolute_principal_axis_angle_error_deg",
                        absolute=True,
                        scale=180.0 / np.pi,
                    ),
                }
            )
    return output


def _write_held_out_maps(path: Path, run: ObservableBenchmarkRun) -> None:
    source_config = run.source.config
    readout_names = np.asarray(READOUT_NAMES)
    fluences = np.asarray(
        source_config["ensemble"]["held_out_fluence_mw_us"],
        dtype=float,
    )
    morphology_names = np.asarray(
        source_config["synthetic_morphologies"]["held_out_names"]
    )
    realization_indices = np.arange(
        int(source_config["ensemble"]["held_out_realizations_per_morphology"]),
        dtype=int,
    )
    readout_index = {str(value): index for index, value in enumerate(readout_names)}
    fluence_index = {float(value): index for index, value in enumerate(fluences)}
    morphology_index = {
        str(value): index for index, value in enumerate(morphology_names)
    }
    realization_index = {
        int(value): index for index, value in enumerate(realization_indices)
    }
    leading_shape = (
        readout_names.size,
        fluences.size,
        morphology_names.size,
        realization_indices.size,
    )
    grid_shape = run.integration_support.shape
    parameter_counts = np.asarray(
        [
            _selected_candidate(run.context, run.source.summary, readout).model.parameter_count
            for readout in READOUT_NAMES
        ],
        dtype=int,
    )
    maximum_parameter_count = int(np.max(parameter_counts))
    reconstructed = np.full((*leading_shape, *grid_shape), np.nan, dtype=float)
    coefficients = np.full((*leading_shape, maximum_parameter_count), np.nan, dtype=float)
    seeds = np.full(leading_shape, -1, dtype=np.int64)
    filled = np.zeros(leading_shape, dtype=bool)
    truth = np.stack(
        [
            next(
                morphology.column_density_m2
                for morphology in run.context.morphology_split.held_out
                if morphology.name == name
            )
            for name in morphology_names
        ]
    )
    for trial in run.trials:
        index = (
            readout_index[trial.readout],
            fluence_index[trial.fluence_mw_us],
            morphology_index[trial.morphology_name],
            realization_index[trial.realization_index],
        )
        if filled[index]:
            raise ValueError(f"duplicate observable replay trial while writing maps: {index}")
        if not np.array_equal(
            trial.truth_column_density_m2,
            truth[index[2]],
        ):
            raise ValueError("trial truth map disagrees with frozen analytic morphology")
        reconstructed[index] = trial.reconstructed_column_density_m2
        coefficient_count = trial.fitted_coefficients.size
        if coefficient_count != parameter_counts[index[0]]:
            raise ValueError("trial coefficient count disagrees with frozen candidate")
        coefficients[index][:coefficient_count] = trial.fitted_coefficients
        seeds[index] = trial.seed
        filled[index] = True
    if not np.all(filled):
        raise ValueError("held-out map tensor has unfilled replay coordinates")
    np.savez_compressed(
        path,
        y_grid_m=run.integration_support.y_grid_m,
        z_grid_m=run.integration_support.z_grid_m,
        cell_area_m2=run.integration_support.cell_area_m2,
        support_mask=run.integration_support.support_mask,
        readout_names=readout_names,
        fluence_mw_us=fluences,
        morphology_names=morphology_names,
        realization_indices=realization_indices,
        candidate_labels=np.asarray(
            [
                _selected_candidate(run.context, run.source.summary, readout).label
                for readout in READOUT_NAMES
            ]
        ),
        coefficient_count_by_readout=parameter_counts,
        seeds=seeds,
        truth_column_density_m2=truth,
        reconstructed_column_density_m2=reconstructed,
        fitted_coefficients=coefficients,
        source_run_id=np.asarray(run.source.run_id),
        tensor_axis_order=np.asarray(
            "readout,fluence,morphology,realization,z_index,y_index"
        ),
    )


def _observable_raw_hashes(output_directory: Path) -> dict[str, str]:
    hashes: dict[str, str] = {}
    for filename in OBSERVABLE_RAW_ARTIFACT_ROLES:
        path = output_directory / filename
        if not path.is_file():
            raise FileNotFoundError(path)
        hashes[filename] = file_sha256(path)
    return hashes


def _deterministic_observable_run_id(
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


def _seal_observable_artifact_set(
    output_directory: Path,
    *,
    label: str,
    source_run_id: str,
    summary: dict[str, Any],
) -> tuple[str, Path, Path]:
    """Seal already-written raw artifacts and add run identity to the summary."""

    config_path = output_directory / STUDY_CONFIG_SNAPSHOT_FILENAME
    config_sha256 = file_sha256(config_path)
    raw_hashes = _observable_raw_hashes(output_directory)
    run_id = _deterministic_observable_run_id(
        config_sha256,
        source_run_id,
        raw_hashes,
    )
    summary_payload = dict(summary)
    summary_payload["run_id"] = run_id
    summary_path = output_directory / OBSERVABLE_SUMMARY_FILENAME
    write_json(summary_path, summary_payload)
    artifact_records = {
        filename: {
            "role": OBSERVABLE_RAW_ARTIFACT_ROLES[filename],
            "sha256": digest,
        }
        for filename, digest in raw_hashes.items()
    }
    artifact_records[summary_path.name] = {
        "role": "run_summary",
        "sha256": file_sha256(summary_path),
    }
    manifest_path = output_directory / OBSERVABLE_MANIFEST_FILENAME
    write_json(
        manifest_path,
        {
            "schema_version": OBSERVABLE_MANIFEST_SCHEMA_VERSION,
            "label": label,
            "run_id": run_id,
            "source_run_id": source_run_id,
            "config": {
                "artifact": config_path.name,
                "sha256": config_sha256,
            },
            "artifacts": artifact_records,
        },
    )
    return run_id, summary_path, manifest_path


def _observable_summary_payload(
    run: ObservableBenchmarkRun,
    trial_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    source_config = run.source.config
    density = source_config["density_basis"]
    parameter_counts = {
        readout: _selected_candidate(
            run.context,
            run.source.summary,
            readout,
        ).model.parameter_count
        for readout in READOUT_NAMES
    }
    return {
        "label": run.study_config["label"],
        "purpose": run.study_config["purpose"],
        "source_benchmark": {
            "label": run.source.manifest["label"],
            "run_id": run.source.run_id,
            "directory": run.study_config["source_benchmark"]["directory"],
            "manifest_sha256": file_sha256(run.source.manifest_path),
            "config_sha256": file_sha256(run.source.config_path),
            "held_out_trials_sha256": file_sha256(run.source.held_out_trials_path),
        },
        "replay": {
            "trial_count": len(run.trials),
            "all_legacy_metrics_verified": all(
                bool(row["legacy_metrics_verified"]) for row in trial_rows
            ),
            "readout_names": list(READOUT_NAMES),
            "fluence_mw_us": source_config["ensemble"]["held_out_fluence_mw_us"],
            "morphology_names": source_config["synthetic_morphologies"]["held_out_names"],
            "realization_indices": list(
                range(
                    int(
                        source_config["ensemble"][
                            "held_out_realizations_per_morphology"
                        ]
                    )
                )
            ),
            "selected_candidates": {
                readout: run.source.summary["selected_candidates"][readout]["label"]
                for readout in READOUT_NAMES
            },
            "parameter_count_by_readout": parameter_counts,
            "elapsed_seconds": run.elapsed_seconds,
        },
        "integration_support": {
            "kind": run.study_config["integration_support"]["kind"],
            "grid_shape": list(run.integration_support.shape),
            "support_half_width_y_um": float(density["support_half_width_y_um"]),
            "support_half_width_z_um": float(density["support_half_width_z_um"]),
            "supported_cell_count": int(
                np.count_nonzero(run.integration_support.support_mask)
            ),
            "physical_area_m2": run.integration_support.physical_area_m2,
            "angle_anisotropy_threshold": run.study_config["integration_support"][
                "angle_anisotropy_threshold"
            ],
            "minimum_integrated_response": run.study_config[
                "integration_support"
            ]["minimum_integrated_response"],
            "camera_fit_roi_is_not_integration_support": True,
        },
        "maps_schema": {
            "artifact": OBSERVABLE_MAPS_FILENAME,
            "truth_column_density_m2_shape": [
                len(source_config["synthetic_morphologies"]["held_out_names"]),
                *run.integration_support.shape,
            ],
            "reconstructed_column_density_m2_axis_order": [
                "readout",
                "fluence",
                "morphology",
                "realization",
                "z_index",
                "y_index",
            ],
            "reconstructed_column_density_m2_shape": [
                len(READOUT_NAMES),
                len(source_config["ensemble"]["held_out_fluence_mw_us"]),
                len(source_config["synthetic_morphologies"]["held_out_names"]),
                int(source_config["ensemble"]["held_out_realizations_per_morphology"]),
                *run.integration_support.shape,
            ],
            "coefficient_axis_is_padded_to_maximum_parameter_count": True,
        },
        "tables": {
            "trials": OBSERVABLE_TRIALS_FILENAME,
            "aggregate": OBSERVABLE_AGGREGATE_FILENAME,
            "legacy_verification": LEGACY_VERIFICATION_FILENAME,
        },
        "claims_boundary": run.study_config["claims_boundary"],
    }


def write_observable_benchmark_run(
    run: ObservableBenchmarkRun,
    study_config_path: Path,
    repository_root: Path,
    *,
    generation_provenance: dict[str, Any],
) -> dict[str, Path]:
    """Write and seal one independent observable held-out result family."""

    repository_root = repository_root.resolve()
    configured_output = Path(str(run.study_config["output_directory"]))
    output_directory = (
        configured_output.resolve()
        if configured_output.is_absolute()
        else (repository_root / configured_output).resolve()
    )
    if (
        output_directory == run.source.directory
        or run.source.directory in output_directory.parents
    ):
        raise ValueError(
            "observable outputs must not be written into the sealed source benchmark"
        )
    output_directory.mkdir(parents=True, exist_ok=True)
    config_snapshot = output_directory / STUDY_CONFIG_SNAPSHOT_FILENAME
    source_config_snapshot = output_directory / SOURCE_CONFIG_SNAPSHOT_FILENAME
    source_manifest_snapshot = output_directory / SOURCE_MANIFEST_SNAPSHOT_FILENAME
    shutil.copyfile(study_config_path, config_snapshot)
    shutil.copyfile(run.source.config_path, source_config_snapshot)
    shutil.copyfile(run.source.manifest_path, source_manifest_snapshot)
    if load_json(config_snapshot) != run.study_config:
        raise ValueError("observable config snapshot differs from in-memory study config")
    if load_json(source_config_snapshot) != run.source.config:
        raise ValueError("source config snapshot differs from validated source config")
    provenance_path = output_directory / GENERATION_PROVENANCE_FILENAME
    write_json(provenance_path, generation_provenance)
    trial_rows = [_trial_row(trial) for trial in run.trials]
    aggregate_rows = aggregate_observable_trial_rows(trial_rows)
    verification_rows = [
        row
        for trial in run.trials
        for row in trial.legacy_verification_rows
    ]
    trial_path = output_directory / OBSERVABLE_TRIALS_FILENAME
    aggregate_path = output_directory / OBSERVABLE_AGGREGATE_FILENAME
    verification_path = output_directory / LEGACY_VERIFICATION_FILENAME
    maps_path = output_directory / OBSERVABLE_MAPS_FILENAME
    write_rows(trial_path, trial_rows)
    write_rows(aggregate_path, aggregate_rows)
    write_rows(verification_path, verification_rows)
    _write_held_out_maps(maps_path, run)
    summary = _observable_summary_payload(run, trial_rows)
    _, summary_path, manifest_path = _seal_observable_artifact_set(
        output_directory,
        label=str(run.study_config["label"]),
        source_run_id=run.source.run_id,
        summary=summary,
    )
    return {
        "trials": trial_path,
        "aggregate": aggregate_path,
        "legacy_verification": verification_path,
        "maps": maps_path,
        "summary": summary_path,
        "manifest": manifest_path,
        "provenance": provenance_path,
    }


def run_and_write_observable_benchmark(
    study_config: dict[str, Any],
    study_config_path: Path,
    repository_root: Path,
    *,
    progress: Any | None = None,
) -> dict[str, Path]:
    """Validate, replay and seal the declared observable benchmark."""

    repository_root = repository_root.resolve()
    study_config_path = study_config_path.resolve()
    config_sha256 = file_sha256(study_config_path)
    source = validate_source_benchmark_artifacts(study_config, repository_root)
    before = capture_reconstruction_provenance(
        repository_root,
        entry_points=(
            Path("scripts/generate_reconstruction_observable_benchmark.py"),
        ),
    )
    run = run_observable_benchmark_study(
        study_config,
        source,
        progress=progress,
    )
    if file_sha256(study_config_path) != config_sha256:
        raise RuntimeError("observable study config changed during held-out replay")
    after = capture_reconstruction_provenance(
        repository_root,
        entry_points=(
            Path("scripts/generate_reconstruction_observable_benchmark.py"),
        ),
    )
    for field_name in ("source_files_sha256", "runtime_versions"):
        if before[field_name] != after[field_name]:
            raise RuntimeError(
                f"observable replay {field_name} changed during numerical generation"
            )
    generation_provenance = {
        "study_config_source": str(study_config_path),
        "study_config_source_sha256": config_sha256,
        "source_benchmark": {
            "run_id": source.run_id,
            "manifest": str(source.manifest_path),
            "manifest_sha256": file_sha256(source.manifest_path),
            "config_sha256": file_sha256(source.config_path),
            "held_out_trials_sha256": file_sha256(source.held_out_trials_path),
        },
        "generation_state_before_replay": before,
        "generation_state_after_replay": after,
        "legacy_metric_tolerances": {
            "relative": study_config["replay_verification"]["relative_tolerance"],
            "absolute": study_config["replay_verification"]["absolute_tolerance"],
        },
        "observable_contract": {
            "fixed_common_support": True,
            "cell_area_weighted": True,
            "common_response_scale_gain_and_optical_transfer_stability": (
                "assumptions required for relative integrated-signal interpretation; "
                "not verified by moment extraction"
            ),
        },
    }
    return write_observable_benchmark_run(
        run,
        study_config_path,
        repository_root,
        generation_provenance=generation_provenance,
    )


__all__ = [
    "GENERATION_PROVENANCE_FILENAME",
    "LEGACY_VERIFICATION_FILENAME",
    "OBSERVABLE_AGGREGATE_FILENAME",
    "OBSERVABLE_MANIFEST_FILENAME",
    "OBSERVABLE_MAPS_FILENAME",
    "OBSERVABLE_RAW_ARTIFACT_ROLES",
    "OBSERVABLE_SUMMARY_FILENAME",
    "OBSERVABLE_TRIALS_FILENAME",
    "ObservableBenchmarkRun",
    "ObservableReplayTrial",
    "SourceBenchmarkArtifacts",
    "aggregate_observable_trial_rows",
    "build_observable_integration_support",
    "run_and_write_observable_benchmark",
    "run_observable_benchmark_study",
    "validate_source_benchmark_artifacts",
    "verify_legacy_trial_metrics",
    "write_observable_benchmark_run",
]
