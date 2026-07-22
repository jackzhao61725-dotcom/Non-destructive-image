"""Targeted curvature-weight range check for a frozen reconstruction basis.

This study is deliberately narrower than the morphology benchmark.  It reuses
that benchmark's calibration morphologies, raw noisy observations and fitting
contract, changes only the declared curvature weights, and records whether the
selected weight was bracketed by sufficiently stronger regularisation.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path
import shutil
from time import perf_counter
from typing import Any, cast

import numpy as np

from ..ensemble import (
    CandidateEnsembleAssessment,
    assess_reconstruction_candidate,
    generate_noisy_observation_ensemble,
    make_dark_field_candidate_initialiser,
    make_linear_candidate_initialiser,
)
from ..measurements import DarkFieldFaradayMeasurement
from .io import file_sha256, write_json, write_rows
from .morphology import (
    MorphologyStudyContext,
    ProgressCallback,
    ReadoutName,
    build_morphology_study_context,
    make_study_measurement,
)
from .provenance import capture_reconstruction_provenance


@dataclass(frozen=True)
class CurvatureRangeCheckRun:
    """One completed calibration-only curvature-weight range check."""

    context: MorphologyStudyContext = field(repr=False)
    check_config: dict[str, Any] = field(repr=False)
    source_benchmark_label: str
    candidate_rows: tuple[dict[str, Any], ...]
    trial_rows: tuple[dict[str, Any], ...]
    elapsed_seconds: float


def _progress(callback: ProgressCallback | None, message: str) -> None:
    if callback is not None:
        callback(message)


def _checked_weights(check_config: dict[str, Any]) -> tuple[float, ...]:
    weights = tuple(float(value) for value in check_config["curvature_weights_um2"])
    if not weights:
        raise ValueError("curvature range check requires at least one weight")
    if any(not np.isfinite(weight) or weight < 0 for weight in weights):
        raise ValueError("curvature weights must be finite and non-negative")
    if len(weights) != len(set(weights)):
        raise ValueError("curvature weights must be unique")
    if tuple(sorted(weights)) != weights:
        raise ValueError("curvature weights must be in increasing order")
    return weights


def build_curvature_range_check_context(
    benchmark_config: dict[str, Any],
    check_config: dict[str, Any],
) -> MorphologyStudyContext:
    """Reuse a benchmark contract while replacing only basis and weight range."""

    derived = deepcopy(benchmark_config)
    basis_label = str(check_config["basis_label"])
    matching = [
        basis
        for basis in derived["density_basis"]["candidates"]
        if basis["basis_label"] == basis_label
    ]
    if len(matching) != 1:
        raise ValueError(
            f"expected exactly one benchmark basis labelled {basis_label!r}"
        )
    weights = _checked_weights(check_config)
    reference_weight = float(check_config["reference_weight_um2"])
    if reference_weight not in weights:
        raise ValueError("reference weight must be included in the checked range")
    source_weights = {float(value) for value in matching[0]["curvature_weights_um2"]}
    if reference_weight not in source_weights:
        raise ValueError("reference weight is not declared by the source benchmark")

    basis = deepcopy(matching[0])
    basis["curvature_weights_um2"] = list(weights)
    derived["density_basis"]["candidates"] = [basis]
    derived["label"] = str(check_config["label"])
    derived["output_directory"] = str(check_config["output_directory"])
    return build_morphology_study_context(derived)


def _candidate_weight(assessment: CandidateEnsembleAssessment) -> float:
    regularisation = assessment.candidate.regularisation
    return 0.0 if regularisation is None else float(regularisation.weight_um2)


def _candidate_row(
    readout: ReadoutName,
    basis_label: str,
    assessment: CandidateEnsembleAssessment,
) -> dict[str, Any]:
    summary = assessment.summary
    return {
        "readout": readout,
        "basis_label": basis_label,
        "candidate": summary.candidate_label,
        "curvature_weight_um2": _candidate_weight(assessment),
        "parameter_count": summary.parameter_count,
        "trial_count": summary.trial_count,
        "success_fraction": summary.success_fraction,
        "full_rank_fraction": summary.full_rank_fraction,
        "median_supported_band_relative_l2_error": summary.median_supported_band_error,
        "upper_quartile_supported_band_relative_l2_error": summary.upper_quartile_supported_band_error,
        "median_absolute_integrated_density_relative_error": summary.median_absolute_integrated_density_error,
        "median_data_jacobian_condition": summary.median_data_jacobian_condition,
        "maximum_data_jacobian_condition": summary.maximum_data_jacobian_condition,
    }


def _trial_rows(
    readout: ReadoutName,
    assessment: CandidateEnsembleAssessment,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    weight = _candidate_weight(assessment)
    for trial in assessment.trials:
        metrics = trial.metrics
        rows.append(
            {
                "readout": readout,
                "candidate": trial.candidate_label,
                "curvature_weight_um2": weight,
                "morphology": trial.morphology_name,
                "realization_index": trial.realization_index,
                "seed": trial.observation_seed,
                "success": trial.success,
                "data_jacobian_rank": trial.data_jacobian_rank,
                "data_jacobian_condition": trial.data_jacobian_condition,
                "parameter_count": trial.parameter_count,
                "supported_band_relative_l2_error": (
                    metrics.supported_band_relative_l2_error
                    if metrics is not None
                    else ""
                ),
                "integrated_density_relative_error": (
                    metrics.integrated_density_relative_error
                    if metrics is not None
                    else ""
                ),
                "message": trial.message,
            }
        )
    return rows


def run_curvature_range_check(
    benchmark_config: dict[str, Any],
    check_config: dict[str, Any],
    *,
    progress: ProgressCallback | None = None,
) -> CurvatureRangeCheckRun:
    """Run the declared weights on the source benchmark's shared calibration set."""

    started = perf_counter()
    context = build_curvature_range_check_context(benchmark_config, check_config)
    readout = str(check_config["readout"])
    if readout not in {"dual_port", "dark_field"}:
        raise ValueError("readout must be 'dual_port' or 'dark_field'")
    typed_readout = cast(ReadoutName, readout)
    mode_index = {"dual_port": 0, "dark_field": 1}[typed_readout]
    ensemble = context.config["ensemble"]
    base_seed = int(ensemble["calibration_seed"]) + 10000 * mode_index
    measurement = make_study_measurement(context, typed_readout)
    observations = generate_noisy_observation_ensemble(
        measurement,
        context.morphology_split.calibration,
        realizations_per_morphology=int(
            ensemble["calibration_realizations_per_morphology"]
        ),
        base_seed=base_seed,
    )
    fit = context.config["fit"]
    if typed_readout == "dual_port":
        initialiser = make_linear_candidate_initialiser(
            measurement,
            ridge_strength=float(fit["dual_port_initialisation_ridge_strength"]),
        )
    else:
        if not isinstance(measurement, DarkFieldFaradayMeasurement):
            raise TypeError("dark-field range check requires a dark-field measurement")
        initialiser = make_dark_field_candidate_initialiser(
            measurement,
            smooth_bounds=context.smooth_bounds,
            projection_ridge_strength=float(
                fit["dark_field_projection_ridge_strength"]
            ),
        )

    basis_label = str(check_config["basis_label"])
    candidate_rows: list[dict[str, Any]] = []
    trial_rows: list[dict[str, Any]] = []
    for index, candidate in enumerate(context.candidates, start=1):
        _progress(
            progress,
            f"{typed_readout}: range candidate {index}/{len(context.candidates)} "
            f"{candidate.label}",
        )
        assessment = assess_reconstruction_candidate(
            measurement,
            candidate,
            observations,
            initialise=initialiser,
        )
        candidate_rows.append(_candidate_row(typed_readout, basis_label, assessment))
        trial_rows.extend(_trial_rows(typed_readout, assessment))

    anchor = check_config.get("consistency_anchor")
    if anchor is not None:
        anchor_weight = float(anchor["curvature_weight_um2"])
        anchor_rows = [
            row
            for row in candidate_rows
            if float(row["curvature_weight_um2"]) == anchor_weight
        ]
        if len(anchor_rows) != 1:
            raise RuntimeError("consistency anchor does not identify one range result")
        observed = float(
            anchor_rows[0]["median_supported_band_relative_l2_error"]
        )
        expected = float(anchor["expected_median_supported_band_relative_l2_error"])
        tolerance = float(anchor["absolute_tolerance"])
        if not np.isclose(observed, expected, rtol=0.0, atol=tolerance):
            raise RuntimeError(
                "source-benchmark consistency anchor failed: "
                f"observed {observed:.17g}, expected {expected:.17g} "
                f"within {tolerance:.3g}"
            )

    return CurvatureRangeCheckRun(
        context=context,
        check_config=check_config,
        source_benchmark_label=str(benchmark_config["label"]),
        candidate_rows=tuple(candidate_rows),
        trial_rows=tuple(trial_rows),
        elapsed_seconds=perf_counter() - started,
    )


def _metadata_path(path: Path, repository_root: Path) -> str:
    try:
        return str(path.resolve().relative_to(repository_root.resolve()))
    except ValueError:
        return str(path.resolve())


def write_curvature_range_check_run(
    run: CurvatureRangeCheckRun,
    check_config_path: Path,
    benchmark_config_path: Path,
    repository_root: Path,
    *,
    generation_check_config_sha256: str,
    generation_benchmark_config_sha256: str,
    generation_provenance: dict[str, object],
) -> dict[str, Path]:
    """Serialize one range check with its two input-config hashes."""

    repository_root = repository_root.resolve()
    check_config_path = check_config_path.resolve()
    benchmark_config_path = benchmark_config_path.resolve()
    if not check_config_path.is_file():
        raise FileNotFoundError(check_config_path)
    if not benchmark_config_path.is_file():
        raise FileNotFoundError(benchmark_config_path)
    config_provenance = {
        "check_config_sha256": file_sha256(check_config_path),
        "source_benchmark_config_sha256": file_sha256(benchmark_config_path),
    }
    if config_provenance["check_config_sha256"] != generation_check_config_sha256:
        raise RuntimeError("range-check config changed during the numerical study")
    if (
        config_provenance["source_benchmark_config_sha256"]
        != generation_benchmark_config_sha256
    ):
        raise RuntimeError("source benchmark config changed during the range check")
    post_run_provenance = capture_reconstruction_provenance(
        repository_root,
        entry_points=(Path("scripts/run_reconstruction_curvature_range_check.py"),),
    )
    for field in ("source_files_sha256", "runtime_versions"):
        if post_run_provenance[field] != generation_provenance.get(field):
            raise RuntimeError(f"range-check {field} changed during the numerical study")
    configured_output = Path(str(run.check_config["output_directory"]))
    output_directory = (
        configured_output
        if configured_output.is_absolute()
        else repository_root / configured_output
    )
    output_directory.mkdir(parents=True, exist_ok=True)

    check_snapshot_path = output_directory / "range_check_config.json"
    benchmark_snapshot_path = output_directory / "source_benchmark_config.json"
    shutil.copyfile(check_config_path, check_snapshot_path)
    shutil.copyfile(benchmark_config_path, benchmark_snapshot_path)
    if file_sha256(check_snapshot_path) != config_provenance["check_config_sha256"]:
        raise RuntimeError("range-check config snapshot differs from its source")
    if (
        file_sha256(benchmark_snapshot_path)
        != config_provenance["source_benchmark_config_sha256"]
    ):
        raise RuntimeError("benchmark config snapshot differs from its source")

    candidate_path = output_directory / "calibration_range_check.csv"
    trial_path = output_directory / "calibration_trials.csv"
    summary_path = output_directory / "range_check_summary.json"
    metadata_path = output_directory / "metadata.json"
    write_rows(candidate_path, list(run.candidate_rows))
    write_rows(trial_path, list(run.trial_rows))

    reference_weight = float(run.check_config["reference_weight_um2"])
    rows_by_weight = {
        float(row["curvature_weight_um2"]): row for row in run.candidate_rows
    }
    reference_error = float(
        rows_by_weight[reference_weight][
            "median_supported_band_relative_l2_error"
        ]
    )
    stronger = {
        weight: float(row["median_supported_band_relative_l2_error"])
        for weight, row in rows_by_weight.items()
        if weight > reference_weight
    }
    all_stronger_worse = bool(
        stronger and all(error > reference_error for error in stronger.values())
    )
    write_json(
        summary_path,
        {
            "label": run.check_config["label"],
            "purpose": run.check_config["purpose"],
            "source_benchmark_label": run.source_benchmark_label,
            "readout": run.check_config["readout"],
            "basis_label": run.check_config["basis_label"],
            "calibration_seed": run.context.config["ensemble"]["calibration_seed"],
            "calibration_realizations_per_morphology": run.context.config["ensemble"][
                "calibration_realizations_per_morphology"
            ],
            "calibration_morphology_names": [
                case.name for case in run.context.morphology_split.calibration
            ],
            "reference_weight_um2": reference_weight,
            "reference_error": reference_error,
            "checked_results": list(run.candidate_rows),
            "conclusion": {
                "all_checked_stronger_weights_increase_median_error": all_stronger_worse,
                "scope": (
                    "calibration-set bracket for the frozen basis; not a new "
                    "model selection or held-out performance claim"
                ),
            },
            "study_elapsed_seconds": run.elapsed_seconds,
        },
    )
    output_hashes = {
        path.name: file_sha256(path)
        for path in (
            candidate_path,
            trial_path,
            summary_path,
            check_snapshot_path,
            benchmark_snapshot_path,
        )
    }
    write_json(
        metadata_path,
        {
            "label": run.check_config["label"],
            "entry_point": "scripts/run_reconstruction_curvature_range_check.py",
            "implementation": _metadata_path(Path(__file__), repository_root),
            "check_config": _metadata_path(check_config_path, repository_root),
            "check_config_sha256": config_provenance["check_config_sha256"],
            "source_benchmark_config": _metadata_path(
                benchmark_config_path, repository_root
            ),
            "source_benchmark_config_sha256": config_provenance[
                "source_benchmark_config_sha256"
            ],
            **generation_provenance,
            "study_elapsed_seconds": run.elapsed_seconds,
            "output_files_sha256": output_hashes,
        },
    )
    return {
        "candidate_summary": candidate_path,
        "calibration_trials": trial_path,
        "summary": summary_path,
        "check_config_snapshot": check_snapshot_path,
        "benchmark_config_snapshot": benchmark_snapshot_path,
        "metadata": metadata_path,
    }
