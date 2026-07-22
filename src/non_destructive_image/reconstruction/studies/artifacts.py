"""Serialize morphology benchmark results as a sealed artifact set.

The benchmark workflow returns a :class:`MorphologyBenchmarkRun` in memory.  This
module is the only layer that knows how that result is laid out on disk.  Every
completed run is sealed by a manifest so a later reporting process can prove that
all tables, arrays, summary data and the frozen configuration belong together.
"""

from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import shutil
from typing import Any

import numpy as np

from .io import file_sha256, load_json, write_json, write_rows
from .morphology import MorphologyBenchmarkRun
from .provenance import capture_reconstruction_provenance


REPRESENTATIVE_ARRAYS_FILENAME = "representative_reconstruction.npz"
CONFIG_SNAPSHOT_FILENAME = "study_config.json"
MANIFEST_FILENAME = "artifact_manifest.json"
MANIFEST_SCHEMA_VERSION = 1

RAW_ARTIFACT_ROLES = {
    "calibration_candidate_summary.csv": "calibration_candidates",
    "held_out_trials.csv": "held_out_trials",
    "held_out_summary.csv": "held_out_summary",
    "reduced_grid_convergence.csv": "grid_convergence",
    REPRESENTATIVE_ARRAYS_FILENAME: "representative_arrays",
}
SUMMARY_FILENAME = "benchmark_summary.json"


def _path_for_metadata(path: Path, repository_root: Path) -> str:
    """Return a stable repository-relative path when one is available."""

    try:
        return str(path.resolve().relative_to(repository_root.resolve()))
    except ValueError:
        return str(path.resolve())


def _selected_candidates(run: MorphologyBenchmarkRun) -> dict[str, dict[str, Any]]:
    selected: dict[str, dict[str, Any]] = {}
    for readout, choice in run.selected.items():
        regularisation = choice.candidate.regularisation
        selected[readout] = {
            "label": choice.candidate.label,
            "parameter_count": choice.candidate.model.parameter_count,
            "curvature_weight_um2": (
                0.0 if regularisation is None else regularisation.weight_um2
            ),
            "regularisation_density_scale_m2": (
                None if regularisation is None else regularisation.density_scale_m2
            ),
            "regularisation_boundary_policy": (
                None if regularisation is None else regularisation.boundary_policy
            ),
            "calibration_morphology_names": choice.calibration_morphology_names,
        }
    return selected


def _representative_summary(run: MorphologyBenchmarkRun) -> dict[str, Any]:
    representative = run.representative
    return {
        "artifact": REPRESENTATIVE_ARRAYS_FILENAME,
        "morphology_name": representative.morphology_name,
        "fluence_mw_us": representative.fluence_mw_us,
        "realization_index": representative.realization_index,
        "seeds": representative.seeds,
        "supported_band_errors": representative.supported_band_errors,
    }


def _write_representative_arrays(path: Path, run: MorphologyBenchmarkRun) -> None:
    representative = run.representative
    np.savez_compressed(
        path,
        y_grid_m=run.context.grid.y_grid_m,
        z_grid_m=run.context.grid.z_grid_m,
        truth_column_density_m2=representative.truth_column_density_m2,
        dual_port_column_density_m2=representative.dual_port_column_density_m2,
        dark_field_column_density_m2=representative.dark_field_column_density_m2,
        dual_port_h_counts_e=representative.dual_port_h_counts_e,
        dual_port_v_counts_e=representative.dual_port_v_counts_e,
        dark_field_counts_e=representative.dark_field_counts_e,
        morphology_name=np.asarray(representative.morphology_name),
        fluence_mw_us=np.asarray(representative.fluence_mw_us),
        realization_index=np.asarray(representative.realization_index),
        dual_port_seed=np.asarray(representative.seeds["dual_port"]),
        dark_field_seed=np.asarray(representative.seeds["dark_field"]),
        dual_port_supported_band_error=np.asarray(
            representative.supported_band_errors["dual_port"]
        ),
        dark_field_supported_band_error=np.asarray(
            representative.supported_band_errors["dark_field"]
        ),
    )


def _artifact_hashes(output_directory: Path) -> dict[str, str]:
    """Hash the raw artifacts used to identify one numerical run."""

    hashes: dict[str, str] = {}
    for filename in RAW_ARTIFACT_ROLES:
        path = output_directory / filename
        if not path.is_file():
            raise FileNotFoundError(f"required raw artifact does not exist: {path}")
        hashes[filename] = file_sha256(path)
    return hashes


def _deterministic_run_id(
    config_sha256: str,
    raw_artifact_hashes: dict[str, str],
) -> str:
    """Identify a run from its frozen configuration and numerical data."""

    payload = json.dumps(
        {
            "config_sha256": config_sha256,
            "raw_artifact_sha256": raw_artifact_hashes,
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _manifest_payload(
    *,
    label: str,
    run_id: str,
    config_sha256: str,
    output_directory: Path,
    raw_hashes: dict[str, str],
) -> dict[str, Any]:
    artifact_records: dict[str, dict[str, str]] = {
        CONFIG_SNAPSHOT_FILENAME: {
            "role": "configuration_snapshot",
            "sha256": config_sha256,
        }
    }
    for filename, digest in raw_hashes.items():
        artifact_records[filename] = {
            "role": RAW_ARTIFACT_ROLES[filename],
            "sha256": digest,
        }
    artifact_records[SUMMARY_FILENAME] = {
        "role": "run_summary",
        "sha256": file_sha256(output_directory / SUMMARY_FILENAME),
    }
    return {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "label": label,
        "run_id": run_id,
        "config": {
            "artifact": CONFIG_SNAPSHOT_FILENAME,
            "sha256": config_sha256,
        },
        "artifacts": artifact_records,
    }


def _seal_artifact_set(
    *,
    output_directory: Path,
    label: str,
    config_path: Path,
    summary: dict[str, Any],
) -> tuple[str, Path, Path]:
    """Freeze the config, add a run ID to the summary, and write a manifest."""

    config_sha256 = file_sha256(config_path)
    config_snapshot_path = output_directory / CONFIG_SNAPSHOT_FILENAME
    shutil.copyfile(config_path, config_snapshot_path)
    if file_sha256(config_snapshot_path) != config_sha256:
        raise RuntimeError("configuration snapshot differs from the source config")

    raw_hashes = _artifact_hashes(output_directory)
    run_id = _deterministic_run_id(config_sha256, raw_hashes)
    summary = dict(summary)
    summary["run_id"] = run_id
    summary_path = output_directory / SUMMARY_FILENAME
    write_json(summary_path, summary)

    manifest_path = output_directory / MANIFEST_FILENAME
    write_json(
        manifest_path,
        _manifest_payload(
            label=label,
            run_id=run_id,
            config_sha256=config_sha256,
            output_directory=output_directory,
            raw_hashes=raw_hashes,
        ),
    )
    return run_id, manifest_path, config_snapshot_path


def _npz_scalar(arrays: Any, key: str) -> Any:
    if key not in arrays.files:
        raise ValueError(f"representative NPZ is missing metadata scalar: {key}")
    value = np.asarray(arrays[key])
    if value.ndim != 0:
        raise ValueError(f"representative NPZ metadata is not scalar: {key}")
    return value.item()


def _validate_existing_representative_summary(
    output_directory: Path,
    summary: dict[str, Any],
) -> None:
    representative = summary.get("representative")
    if not isinstance(representative, dict):
        raise ValueError("existing summary has no representative record")
    if representative.get("artifact") != REPRESENTATIVE_ARRAYS_FILENAME:
        raise ValueError("existing summary does not name the representative NPZ")
    seeds = representative.get("seeds")
    errors = representative.get("supported_band_errors")
    if not isinstance(seeds, dict) or not isinstance(errors, dict):
        raise ValueError("existing representative seeds and errors must be mappings")

    path = output_directory / REPRESENTATIVE_ARRAYS_FILENAME
    with np.load(path, allow_pickle=False) as arrays:
        if str(_npz_scalar(arrays, "morphology_name")) != str(
            representative.get("morphology_name")
        ):
            raise ValueError("representative morphology differs between summary and NPZ")
        if int(_npz_scalar(arrays, "realization_index")) != int(
            representative.get("realization_index")
        ):
            raise ValueError(
                "representative realization index differs between summary and NPZ"
            )
        if not np.isclose(
            float(_npz_scalar(arrays, "fluence_mw_us")),
            float(representative.get("fluence_mw_us")),
            rtol=1e-12,
            atol=1e-12,
        ):
            raise ValueError("representative fluence differs between summary and NPZ")
        for readout in ("dual_port", "dark_field"):
            if readout not in seeds or readout not in errors:
                raise ValueError(f"existing summary lacks {readout} seed or error")
            if int(_npz_scalar(arrays, f"{readout}_seed")) != int(seeds[readout]):
                raise ValueError(
                    f"representative {readout} seed differs between summary and NPZ"
                )
            if not np.isclose(
                float(_npz_scalar(arrays, f"{readout}_supported_band_error")),
                float(errors[readout]),
                rtol=1e-12,
                atol=1e-12,
            ):
                raise ValueError(
                    f"representative {readout} error differs between summary and NPZ"
                )


def seal_existing_morphology_benchmark_artifacts(
    output_directory: Path,
    config_path: Path,
) -> dict[str, Path]:
    """Seal an earlier complete run without repeating any reconstruction fit.

    This migration helper refuses to proceed unless the pre-existing metadata
    records the same configuration hash.  It only adds a frozen config snapshot,
    a deterministic run ID, and the manifest; numerical artifacts are not edited.
    """

    output_directory = output_directory.resolve()
    config_path = config_path.resolve()
    metadata_path = output_directory / "metadata.json"
    summary_path = output_directory / SUMMARY_FILENAME
    if not config_path.is_file():
        raise FileNotFoundError(config_path)
    if not metadata_path.is_file():
        raise FileNotFoundError(metadata_path)
    if not summary_path.is_file():
        raise FileNotFoundError(summary_path)

    metadata = load_json(metadata_path)
    config_sha256 = file_sha256(config_path)
    if metadata.get("config_sha256") != config_sha256:
        raise ValueError(
            "existing metadata config hash does not match the supplied config"
        )
    summary = load_json(summary_path)
    label = str(summary.get("label") or metadata.get("label") or "")
    if not label:
        raise ValueError("existing artifacts do not declare a benchmark label")
    _artifact_hashes(output_directory)
    _validate_existing_representative_summary(output_directory, summary)

    run_id, manifest_path, config_snapshot_path = _seal_artifact_set(
        output_directory=output_directory,
        label=label,
        config_path=config_path,
        summary=summary,
    )
    metadata["run_id"] = run_id
    metadata["artifact_manifest"] = manifest_path.name
    output_files = list(metadata.get("output_files") or [])
    for filename in (config_snapshot_path.name, manifest_path.name):
        if filename not in output_files:
            output_files.append(filename)
    metadata["output_files"] = output_files
    write_json(metadata_path, metadata)
    return {
        "summary": summary_path,
        "config_snapshot": config_snapshot_path,
        "manifest": manifest_path,
        "metadata": metadata_path,
    }


def write_morphology_benchmark_run(
    run: MorphologyBenchmarkRun,
    config_path: Path,
    repository_root: Path,
    *,
    generation_config_sha256: str,
    generation_provenance: dict[str, object],
) -> dict[str, Path]:
    """Write one completed benchmark run and return its artifact paths.

    The output directory is taken from ``run.context.config['output_directory']``
    relative to ``repository_root``.  ``config_path`` is frozen and hashed for
    provenance; it must already exist.  No model calculation is performed here.
    """

    repository_root = repository_root.resolve()
    config_path = config_path.resolve()
    if not config_path.is_file():
        raise FileNotFoundError(config_path)
    if file_sha256(config_path) != generation_config_sha256:
        raise RuntimeError("benchmark config changed during the numerical study")
    if load_json(config_path) != run.context.config:
        raise ValueError("benchmark run context does not match the supplied config")
    post_run_provenance = capture_reconstruction_provenance(
        repository_root,
        entry_points=(Path("scripts/generate_reconstruction_morphology_benchmark.py"),),
    )
    for field in ("source_files_sha256", "runtime_versions"):
        if post_run_provenance[field] != generation_provenance.get(field):
            raise RuntimeError(f"benchmark {field} changed during the numerical study")

    configured_output = Path(str(run.context.config["output_directory"]))
    output_directory = (
        configured_output
        if configured_output.is_absolute()
        else repository_root / configured_output
    )
    output_directory.mkdir(parents=True, exist_ok=True)

    candidate_path = output_directory / "calibration_candidate_summary.csv"
    held_out_path = output_directory / "held_out_trials.csv"
    held_out_summary_path = output_directory / "held_out_summary.csv"
    convergence_path = output_directory / "reduced_grid_convergence.csv"
    summary_path = output_directory / SUMMARY_FILENAME
    representative_path = output_directory / REPRESENTATIVE_ARRAYS_FILENAME
    metadata_path = output_directory / "metadata.json"

    write_rows(candidate_path, list(run.candidate_rows))
    write_rows(held_out_path, list(run.held_out_rows))
    write_rows(held_out_summary_path, list(run.held_out_summary_rows))
    write_rows(convergence_path, list(run.convergence_rows))
    _write_representative_arrays(representative_path, run)

    config = run.context.config
    physics = config["physics"]
    ensemble = config["ensemble"]
    summary = {
        "label": config["label"],
        "grid_reduction": asdict(run.context.reduction),
        "grid_convergence_gate_passed": True,
        "study_elapsed_seconds": run.elapsed_seconds,
        "selected_candidates": _selected_candidates(run),
        "held_out_summaries": run.held_out_summaries,
        "representative": _representative_summary(run),
        "interpretation": {
            "status": "synthetic calibration and held-out assessment only",
            "truth_blind_fit": True,
            "faraday_absolute_scale": physics["faraday_scale_status"],
            "publication_uncertainty_boundary": ensemble["status"],
        },
    }
    run_id, manifest_path, config_snapshot_path = _seal_artifact_set(
        output_directory=output_directory,
        label=str(config["label"]),
        config_path=config_path,
        summary=summary,
    )

    output_files = [
        candidate_path.name,
        held_out_path.name,
        held_out_summary_path.name,
        convergence_path.name,
        summary_path.name,
        representative_path.name,
        config_snapshot_path.name,
        manifest_path.name,
    ]
    write_json(
        metadata_path,
        {
            "label": config["label"],
            "run_id": run_id,
            "artifact_manifest": manifest_path.name,
            "entry_point": "scripts/generate_reconstruction_morphology_benchmark.py",
            "writer": _path_for_metadata(Path(__file__), repository_root),
            "config": _path_for_metadata(config_path, repository_root),
            "config_sha256": file_sha256(config_path),
            **generation_provenance,
            "study_elapsed_seconds": run.elapsed_seconds,
            "output_files": output_files,
        },
    )
    return {
        "candidate_summary": candidate_path,
        "held_out_trials": held_out_path,
        "held_out_summary": held_out_summary_path,
        "grid_convergence": convergence_path,
        "summary": summary_path,
        "representative_arrays": representative_path,
        "config_snapshot": config_snapshot_path,
        "manifest": manifest_path,
        "metadata": metadata_path,
    }
