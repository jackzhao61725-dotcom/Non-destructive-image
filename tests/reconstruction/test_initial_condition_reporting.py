from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from non_destructive_image.reconstruction.studies.initial_condition_reporting import (
    _shared_density_limit,
    generate_initial_condition_suite_figures,
    verify_initial_condition_artifact_set,
)
from non_destructive_image.reconstruction.studies.initial_condition_suite import (
    CONDITION_SUMMARY_FILENAME,
    MANIFEST_FILENAME,
    MAPS_FILENAME,
    PROVENANCE_FILENAME,
    RAW_ARTIFACT_ROLES,
    SOURCE_CONFIG_FILENAME,
    SOURCE_MANIFEST_FILENAME,
    STUDY_CONFIG_FILENAME,
    SUMMARY_FILENAME,
    TRIALS_FILENAME,
    _deterministic_run_id,
)
from non_destructive_image.reconstruction.studies.io import (
    file_sha256,
    load_json,
    write_json,
    write_rows,
)
from scripts.plot_dpfi_initial_condition_suite import DEFAULT_CONFIG as PLOT_CONFIG
from scripts.run_dpfi_initial_condition_suite import DEFAULT_CONFIG as RUN_CONFIG


def _write_fake_sealed_suite(directory: Path) -> None:
    label = "fake_dpfi_initial_condition_suite"
    source_run_id = "a" * 64
    condition_ids = ("condition_alpha", "condition_beta")
    fluences = (90.0, 150.0, 300.0)
    realizations = (0,)
    config = {
        "label": label,
        "output_directory": str(directory),
        "initial_conditions": [{"id": identifier} for identifier in condition_ids],
        "ensemble": {
            "fluence_mw_us": list(fluences),
            "realizations_per_condition": 1,
        },
        "figure": {
            "representative_fluence_mw_us": 90.0,
            "formats": ["png"],
            "overview_basename": "condition_overview_F90",
            "metric_basename": "recovery_metrics_by_condition",
        },
    }
    write_json(directory / STUDY_CONFIG_FILENAME, config)
    write_json(directory / SOURCE_CONFIG_FILENAME, {"label": "source"})
    write_json(
        directory / SOURCE_MANIFEST_FILENAME,
        {"schema_version": 1, "run_id": source_run_id},
    )
    write_json(directory / PROVENANCE_FILENAME, {"source_run_id": source_run_id})

    leading_shape = (len(fluences), len(condition_ids), len(realizations))
    trial_rows: list[dict[str, object]] = []
    seeds = np.empty(leading_shape, dtype=np.int64)
    fit_success = np.ones(leading_shape, dtype=bool)
    for fluence_index, fluence in enumerate(fluences):
        for condition_index, identifier in enumerate(condition_ids):
            seed = 1000 + 100 * fluence_index + condition_index
            seeds[fluence_index, condition_index, 0] = seed
            trial_rows.append(
                {
                    "condition_id": identifier,
                    "fluence_mw_us": fluence,
                    "realization_index": 0,
                    "seed": seed,
                    "fit_success": True,
                    "c_A": 0.25 + 0.1 * fluence_index + 0.05 * condition_index,
                    "c_r": 0.35 + 0.1 * fluence_index + 0.05 * condition_index,
                    "c_w": 0.45 + 0.1 * fluence_index + 0.05 * condition_index,
                }
            )
    write_rows(directory / TRIALS_FILENAME, trial_rows)
    write_rows(
        directory / CONDITION_SUMMARY_FILENAME,
        [
            {
                "condition_id": identifier,
                "trial_count": len(fluences),
                "successful_fit_count": len(fluences),
            }
            for identifier in condition_ids
        ],
    )

    y_axis = np.linspace(-3e-6, 3e-6, 6)
    z_axis = np.linspace(-2e-6, 2e-6, 4)
    y_grid, z_grid = np.meshgrid(y_axis, z_axis)
    support = np.ones_like(y_grid, dtype=bool)
    truth = np.stack(
        [
            (index + 1.0) * 1e14 * np.exp(-((y_grid / 2e-6) ** 2 + (z_grid / 1e-6) ** 2))
            for index in range(len(condition_ids))
        ]
    )
    reconstructed = np.empty((*leading_shape, *y_grid.shape), dtype=float)
    for fluence_index in range(len(fluences)):
        reconstructed[fluence_index, :, 0] = truth * (0.9 + 0.05 * fluence_index)
    np.savez_compressed(
        directory / MAPS_FILENAME,
        y_grid_m=y_grid,
        z_grid_m=z_grid,
        cell_area_m2=np.full_like(y_grid, 1e-12),
        support_mask=support,
        fluence_mw_us=np.asarray(fluences),
        condition_ids=np.asarray(condition_ids),
        realization_indices=np.asarray(realizations),
        candidate_label=np.asarray("frozen_candidate"),
        seeds=seeds,
        fit_success=fit_success,
        truth_column_density_m2=truth,
        reconstructed_column_density_m2=reconstructed,
        fitted_coefficients=np.ones((*leading_shape, 2)),
        dual_port_h_counts_e=np.ones((*leading_shape, 2, 3)),
        dual_port_v_counts_e=np.ones((*leading_shape, 2, 3)),
        source_run_id=np.asarray(source_run_id),
        tensor_axis_order=np.asarray(
            "fluence,condition,realization,z_index,y_index"
        ),
    )

    raw_hashes = {
        filename: file_sha256(directory / filename) for filename in RAW_ARTIFACT_ROLES
    }
    config_sha256 = raw_hashes[STUDY_CONFIG_FILENAME]
    run_id = _deterministic_run_id(config_sha256, source_run_id, raw_hashes)
    summary = {
        "label": label,
        "run_id": run_id,
        "source_benchmark": {
            "run_id": source_run_id,
            "manifest_sha256": file_sha256(directory / SOURCE_MANIFEST_FILENAME),
            "config_sha256": file_sha256(directory / SOURCE_CONFIG_FILENAME),
        },
        "suite": {
            "readout": "dual_port",
            "trial_count": len(trial_rows),
            "condition_ids": list(condition_ids),
            "fluence_mw_us": list(fluences),
            "realization_indices": list(realizations),
            "successful_fit_count": len(trial_rows),
            "selected_candidate": "frozen_candidate",
            "parameter_count": 2,
        },
        "observable_usability": {
            "scores_are_reported_independently": True,
            "overall_image_usability_score": None,
        },
        "integration_support": {
            "grid_shape": list(y_grid.shape),
            "supported_cell_count": int(np.count_nonzero(support)),
        },
        "maps_schema": {
            "reconstruction_axis_order": [
                "fluence",
                "condition",
                "realization",
                "z_index",
                "y_index",
            ]
        },
        "claims_boundary": {"synthetic_only": True},
    }
    write_json(directory / SUMMARY_FILENAME, summary)
    artifacts = {
        filename: {"role": RAW_ARTIFACT_ROLES[filename], "sha256": digest}
        for filename, digest in raw_hashes.items()
    }
    artifacts[SUMMARY_FILENAME] = {
        "role": "run_summary",
        "sha256": file_sha256(directory / SUMMARY_FILENAME),
    }
    write_json(
        directory / MANIFEST_FILENAME,
        {
            "schema_version": 1,
            "label": label,
            "run_id": run_id,
            "source_run_id": source_run_id,
            "config": {
                "artifact": STUDY_CONFIG_FILENAME,
                "sha256": config_sha256,
            },
            "artifacts": artifacts,
        },
    )


def test_sealed_suite_renderer_verifies_then_hashes_derived_figures(
    tmp_path: Path,
) -> None:
    _write_fake_sealed_suite(tmp_path)
    source_manifest_hash = file_sha256(tmp_path / MANIFEST_FILENAME)

    verified = verify_initial_condition_artifact_set(tmp_path)
    outputs = generate_initial_condition_suite_figures(tmp_path)

    assert verified.manifest["run_id"] == verified.summary["run_id"]
    assert set(outputs) == {"overview", "metrics", "metadata"}
    assert all(path.is_file() for path in outputs.values())
    metadata = load_json(outputs["metadata"])
    assert metadata["source_manifest_sha256"] == source_manifest_hash
    assert metadata["metric_coefficients"] == ["c_A", "c_r", "c_w"]
    assert metadata["metric_threshold"] == 1.0
    assert "overall_image_usability_score" not in metadata
    assert set(metadata["output_files"]) == {
        "condition_overview_F90.png",
        "recovery_metrics_by_condition.png",
    }
    for filename, expected_hash in metadata["rendered_output_sha256"].items():
        assert file_sha256(tmp_path / filename) == expected_hash
    assert file_sha256(tmp_path / MANIFEST_FILENAME) == source_manifest_hash
    assert "figure_metadata.json" not in verified.manifest["artifacts"]


def test_suite_renderer_rejects_tampering_before_plotting(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _write_fake_sealed_suite(tmp_path)
    with (tmp_path / TRIALS_FILENAME).open("a", encoding="utf-8") as handle:
        handle.write("\n")

    def plotting_would_be_a_failure(*args: object, **kwargs: object) -> None:
        raise AssertionError("plotting began before artifact verification")

    monkeypatch.setattr(
        "non_destructive_image.reconstruction.studies.initial_condition_reporting."
        "_plot_condition_overview",
        plotting_would_be_a_failure,
    )
    with pytest.raises(ValueError, match="artifact hash mismatch"):
        generate_initial_condition_suite_figures(tmp_path)


def test_truth_and_reconstruction_use_one_condition_scale() -> None:
    assert _shared_density_limit(
        np.asarray([[1.0, 5.0]]),
        np.asarray([[2.0, 7.0]]),
    ) == 7.0


def test_initial_condition_scripts_share_the_active_config() -> None:
    assert PLOT_CONFIG == RUN_CONFIG
