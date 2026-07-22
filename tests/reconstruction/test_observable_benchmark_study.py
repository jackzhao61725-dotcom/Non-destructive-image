from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from non_destructive_image.reconstruction.benchmark import DensityRecoveryMetrics
from non_destructive_image.reconstruction.studies.io import (
    file_sha256,
    load_json,
    load_rows,
    write_json,
    write_rows,
)
from non_destructive_image.reconstruction.studies.morphology import (
    build_morphology_study_context,
)
from non_destructive_image.reconstruction.studies.observable_benchmark import (
    GENERATION_PROVENANCE_FILENAME,
    LEGACY_VERIFICATION_FILENAME,
    OBSERVABLE_AGGREGATE_FILENAME,
    OBSERVABLE_MAPS_FILENAME,
    OBSERVABLE_TRIALS_FILENAME,
    SOURCE_CONFIG_SNAPSHOT_FILENAME,
    SOURCE_MANIFEST_SNAPSHOT_FILENAME,
    STUDY_CONFIG_SNAPSHOT_FILENAME,
    _seal_observable_artifact_set,
    build_observable_integration_support,
    validate_source_benchmark_artifacts,
    verify_legacy_trial_metrics,
)
from non_destructive_image.reconstruction.studies.observable_reporting import (
    generate_observable_benchmark_figures,
    verify_observable_artifact_set,
)
from scripts.generate_reconstruction_observable_benchmark import (
    DEFAULT_CONFIG as GENERATION_CONFIG,
)
from scripts.plot_reconstruction_observable_benchmark import (
    DEFAULT_CONFIG as PLOT_CONFIG,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
OBSERVABLE_CONFIG = (
    REPO_ROOT / "configs" / "reconstruction_observables_v1_orca_fusion_m10.json"
)
ACTIVE_OBSERVABLE_RESULTS = (
    REPO_ROOT / "results" / "reconstruction_observables_v1_orca_fusion_m10"
)
ACTIVE_OBSERVABLE_RUN_ID = (
    "81f7c045d0745d266dbcf73cb1e908a1c016aff12b2dc53df155b6815540635b"
)


def test_active_observable_source_is_sealed_complete_and_support_matched() -> None:
    config = load_json(OBSERVABLE_CONFIG)
    source = validate_source_benchmark_artifacts(config, REPO_ROOT)

    assert source.run_id == config["source_benchmark"]["expected_run_id"]
    assert len(source.held_out_rows) == 60
    assert {
        float(row["fluence_mw_us"]) for row in source.held_out_rows
    } == {30.0, 90.0, 150.0}
    assert {row["readout"] for row in source.held_out_rows} == {
        "dual_port",
        "dark_field",
    }

    context = build_morphology_study_context(source.config)
    support = build_observable_integration_support(config, source, context)
    assert support.shape == (306, 306)
    assert np.count_nonzero(support.support_mask) == 7425
    assert support.physical_area_m2 == pytest.approx(7.929642445213381e-10)


def test_observable_scripts_share_the_active_config() -> None:
    assert GENERATION_CONFIG == OBSERVABLE_CONFIG
    assert PLOT_CONFIG == OBSERVABLE_CONFIG


def _matching_legacy_source_row(metrics: DensityRecoveryMetrics) -> dict[str, str]:
    return {
        "readout": "dual_port",
        "fluence_mw_us": "90.0",
        "candidate": "frozen",
        "morphology": "held_out",
        "realization_index": "0",
        "seed": "123",
        "success": "True",
        "data_jacobian_rank": "2",
        "data_jacobian_condition": "3.0",
        "parameter_count": "2",
        "full_map_relative_l2_error": str(metrics.full_map_relative_l2_error),
        "supported_band_relative_l2_error": str(
            metrics.supported_band_relative_l2_error
        ),
        "integrated_density_relative_error": str(
            metrics.integrated_density_relative_error
        ),
        "centroid_y_error_um": str(metrics.centroid_y_error_um),
        "centroid_z_error_um": str(metrics.centroid_z_error_um),
        "rms_y_relative_error": str(metrics.rms_y_relative_error),
        "rms_z_relative_error": str(metrics.rms_z_relative_error),
        "message": "ok",
    }


def test_legacy_metric_closure_checks_every_stored_quantity() -> None:
    metrics = DensityRecoveryMetrics(0.1, 0.2, -0.03, 0.4, -0.5, 0.06, -0.07)
    fit = SimpleNamespace(
        coefficients=np.ones(2),
        diagnostics=SimpleNamespace(
            data_jacobian_rank=2,
            data_jacobian_condition=3.0,
        ),
    )
    source_row = _matching_legacy_source_row(metrics)
    rows = verify_legacy_trial_metrics(
        source_row,
        fit,
        metrics,
        relative_tolerance=1e-10,
        absolute_tolerance=1e-12,
    )
    assert len(rows) == 10
    assert all(row["within_tolerance"] for row in rows)

    source_row["rms_z_relative_error"] = "0.5"
    with pytest.raises(ValueError, match="rms_z_relative_error"):
        verify_legacy_trial_metrics(
            source_row,
            fit,
            metrics,
            relative_tolerance=1e-10,
            absolute_tolerance=1e-12,
        )


def _aggregate_row(readout: str, fluence: float) -> dict[str, object]:
    row: dict[str, object] = {
        "readout": readout,
        "fluence_mw_us": fluence,
        "trial_count": 1,
        "legacy_verified_trials": 1,
        "reconstructed_moment_supported_trials": 1,
        "principal_axis_angle_supported_pairs": 1,
    }
    stems = (
        "absolute_integrated_response_relative_error",
        "centroid_position_error_um",
        "covariance_tensor_relative_frobenius_error",
        "absolute_major_rms_width_relative_error",
        "absolute_minor_rms_width_relative_error",
        "absolute_aspect_ratio_relative_error",
        "absolute_principal_axis_angle_error_deg",
    )
    for index, stem in enumerate(stems, start=1):
        value = 0.01 * index + fluence * 1e-5
        row[f"supported_{stem}_trials"] = 1
        row[f"minimum_{stem}"] = 0.8 * value
        row[f"median_{stem}"] = value
        row[f"maximum_{stem}"] = 1.2 * value
    return row


def _write_fake_sealed_observable_result(directory: Path) -> None:
    source_run_id = "a" * 64
    config = {
        "label": "observable_test",
        "figure": {
            "basename": "observable_recovery_vs_fluence",
            "formats": ["png", "pdf", "svg"],
            "finite_ensemble_display": "finite range, not a confidence interval",
        },
    }
    write_json(directory / STUDY_CONFIG_SNAPSHOT_FILENAME, config)
    write_json(directory / SOURCE_CONFIG_SNAPSHOT_FILENAME, {"label": "source"})
    write_json(
        directory / SOURCE_MANIFEST_SNAPSHOT_FILENAME,
        {"schema_version": 1, "run_id": source_run_id},
    )
    write_json(
        directory / GENERATION_PROVENANCE_FILENAME,
        {"source_run_id": source_run_id},
    )
    write_rows(
        directory / OBSERVABLE_TRIALS_FILENAME,
        [{"readout": "dual_port", "fluence_mw_us": 30.0}],
    )
    write_rows(
        directory / LEGACY_VERIFICATION_FILENAME,
        [{"metric": "full_map_relative_l2_error", "within_tolerance": True}],
    )
    aggregate_rows = [
        _aggregate_row(readout, fluence)
        for readout in ("dual_port", "dark_field")
        for fluence in (30.0, 90.0, 150.0)
    ]
    write_rows(directory / OBSERVABLE_AGGREGATE_FILENAME, aggregate_rows)

    y_axis = np.linspace(-1e-6, 1e-6, 4)
    z_axis = np.linspace(-0.5e-6, 0.5e-6, 3)
    y_grid, z_grid = np.meshgrid(y_axis, z_axis)
    leading_shape = (2, 3, 1, 1)
    truth = np.ones((1, 3, 4))
    reconstructed = np.ones((*leading_shape, 3, 4))
    np.savez_compressed(
        directory / OBSERVABLE_MAPS_FILENAME,
        y_grid_m=y_grid,
        z_grid_m=z_grid,
        cell_area_m2=np.ones_like(y_grid),
        support_mask=np.ones_like(y_grid, dtype=bool),
        readout_names=np.asarray(["dual_port", "dark_field"]),
        fluence_mw_us=np.asarray([30.0, 90.0, 150.0]),
        morphology_names=np.asarray(["held_out"]),
        realization_indices=np.asarray([0]),
        candidate_labels=np.asarray(["frozen", "frozen"]),
        coefficient_count_by_readout=np.asarray([2, 2]),
        seeds=np.ones(leading_shape, dtype=np.int64),
        truth_column_density_m2=truth,
        reconstructed_column_density_m2=reconstructed,
        fitted_coefficients=np.ones((*leading_shape, 2)),
        source_run_id=np.asarray(source_run_id),
        tensor_axis_order=np.asarray(
            "readout,fluence,morphology,realization,z_index,y_index"
        ),
    )
    summary = {
        "label": "observable_test",
        "source_benchmark": {
            "run_id": source_run_id,
            "label": "source",
            "manifest_sha256": file_sha256(
                directory / SOURCE_MANIFEST_SNAPSHOT_FILENAME
            ),
            "config_sha256": file_sha256(
                directory / SOURCE_CONFIG_SNAPSHOT_FILENAME
            ),
        },
        "replay": {
            "all_legacy_metrics_verified": True,
            "readout_names": ["dual_port", "dark_field"],
            "fluence_mw_us": [30.0, 90.0, 150.0],
            "morphology_names": ["held_out"],
            "realization_indices": [0],
        },
        "integration_support": {
            "grid_shape": [3, 4],
            "angle_anisotropy_threshold": 0.05,
        },
        "maps_schema": {
            "truth_column_density_m2_shape": [1, 3, 4],
            "reconstructed_column_density_m2_shape": [2, 3, 1, 1, 3, 4],
        },
        "claims_boundary": {"synthetic_only": True},
    }
    _seal_observable_artifact_set(
        directory,
        label="observable_test",
        source_run_id=source_run_id,
        summary=summary,
    )


def test_sealed_observable_renderer_reads_no_live_fit_state(tmp_path: Path) -> None:
    _write_fake_sealed_observable_result(tmp_path)
    verified = verify_observable_artifact_set(tmp_path)
    outputs = generate_observable_benchmark_figures(tmp_path)

    assert verified.manifest["run_id"] == verified.summary["run_id"]
    assert outputs["figure"].is_file()
    assert outputs["metadata"].is_file()
    metadata = load_json(outputs["metadata"])
    assert metadata["source_run_id"] == verified.manifest["run_id"]
    assert set(metadata["output_files"]) == {
        "observable_recovery_vs_fluence.png",
        "observable_recovery_vs_fluence.pdf",
        "observable_recovery_vs_fluence.svg",
    }


def test_observable_renderer_rejects_tampered_aggregate(tmp_path: Path) -> None:
    _write_fake_sealed_observable_result(tmp_path)
    with (tmp_path / OBSERVABLE_AGGREGATE_FILENAME).open(
        "a",
        encoding="utf-8",
    ) as handle:
        handle.write("\n")
    with pytest.raises(ValueError, match="artifact hash mismatch"):
        generate_observable_benchmark_figures(tmp_path)


def test_active_observable_result_is_sealed_and_numerically_regressed() -> None:
    verified = verify_observable_artifact_set(ACTIVE_OBSERVABLE_RESULTS)
    assert verified.manifest["run_id"] == ACTIVE_OBSERVABLE_RUN_ID
    assert verified.summary["replay"]["trial_count"] == 60
    assert verified.summary["replay"]["all_legacy_metrics_verified"] is True

    rows = load_rows(verified.paths[OBSERVABLE_AGGREGATE_FILENAME])
    by_key = {
        (row["readout"], float(row["fluence_mw_us"])): row for row in rows
    }
    expected = {
        ("dual_port", 30.0): (
            0.05644627634330712,
            0.19532405394126215,
            0.1386587589263497,
            0.0661237214771,
            0.15457184495396437,
            0.09187294301290466,
            0.37767424254995535,
        ),
        ("dual_port", 90.0): (
            0.04144345212887979,
            0.1522251165861509,
            0.06466259656170287,
            0.03088574654162979,
            0.12711640051496265,
            0.09410062626413762,
            0.35465895777210987,
        ),
        ("dual_port", 150.0): (
            0.0441465893069759,
            0.11642117487956564,
            0.05196375289793523,
            0.02451563233112142,
            0.12427271987262956,
            0.09044909759857178,
            0.24213549065651668,
        ),
        ("dark_field", 30.0): (
            0.3193019858371664,
            1.4331271307966955,
            0.44673460703265044,
            0.20140631754819782,
            0.3395761866772363,
            0.12643500633341404,
            1.0621719876140023,
        ),
        ("dark_field", 90.0): (
            0.18351783647574582,
            0.47353660166626105,
            0.29081322149794697,
            0.1347705078263165,
            0.22095662477609523,
            0.07491259785130111,
            0.8949778504327086,
        ),
        ("dark_field", 150.0): (
            0.16154158203133773,
            0.2790659694591533,
            0.31927521455486596,
            0.14766156407615305,
            0.19144497373682356,
            0.048387295022380705,
            0.9667000319622416,
        ),
    }
    columns = (
        "median_absolute_integrated_response_relative_error",
        "median_centroid_position_error_um",
        "median_covariance_tensor_relative_frobenius_error",
        "median_absolute_major_rms_width_relative_error",
        "median_absolute_minor_rms_width_relative_error",
        "median_absolute_aspect_ratio_relative_error",
        "median_absolute_principal_axis_angle_error_deg",
    )
    assert set(by_key) == set(expected)
    for key, expected_values in expected.items():
        row = by_key[key]
        assert int(row["trial_count"]) == 10
        assert int(row["legacy_verified_trials"]) == 10
        assert int(row["reconstructed_moment_supported_trials"]) == 10
        assert int(row["principal_axis_angle_supported_pairs"]) == 10
        actual_values = tuple(float(row[column]) for column in columns)
        assert actual_values == pytest.approx(expected_values, rel=1e-12, abs=1e-12)

    figure_metadata = load_json(ACTIVE_OBSERVABLE_RESULTS / "figure_metadata.json")
    assert figure_metadata["source_run_id"] == ACTIVE_OBSERVABLE_RUN_ID
    for filename, expected_digest in figure_metadata[
        "rendered_output_sha256"
    ].items():
        assert file_sha256(ACTIVE_OBSERVABLE_RESULTS / filename) == expected_digest
