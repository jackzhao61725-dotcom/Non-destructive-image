from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from non_destructive_image.reconstruction.studies import artifacts as artifact_contract
from non_destructive_image.reconstruction.studies import (
    aggregate_held_out_trials,
    build_morphology_study_context,
    generate_morphology_benchmark_figures,
    file_sha256,
    load_json,
    make_study_measurement,
    seal_existing_morphology_benchmark_artifacts,
    write_json,
    write_rows,
)
from non_destructive_image.reconstruction.studies.morphology import _grid_convergence


REPO_ROOT = Path(__file__).resolve().parents[2]
STUDY_CONFIG = REPO_ROOT / "configs" / "reconstruction_morphology_benchmark_v3.json"
ORCA_STUDY_CONFIG = (
    REPO_ROOT
    / "configs"
    / "reconstruction_morphology_benchmark_v4_orca_fusion_m10.json"
)


def test_declared_study_context_resolves_candidates_and_fluence_scale() -> None:
    config = load_json(STUDY_CONFIG)
    assert "source_configs" not in config
    assert config["configuration_contract"]["inheritance"].startswith("none;")
    context = build_morphology_study_context(config)
    assert [candidate.model.parameter_count for candidate in context.candidates] == (
        [21] * 5 + [50] * 5 + [85] * 5
    )
    assert len({candidate.label for candidate in context.candidates}) == 15
    low = make_study_measurement(context, "dual_port", fluence_mw_us=30.0)
    high = make_study_measurement(context, "dual_port", fluence_mw_us=150.0)
    assert np.isclose(
        high.detector.photoelectrons_per_i0_pixel
        / low.detector.photoelectrons_per_i0_pixel,
        5.0,
    )


def test_orca_study_context_uses_physical_camera_contract() -> None:
    context = build_morphology_study_context(load_json(ORCA_STUDY_CONFIG))

    assert context.canonical_grid.sampling_mode == "physical_pixel"
    assert context.grid.sampling_mode == "physical_pixel"
    assert context.canonical_grid.camera_shape == context.grid.camera_shape == (153, 153)
    assert context.grid.y_grid_m.shape == (306, 306)
    assert context.reference_count_scale_e == pytest.approx(220.58087277528466)
    assert context.read_noise_e == pytest.approx(1.4)
    assert context.reduction.maximum_camera_coordinate_mismatch_m < 1e-15


def test_historical_orca_reconstruction_is_isolated_from_active_forward_config() -> None:
    reconstruction = load_json(ORCA_STUDY_CONFIG)
    forward = load_json(REPO_ROOT / "configs" / "dissertation_v3_orca_fusion.json")
    figure = load_json(REPO_ROOT / "configs" / "figure_5_1.json")

    # The sealed v4 benchmark remains method evidence under its original
    # optical and detector contract. Active Chapter 4/5 figures use a new
    # physical scenario and must not silently rewrite that provenance.
    assert reconstruction["physics"]["numerical_aperture"] == pytest.approx(0.08)
    assert reconstruction["physics"]["kappa_F"] == pytest.approx(1.0)
    assert reconstruction["detector"][
        "read_noise_electrons_per_pixel_per_readout"
    ] == pytest.approx(1.4)
    assert reconstruction["detector"]["reference_fluence_mw_us"] == pytest.approx(90.0)
    assert reconstruction["detector"][
        "photoelectrons_per_i0_pixel_at_reference"
    ] == pytest.approx(220.58087277528466)

    assert figure["fixed_context"]["numerical_aperture"] == pytest.approx(0.13)
    assert figure["fixed_context"]["kappa_F"] == pytest.approx(-45.0 / 91.0)
    assert forward["camera_recovery"]["read_noise_electrons"] == pytest.approx(0.7)
    assert figure["scan"]["reference_fluence_mw_us"] == pytest.approx(300.0)
    assert figure["canonical_checks"][
        "photoelectrons_per_incident_I0_pixel"
    ] == pytest.approx(735.2695759176155)

    assert reconstruction["physics"]["wavelength_m"] == forward["atom"][
        "transition_wavelength_m"
    ]
    assert reconstruction["grid"]["canonical_ngrid"] == forward["grid"]["ngrid"]
    assert reconstruction["grid"]["canonical_field_of_view_m"] == forward["grid"][
        "field_of_view_m"
    ]
    assert reconstruction["grid"]["camera_pixel_size_m"] == forward[
        "imaging_geometry"
    ]["object_plane_pixel_m"]
    assert reconstruction["grid"]["camera_output_shape"] == forward[
        "camera_recovery"
    ]["camera_output_shape"]
    assert reconstruction["detector"]["effective_magnification"] == forward[
        "imaging_geometry"
    ]["magnification"]


def test_orca_reduced_grid_passes_declared_signal_level_gate() -> None:
    context = build_morphology_study_context(load_json(ORCA_STUDY_CONFIG))
    rows = _grid_convergence(context)

    assert len(rows) == 10
    assert max(row["dual_port_signal_relative_l2_error"] for row in rows) < 0.01
    assert max(row["dark_field_relative_l2_error"] for row in rows) < 0.01
    assert max(abs(row["dark_field_peak_relative_error"]) for row in rows) < 0.01


def test_held_out_summary_keeps_failed_trials_in_denominator() -> None:
    rows = [
        {
            "readout": "dual_port",
            "fluence_mw_us": 90.0,
            "success": True,
            "supported_band_relative_l2_error": 0.2,
            "integrated_density_relative_error": -0.1,
        },
        {
            "readout": "dual_port",
            "fluence_mw_us": 90.0,
            "success": False,
            "supported_band_relative_l2_error": "",
            "integrated_density_relative_error": "",
        },
    ]
    summary = aggregate_held_out_trials(rows)
    assert summary == (
        {
            "readout": "dual_port",
            "fluence_mw_us": 90.0,
            "trial_count": 2,
            "successful_trials": 1,
            "success_fraction": 0.5,
            "median_supported_band_relative_l2_error": 0.2,
            "minimum_supported_band_relative_l2_error": 0.2,
            "maximum_supported_band_relative_l2_error": 0.2,
            "median_absolute_integrated_density_relative_error": 0.1,
        },
    )


def _write_fake_sealed_benchmark(tmp_path: Path) -> None:
    summary_rows = []
    for readout, offset in (("dual_port", 0.0), ("dark_field", 0.1)):
        for fluence, error in ((30.0, 0.4), (90.0, 0.2), (150.0, 0.15)):
            summary_rows.append(
                {
                    "readout": readout,
                    "fluence_mw_us": fluence,
                    "trial_count": 2,
                    "successful_trials": 2,
                    "success_fraction": 1.0,
                    "median_supported_band_relative_l2_error": error + offset,
                    "minimum_supported_band_relative_l2_error": error + offset,
                    "maximum_supported_band_relative_l2_error": error + offset,
                    "median_absolute_integrated_density_relative_error": 0.05,
                }
            )
    write_rows(tmp_path / "held_out_summary.csv", summary_rows)
    write_rows(tmp_path / "calibration_candidate_summary.csv", [{"candidate": "a"}])
    write_rows(tmp_path / "held_out_trials.csv", [{"trial": 0, "success": True}])
    write_rows(tmp_path / "reduced_grid_convergence.csv", [{"passed": True}])
    write_json(
        tmp_path / "benchmark_summary.json",
        {
            "label": "synthetic_test_benchmark",
            "representative": {
                "artifact": "representative_reconstruction.npz",
                "morphology_name": "synthetic_test",
                "fluence_mw_us": 90.0,
                "realization_index": 0,
                "seeds": {"dual_port": 11, "dark_field": 12},
                "supported_band_errors": {"dual_port": 0.2, "dark_field": 0.3},
            },
            "interpretation": {"faraday_absolute_scale": "synthetic"},
        },
    )
    y_axis = np.linspace(-30e-6, 30e-6, 9)
    z_axis = np.linspace(-8e-6, 8e-6, 5)
    y_grid, z_grid = np.meshgrid(y_axis, z_axis)
    truth = 5e14 * np.exp(-((y_grid / 15e-6) ** 2 + (z_grid / 4e-6) ** 2))
    np.savez_compressed(
        tmp_path / "representative_reconstruction.npz",
        y_grid_m=y_grid,
        z_grid_m=z_grid,
        truth_column_density_m2=truth,
        dual_port_column_density_m2=0.95 * truth,
        dark_field_column_density_m2=0.85 * truth,
        morphology_name=np.asarray("synthetic_test"),
        fluence_mw_us=np.asarray(90.0),
        realization_index=np.asarray(0),
        dual_port_seed=np.asarray(11),
        dark_field_seed=np.asarray(12),
        dual_port_supported_band_error=np.asarray(0.2),
        dark_field_supported_band_error=np.asarray(0.3),
    )
    config_path = tmp_path / "input_config.json"
    write_json(config_path, {"label": "synthetic_test_benchmark"})
    write_json(
        tmp_path / "metadata.json",
        {
            "label": "synthetic_test_benchmark",
            "config_sha256": file_sha256(config_path),
            "output_files": [],
        },
    )
    seal_existing_morphology_benchmark_artifacts(tmp_path, config_path)


def test_figure_generation_reads_one_sealed_artifact_set(tmp_path: Path) -> None:
    _write_fake_sealed_benchmark(tmp_path)
    outputs = generate_morphology_benchmark_figures(tmp_path)
    assert outputs["quality_curve"].is_file()
    assert outputs["representative_map"].is_file()
    assert outputs["metadata"].is_file()
    figure_metadata = load_json(outputs["metadata"])
    manifest = load_json(tmp_path / "artifact_manifest.json")
    assert figure_metadata["source_run_id"] == manifest["run_id"]
    assert figure_metadata["source_artifact_sha256"] == {
        name: record["sha256"] for name, record in manifest["artifacts"].items()
    }


def test_figure_generation_rejects_tampered_artifact(tmp_path: Path) -> None:
    _write_fake_sealed_benchmark(tmp_path)
    with (tmp_path / "held_out_summary.csv").open("a", encoding="utf-8") as handle:
        handle.write("\n")
    with pytest.raises(ValueError, match="artifact hash mismatch"):
        generate_morphology_benchmark_figures(tmp_path)


def test_figure_generation_requires_manifest(tmp_path: Path) -> None:
    _write_fake_sealed_benchmark(tmp_path)
    (tmp_path / "artifact_manifest.json").unlink()
    with pytest.raises(FileNotFoundError, match="artifact_manifest.json"):
        generate_morphology_benchmark_figures(tmp_path)


def test_figure_generation_rejects_cross_file_representative_mismatch(
    tmp_path: Path,
) -> None:
    _write_fake_sealed_benchmark(tmp_path)
    representative_path = tmp_path / "representative_reconstruction.npz"
    with np.load(representative_path, allow_pickle=False) as arrays:
        payload = {key: np.array(arrays[key], copy=True) for key in arrays.files}
    payload["dual_port_seed"] = np.asarray(999)
    np.savez_compressed(representative_path, **payload)

    manifest_path = tmp_path / "artifact_manifest.json"
    summary_path = tmp_path / "benchmark_summary.json"
    manifest = load_json(manifest_path)
    manifest["artifacts"][representative_path.name]["sha256"] = file_sha256(
        representative_path
    )
    raw_hashes = {
        filename: manifest["artifacts"][filename]["sha256"]
        for filename in artifact_contract.RAW_ARTIFACT_ROLES
    }
    run_id = artifact_contract._deterministic_run_id(
        manifest["config"]["sha256"], raw_hashes
    )
    summary = load_json(summary_path)
    summary["run_id"] = run_id
    write_json(summary_path, summary)
    manifest["run_id"] = run_id
    manifest["artifacts"][summary_path.name]["sha256"] = file_sha256(summary_path)
    write_json(manifest_path, manifest)

    with pytest.raises(ValueError, match="dual_port seed differs"):
        generate_morphology_benchmark_figures(tmp_path)
