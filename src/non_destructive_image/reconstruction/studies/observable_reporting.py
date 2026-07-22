"""Render Chapter 7 observable recovery evidence from sealed artifacts only."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from .io import file_sha256, load_json, load_rows, write_json
from .observable_benchmark import (
    OBSERVABLE_AGGREGATE_FILENAME,
    OBSERVABLE_MANIFEST_FILENAME,
    OBSERVABLE_MANIFEST_SCHEMA_VERSION,
    OBSERVABLE_MAPS_FILENAME,
    OBSERVABLE_RAW_ARTIFACT_ROLES,
    OBSERVABLE_SUMMARY_FILENAME,
    SOURCE_CONFIG_SNAPSHOT_FILENAME,
    SOURCE_MANIFEST_SNAPSHOT_FILENAME,
    STUDY_CONFIG_SNAPSHOT_FILENAME,
    _deterministic_observable_run_id,
    _safe_artifact_path,
)


FIGURE_METADATA_FILENAME = "figure_metadata.json"


@dataclass(frozen=True)
class VerifiedObservableArtifactSet:
    """Paths and metadata proven to belong to one observable numerical run."""

    directory: Path
    manifest: dict[str, Any] = field(repr=False)
    summary: dict[str, Any] = field(repr=False)
    config: dict[str, Any] = field(repr=False)
    paths: dict[str, Path] = field(repr=False)


def _validate_maps_schema(path: Path, summary: dict[str, Any]) -> None:
    required = {
        "y_grid_m",
        "z_grid_m",
        "cell_area_m2",
        "support_mask",
        "readout_names",
        "fluence_mw_us",
        "morphology_names",
        "realization_indices",
        "candidate_labels",
        "coefficient_count_by_readout",
        "seeds",
        "truth_column_density_m2",
        "reconstructed_column_density_m2",
        "fitted_coefficients",
        "source_run_id",
        "tensor_axis_order",
    }
    with np.load(path, allow_pickle=False) as arrays:
        missing = required.difference(arrays.files)
        if missing:
            raise ValueError(
                "observable maps NPZ is missing keys: " + ", ".join(sorted(missing))
            )
        grid_shape = tuple(int(value) for value in summary["integration_support"]["grid_shape"])
        for key in ("y_grid_m", "z_grid_m", "cell_area_m2", "support_mask"):
            if arrays[key].shape != grid_shape:
                raise ValueError(f"observable maps {key} shape disagrees with summary")
        truth_shape = tuple(
            int(value)
            for value in summary["maps_schema"]["truth_column_density_m2_shape"]
        )
        reconstructed_shape = tuple(
            int(value)
            for value in summary["maps_schema"][
                "reconstructed_column_density_m2_shape"
            ]
        )
        if arrays["truth_column_density_m2"].shape != truth_shape:
            raise ValueError("observable truth-map shape disagrees with summary")
        if arrays["reconstructed_column_density_m2"].shape != reconstructed_shape:
            raise ValueError("observable reconstruction-map shape disagrees with summary")
        leading_shape = reconstructed_shape[:4]
        if arrays["seeds"].shape != leading_shape:
            raise ValueError("observable seed tensor does not match replay axes")
        if arrays["fitted_coefficients"].shape[:4] != leading_shape:
            raise ValueError("observable coefficient tensor does not match replay axes")
        if arrays["readout_names"].shape != (leading_shape[0],):
            raise ValueError("observable readout coordinate does not match map tensor")
        if arrays["fluence_mw_us"].shape != (leading_shape[1],):
            raise ValueError("observable fluence coordinate does not match map tensor")
        if arrays["morphology_names"].shape != (leading_shape[2],):
            raise ValueError("observable morphology coordinate does not match map tensor")
        if arrays["realization_indices"].shape != (leading_shape[3],):
            raise ValueError("observable realization coordinate does not match map tensor")
        if arrays["readout_names"].tolist() != summary["replay"]["readout_names"]:
            raise ValueError("observable readout coordinates disagree with summary")
        if not np.array_equal(
            arrays["fluence_mw_us"],
            np.asarray(summary["replay"]["fluence_mw_us"], dtype=float),
        ):
            raise ValueError("observable fluence coordinates disagree with summary")
        if arrays["morphology_names"].tolist() != summary["replay"][
            "morphology_names"
        ]:
            raise ValueError("observable morphology coordinates disagree with summary")
        if arrays["realization_indices"].tolist() != summary["replay"][
            "realization_indices"
        ]:
            raise ValueError("observable realization coordinates disagree with summary")
        axis_order = np.asarray(arrays["tensor_axis_order"])
        expected_axis_order = ",".join(
            summary["maps_schema"].get(
                "reconstructed_column_density_m2_axis_order",
                [
                    "readout",
                    "fluence",
                    "morphology",
                    "realization",
                    "z_index",
                    "y_index",
                ],
            )
        )
        if axis_order.ndim != 0 or str(axis_order.item()) != expected_axis_order:
            raise ValueError("observable tensor-axis declaration disagrees with summary")
        source_run = np.asarray(arrays["source_run_id"])
        if source_run.ndim != 0 or str(source_run.item()) != str(
            summary["source_benchmark"]["run_id"]
        ):
            raise ValueError("observable maps source run_id disagrees with summary")


def verify_observable_artifact_set(
    output_directory: Path,
) -> VerifiedObservableArtifactSet:
    """Verify one complete sealed artifact set before any figure is rendered."""

    output_directory = output_directory.resolve()
    manifest_path = output_directory / OBSERVABLE_MANIFEST_FILENAME
    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)
    manifest = load_json(manifest_path)
    if manifest.get("schema_version") != OBSERVABLE_MANIFEST_SCHEMA_VERSION:
        raise ValueError("unsupported observable artifact-manifest schema version")
    records = manifest.get("artifacts")
    if not isinstance(records, dict):
        raise ValueError("observable artifact manifest has no artifact records")
    paths: dict[str, Path] = {}
    roles: dict[str, str] = {}
    for filename, record in records.items():
        if not isinstance(record, dict):
            raise ValueError(f"invalid observable artifact record for {filename!r}")
        path = _safe_artifact_path(output_directory, str(filename))
        if not path.is_file():
            raise FileNotFoundError(path)
        expected = str(record.get("sha256", ""))
        actual = file_sha256(path)
        if actual != expected:
            raise ValueError(
                f"observable artifact hash mismatch for {filename}: "
                f"expected {expected}, got {actual}"
            )
        paths[str(filename)] = path
        roles[str(filename)] = str(record.get("role", ""))
    required = set(OBSERVABLE_RAW_ARTIFACT_ROLES) | {OBSERVABLE_SUMMARY_FILENAME}
    missing = required.difference(paths)
    if missing:
        raise ValueError(
            "observable artifact manifest is missing records: "
            + ", ".join(sorted(missing))
        )
    config_record = manifest.get("config")
    if not isinstance(config_record, dict):
        raise ValueError("observable artifact manifest has no config record")
    if config_record.get("artifact") != STUDY_CONFIG_SNAPSHOT_FILENAME:
        raise ValueError("observable artifact manifest names an unexpected config")
    config_sha256 = str(config_record.get("sha256", ""))
    if config_sha256 != file_sha256(paths[STUDY_CONFIG_SNAPSHOT_FILENAME]):
        raise ValueError("observable manifest config hash disagrees with config snapshot")
    raw_hashes = {
        filename: file_sha256(paths[filename])
        for filename in OBSERVABLE_RAW_ARTIFACT_ROLES
    }
    expected_run_id = _deterministic_observable_run_id(
        config_sha256,
        str(manifest.get("source_run_id", "")),
        raw_hashes,
    )
    if manifest.get("run_id") != expected_run_id:
        raise ValueError("observable run_id does not match sealed raw artifacts")
    summary = load_json(paths[OBSERVABLE_SUMMARY_FILENAME])
    config = load_json(paths[STUDY_CONFIG_SNAPSHOT_FILENAME])
    if summary.get("run_id") != expected_run_id:
        raise ValueError("observable summary run_id disagrees with manifest")
    labels_match = (
        summary.get("label") == manifest.get("label")
        and config.get("label") == manifest.get("label")
    )
    if not labels_match:
        raise ValueError("observable label differs across config, summary and manifest")
    if summary["source_benchmark"]["run_id"] != manifest.get("source_run_id"):
        raise ValueError("observable source run_id differs across summary and manifest")
    source_manifest = load_json(paths[SOURCE_MANIFEST_SNAPSHOT_FILENAME])
    if source_manifest.get("run_id") != manifest.get("source_run_id"):
        raise ValueError("source-manifest snapshot run_id disagrees with observable lineage")
    if summary["source_benchmark"].get("manifest_sha256") != file_sha256(
        paths[SOURCE_MANIFEST_SNAPSHOT_FILENAME]
    ):
        raise ValueError("source-manifest snapshot hash disagrees with observable summary")
    if summary["source_benchmark"].get("config_sha256") != file_sha256(
        paths[SOURCE_CONFIG_SNAPSHOT_FILENAME]
    ):
        raise ValueError("source-config snapshot hash disagrees with observable summary")
    if not bool(summary["replay"]["all_legacy_metrics_verified"]):
        raise ValueError("observable summary does not certify legacy metric closure")
    _validate_maps_schema(paths[OBSERVABLE_MAPS_FILENAME], summary)
    return VerifiedObservableArtifactSet(
        directory=output_directory,
        manifest=manifest,
        summary=summary,
        config=config,
        paths=paths,
    )


_PANEL_CONTRACTS = (
    (
        "median_absolute_integrated_response_relative_error",
        "minimum_absolute_integrated_response_relative_error",
        "maximum_absolute_integrated_response_relative_error",
        "Integrated response\nabsolute relative error",
    ),
    (
        "median_centroid_position_error_um",
        "minimum_centroid_position_error_um",
        "maximum_centroid_position_error_um",
        r"Centroid position error ($\mu$m)",
    ),
    (
        "median_absolute_major_rms_width_relative_error",
        "minimum_absolute_major_rms_width_relative_error",
        "maximum_absolute_major_rms_width_relative_error",
        "Major rms width\nabsolute relative error",
    ),
    (
        "median_absolute_minor_rms_width_relative_error",
        "minimum_absolute_minor_rms_width_relative_error",
        "maximum_absolute_minor_rms_width_relative_error",
        "Minor rms width\nabsolute relative error",
    ),
    (
        "median_absolute_aspect_ratio_relative_error",
        "minimum_absolute_aspect_ratio_relative_error",
        "maximum_absolute_aspect_ratio_relative_error",
        "Aspect ratio\nabsolute relative error",
    ),
    (
        "median_absolute_principal_axis_angle_error_deg",
        "minimum_absolute_principal_axis_angle_error_deg",
        "maximum_absolute_principal_axis_angle_error_deg",
        "Principal-axis error (deg)\n(supported pairs only)",
    ),
)


def _plot_observable_recovery(
    verified: VerifiedObservableArtifactSet,
) -> list[Path]:
    rows = load_rows(verified.paths[OBSERVABLE_AGGREGATE_FILENAME])
    if not rows:
        raise ValueError("observable aggregate table is empty")
    expected_groups = {
        (readout, float(fluence))
        for readout in verified.summary["replay"]["readout_names"]
        for fluence in verified.summary["replay"]["fluence_mw_us"]
    }
    actual_groups = {
        (row["readout"], float(row["fluence_mw_us"])) for row in rows
    }
    if actual_groups != expected_groups:
        raise ValueError("observable aggregate axes disagree with sealed run summary")
    expected_trials = (
        len(verified.summary["replay"]["morphology_names"])
        * len(verified.summary["replay"]["realization_indices"])
    )
    if any(int(row["trial_count"]) != expected_trials for row in rows):
        raise ValueError("observable aggregate trial denominator is incomplete")
    if any(int(row["legacy_verified_trials"]) != expected_trials for row in rows):
        raise ValueError("observable aggregate contains unverified replay trials")

    figure, axes = plt.subplots(2, 3, figsize=(10.9, 6.7), constrained_layout=True)
    styles = {
        "dual_port": ("#006BA4", "o", "Dual-port (DPFI)"),
        "dark_field": ("#A23B72", "s", "Dark-field (DFFI)"),
    }
    for axis, contract in zip(axes.flat, _PANEL_CONTRACTS, strict=True):
        median_key, minimum_key, maximum_key, ylabel = contract
        for readout in verified.summary["replay"]["readout_names"]:
            group = sorted(
                (row for row in rows if row["readout"] == readout),
                key=lambda row: float(row["fluence_mw_us"]),
            )
            finite = [row for row in group if row.get(median_key, "") != ""]
            if not finite:
                continue
            x = np.asarray([float(row["fluence_mw_us"]) for row in finite])
            median = np.asarray([float(row[median_key]) for row in finite])
            minimum = np.asarray([float(row[minimum_key]) for row in finite])
            maximum = np.asarray([float(row[maximum_key]) for row in finite])
            colour, marker, label = styles[readout]
            axis.errorbar(
                x,
                median,
                yerr=np.vstack((median - minimum, maximum - median)),
                color=colour,
                marker=marker,
                linewidth=1.5,
                capsize=3.0,
                label=label,
            )
        axis.set_xlabel(r"Fluence $F$ (mW $\mu$s)")
        axis.set_ylabel(ylabel)
        axis.set_xticks(verified.summary["replay"]["fluence_mw_us"])
        axis.set_ylim(bottom=0.0)
        axis.grid(alpha=0.22, linewidth=0.7)
    axes[0, 0].legend(frameon=False, loc="best")
    figure.suptitle(
        "Held-out physical-observable recovery from independent frozen fits\n"
        "markers: median; bars: observed finite-ensemble range (not a confidence interval)",
        fontsize=12,
    )
    figure_config = verified.config["figure"]
    basename = str(figure_config["basename"])
    outputs: list[Path] = []
    for extension in figure_config["formats"]:
        path = verified.directory / f"{basename}.{extension}"
        if extension == "png":
            figure.savefig(path, dpi=300)
        else:
            figure.savefig(path)
        if extension == "svg":
            svg_text = path.read_text(encoding="utf-8")
            path.write_text(
                "\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n",
                encoding="utf-8",
            )
        outputs.append(path)
    plt.close(figure)
    return outputs


def generate_observable_benchmark_figures(
    output_directory: Path,
) -> dict[str, Path]:
    """Verify sealed numerical outputs and render the compact Chapter 7 figure."""

    verified = verify_observable_artifact_set(output_directory)
    figure_paths = _plot_observable_recovery(verified)
    metadata_path = verified.directory / FIGURE_METADATA_FILENAME
    write_json(
        metadata_path,
        {
            "source_run_id": verified.manifest["run_id"],
            "source_manifest": OBSERVABLE_MANIFEST_FILENAME,
            "source_manifest_sha256": file_sha256(
                verified.directory / OBSERVABLE_MANIFEST_FILENAME
            ),
            "source_artifact_sha256": {
                filename: record["sha256"]
                for filename, record in verified.manifest["artifacts"].items()
            },
            "source_benchmark": verified.summary["source_benchmark"],
            "source_aggregate_table": OBSERVABLE_AGGREGATE_FILENAME,
            "renderer": Path(__file__).name,
            "renderer_sha256": file_sha256(Path(__file__)),
            "finite_ensemble_display": verified.config["figure"][
                "finite_ensemble_display"
            ],
            "angle_anisotropy_threshold": verified.summary["integration_support"][
                "angle_anisotropy_threshold"
            ],
            "claims_boundary": verified.summary["claims_boundary"],
            "output_files": [path.name for path in figure_paths],
            "rendered_output_sha256": {
                path.name: file_sha256(path) for path in figure_paths
            },
        },
    )
    png = next((path for path in figure_paths if path.suffix == ".png"), figure_paths[0])
    return {
        "figure": png,
        "metadata": metadata_path,
    }


__all__ = [
    "FIGURE_METADATA_FILENAME",
    "VerifiedObservableArtifactSet",
    "generate_observable_benchmark_figures",
    "verify_observable_artifact_set",
]
