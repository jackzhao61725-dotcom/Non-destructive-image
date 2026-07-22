"""Generate benchmark figures exclusively from frozen disk artifacts.

This module deliberately imports no measurement, forward-model or fitting code.
It may format existing tables and arrays, but it cannot resimulate an observation
or reconstruct a density map.  The figures therefore remain tied to the exact run
recorded by :mod:`non_destructive_image.reconstruction.studies.artifacts`.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .io import file_sha256, load_json, load_rows, write_json


_QUALITY_BASENAME = "reconstruction_quality_vs_fluence"
_MAP_BASENAME = "representative_three_peak_reconstruction_F90"
_MANIFEST_FILENAME = "artifact_manifest.json"
_MANIFEST_SCHEMA_VERSION = 1
_CONFIG_SNAPSHOT_FILENAME = "study_config.json"
_SUMMARY_FILENAME = "benchmark_summary.json"
_REPRESENTATIVE_FILENAME = "representative_reconstruction.npz"
_RAW_ARTIFACT_FILENAMES = (
    "calibration_candidate_summary.csv",
    "held_out_trials.csv",
    "held_out_summary.csv",
    "reduced_grid_convergence.csv",
    _REPRESENTATIVE_FILENAME,
)
_REQUIRED_ARTIFACT_FILENAMES = {
    _CONFIG_SNAPSHOT_FILENAME,
    *_RAW_ARTIFACT_FILENAMES,
    _SUMMARY_FILENAME,
}


def _required_file(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"required benchmark artifact does not exist: {path}")
    return path


def _digest_is_valid(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _artifact_path(output_directory: Path, filename: str) -> Path:
    relative = Path(filename)
    if relative.is_absolute() or len(relative.parts) != 1 or relative.name != filename:
        raise ValueError(f"manifest artifact name is not a safe basename: {filename!r}")
    return output_directory / relative


def _expected_run_id(
    config_sha256: str,
    artifact_records: dict[str, Any],
) -> str:
    raw_hashes = {
        filename: str(artifact_records[filename]["sha256"])
        for filename in _RAW_ARTIFACT_FILENAMES
    }
    payload = json.dumps(
        {
            "config_sha256": config_sha256,
            "raw_artifact_sha256": raw_hashes,
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _verified_artifact_set(
    output_directory: Path,
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Path]]:
    """Load a manifest and fail unless every declared source hash is valid."""

    manifest_path = _required_file(output_directory / _MANIFEST_FILENAME)
    manifest = load_json(manifest_path)
    if manifest.get("schema_version") != _MANIFEST_SCHEMA_VERSION:
        raise ValueError("unsupported or missing artifact-manifest schema version")
    run_id = manifest.get("run_id")
    if not _digest_is_valid(run_id):
        raise ValueError("artifact manifest has an invalid run_id")

    artifact_records = manifest.get("artifacts")
    if not isinstance(artifact_records, dict):
        raise ValueError("artifact manifest has no artifact records")
    missing = _REQUIRED_ARTIFACT_FILENAMES.difference(artifact_records)
    if missing:
        raise ValueError(
            "artifact manifest is missing required records: "
            + ", ".join(sorted(missing))
        )

    paths: dict[str, Path] = {}
    for filename, record in artifact_records.items():
        if not isinstance(record, dict) or not _digest_is_valid(record.get("sha256")):
            raise ValueError(f"invalid hash record for artifact {filename!r}")
        path = _required_file(_artifact_path(output_directory, str(filename)))
        actual = file_sha256(path)
        expected = str(record["sha256"])
        if actual != expected:
            raise ValueError(
                f"artifact hash mismatch for {filename}: expected {expected}, got {actual}"
            )
        paths[str(filename)] = path

    config_record = manifest.get("config")
    if not isinstance(config_record, dict):
        raise ValueError("artifact manifest has no config record")
    if config_record.get("artifact") != _CONFIG_SNAPSHOT_FILENAME:
        raise ValueError("manifest config does not name the frozen config snapshot")
    config_sha256 = config_record.get("sha256")
    if not _digest_is_valid(config_sha256):
        raise ValueError("manifest config hash is invalid")
    if artifact_records[_CONFIG_SNAPSHOT_FILENAME]["sha256"] != config_sha256:
        raise ValueError("config record and config-artifact hashes disagree")
    if _expected_run_id(str(config_sha256), artifact_records) != run_id:
        raise ValueError("manifest run_id does not match its config and raw artifacts")

    config_snapshot = load_json(paths[_CONFIG_SNAPSHOT_FILENAME])
    if config_snapshot.get("label") != manifest.get("label"):
        raise ValueError("frozen config label does not match artifact manifest")
    summary = load_json(paths[_SUMMARY_FILENAME])
    if summary.get("label") != manifest.get("label"):
        raise ValueError("benchmark summary label does not match artifact manifest")
    if summary.get("run_id") != run_id:
        raise ValueError("benchmark summary run_id does not match artifact manifest")
    representative = summary.get("representative")
    if not isinstance(representative, dict):
        raise ValueError("benchmark summary has no representative record")
    if representative.get("artifact") != _REPRESENTATIVE_FILENAME:
        raise ValueError("summary representative does not name the sealed NPZ artifact")
    return manifest, summary, paths


def _required_scalar(arrays: Any, key: str) -> Any:
    if key not in arrays.files:
        raise ValueError(f"representative NPZ is missing metadata scalar: {key}")
    value = np.asarray(arrays[key])
    if value.ndim != 0:
        raise ValueError(f"representative NPZ metadata is not scalar: {key}")
    return value.item()


def _validate_representative_metadata(
    arrays: Any,
    representative: dict[str, Any],
) -> None:
    """Prove that summary metadata describes the exact frozen NPZ realization."""

    required_summary_fields = (
        "morphology_name",
        "fluence_mw_us",
        "realization_index",
        "seeds",
        "supported_band_errors",
    )
    missing = [key for key in required_summary_fields if key not in representative]
    if missing:
        raise ValueError(
            "summary representative is missing fields: " + ", ".join(missing)
        )
    if str(_required_scalar(arrays, "morphology_name")) != str(
        representative["morphology_name"]
    ):
        raise ValueError("representative morphology differs between summary and NPZ")
    if int(_required_scalar(arrays, "realization_index")) != int(
        representative["realization_index"]
    ):
        raise ValueError("representative realization index differs between summary and NPZ")
    if not np.isclose(
        float(_required_scalar(arrays, "fluence_mw_us")),
        float(representative["fluence_mw_us"]),
        rtol=1e-12,
        atol=1e-12,
    ):
        raise ValueError("representative fluence differs between summary and NPZ")

    seeds = representative["seeds"]
    errors = representative["supported_band_errors"]
    if not isinstance(seeds, dict) or not isinstance(errors, dict):
        raise ValueError("representative seeds and errors must be mappings")
    for readout in ("dual_port", "dark_field"):
        if readout not in seeds or readout not in errors:
            raise ValueError(f"summary representative lacks {readout} seed or error")
        if int(_required_scalar(arrays, f"{readout}_seed")) != int(seeds[readout]):
            raise ValueError(
                f"representative {readout} seed differs between summary and NPZ"
            )
        if not np.isclose(
            float(_required_scalar(arrays, f"{readout}_supported_band_error")),
            float(errors[readout]),
            rtol=1e-12,
            atol=1e-12,
        ):
            raise ValueError(
                f"representative {readout} error differs between summary and NPZ"
            )


def _save_figure(figure: plt.Figure, directory: Path, basename: str) -> list[Path]:
    outputs: list[Path] = []
    for suffix in ("svg", "png", "pdf"):
        path = directory / f"{basename}.{suffix}"
        figure.savefig(path, dpi=300)
        outputs.append(path)
    return outputs


def _plot_quality_curve(
    rows: list[dict[str, str]],
    output_directory: Path,
    reference_fluence: float,
) -> list[Path]:
    fig, axis = plt.subplots(figsize=(5.3, 3.7), constrained_layout=True)
    styles = {
        "dual_port": ("#1565c0", "o", "Dual-port"),
        "dark_field": ("#d84315", "s", "Dark-field"),
    }
    for readout, (colour, marker, label) in styles.items():
        selected = sorted(
            (
                row
                for row in rows
                if row["readout"] == readout
                and row.get("median_supported_band_relative_l2_error", "") != ""
            ),
            key=lambda row: float(row["fluence_mw_us"]),
        )
        if not selected:
            raise ValueError(f"held-out summary has no rows for {readout}")
        axis.plot(
            [float(row["fluence_mw_us"]) for row in selected],
            [
                float(row["median_supported_band_relative_l2_error"])
                for row in selected
            ],
            color=colour,
            marker=marker,
            linewidth=2.0,
            markersize=5.5,
            label=label,
        )
    axis.axvline(reference_fluence, color="0.45", linestyle="--", linewidth=1.1)
    axis.text(
        reference_fluence + 2.5,
        0.98,
        "reference",
        color="0.38",
        fontsize=8.5,
        transform=axis.get_xaxis_transform(),
        va="top",
    )
    axis.set_xlabel(r"Power--time product $F$ (mW $\mu$s)")
    axis.set_ylabel("Median supported-band relative error")
    axis.set_ylim(bottom=0.0)
    axis.legend(frameon=False)
    axis.grid(alpha=0.22, linewidth=0.7)
    outputs = _save_figure(fig, output_directory, _QUALITY_BASENAME)
    plt.close(fig)
    return outputs


def _plot_reconstruction_maps(
    arrays: Any,
    representative: dict[str, Any],
    output_directory: Path,
) -> tuple[list[Path], dict[str, Any]]:
    required_arrays = (
        "y_grid_m",
        "z_grid_m",
        "truth_column_density_m2",
        "dual_port_column_density_m2",
        "dark_field_column_density_m2",
    )
    missing = [key for key in required_arrays if key not in arrays.files]
    if missing:
        raise ValueError(f"representative NPZ is missing arrays: {', '.join(missing)}")

    errors = dict(representative["supported_band_errors"])
    morphology_name = str(representative["morphology_name"])
    fluence = float(representative["fluence_mw_us"])
    realization = int(representative["realization_index"])

    reconstructions = {
        "truth": np.asarray(arrays["truth_column_density_m2"], dtype=float),
        "dual_port": np.asarray(arrays["dual_port_column_density_m2"], dtype=float),
        "dark_field": np.asarray(arrays["dark_field_column_density_m2"], dtype=float),
    }
    display = (
        ("truth", "Held-out truth"),
        ("dual_port", "Dual-port reconstruction"),
        ("dark_field", "Dark-field reconstruction"),
    )
    common_maximum = max(float(np.max(reconstructions[key])) for key, _ in display) / 1e14
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(8.2, 2.65),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    y_grid_m = np.asarray(arrays["y_grid_m"], dtype=float)
    z_grid_m = np.asarray(arrays["z_grid_m"], dtype=float)
    extent = [
        float(np.min(y_grid_m) * 1e6),
        float(np.max(y_grid_m) * 1e6),
        float(np.min(z_grid_m) * 1e6),
        float(np.max(z_grid_m) * 1e6),
    ]
    image = None
    for axis, (key, title) in zip(axes, display, strict=True):
        image = axis.imshow(
            reconstructions[key] / 1e14,
            origin="lower",
            extent=extent,
            cmap="viridis",
            vmin=0.0,
            vmax=common_maximum,
            interpolation="nearest",
            aspect="auto",
        )
        axis.set_title(title, fontsize=9.5)
        axis.set_xlim(-30.0, 30.0)
        axis.set_ylim(-8.0, 8.0)
        axis.set_xlabel(r"$y$ ($\mu$m)")
        error = errors.get(key)
        if key != "truth" and error is not None:
            axis.text(
                0.03,
                0.05,
                rf"$\epsilon_{{\rm band}}={float(error):.2f}$",
                transform=axis.transAxes,
                fontsize=8.2,
                color="white",
                bbox={
                    "facecolor": "black",
                    "alpha": 0.45,
                    "edgecolor": "none",
                    "pad": 1.5,
                },
            )
    axes[0].set_ylabel(r"$z$ ($\mu$m)")
    if image is None:
        raise RuntimeError("reconstruction map figure did not create an image")
    colourbar = fig.colorbar(image, ax=axes, fraction=0.025, pad=0.02)
    colourbar.set_label(r"Column density ($10^{14}\,\mathrm{m}^{-2}$)")
    outputs = _save_figure(fig, output_directory, _MAP_BASENAME)
    plt.close(fig)
    return outputs, {
        "morphology_name": str(morphology_name),
        "fluence_mw_us": fluence,
        "realization_index": realization,
        "supported_band_errors": errors,
    }


def generate_morphology_benchmark_figures(
    output_directory: Path,
) -> dict[str, Path]:
    """Render quality and representative-map figures from frozen artifacts only."""

    output_directory = output_directory.resolve()
    manifest, summary, paths = _verified_artifact_set(output_directory)
    summary_path = paths[_SUMMARY_FILENAME]
    held_out_summary_path = paths["held_out_summary.csv"]
    representative_path = paths[_REPRESENTATIVE_FILENAME]
    rows = load_rows(held_out_summary_path)
    if not rows:
        raise ValueError("held_out_summary.csv contains no rows")

    representative = dict(summary["representative"])
    reference_fluence = float(representative["fluence_mw_us"])
    with np.load(representative_path, allow_pickle=False) as arrays:
        _validate_representative_metadata(arrays, representative)
        quality_outputs = _plot_quality_curve(
            rows, output_directory, reference_fluence
        )
        map_outputs, representative_metadata = _plot_reconstruction_maps(
            arrays,
            representative,
            output_directory,
        )

    metadata_path = output_directory / "figure_metadata.json"
    rendered_outputs = [*quality_outputs, *map_outputs]
    write_json(
        metadata_path,
        {
            "source_run_id": manifest["run_id"],
            "source_manifest": _MANIFEST_FILENAME,
            "renderer_sha256": file_sha256(Path(__file__)),
            "rendered_output_sha256": {
                path.name: file_sha256(path) for path in rendered_outputs
            },
            "source_artifact_sha256": {
                filename: record["sha256"]
                for filename, record in manifest["artifacts"].items()
            },
            "source_summary": summary_path.name,
            "source_held_out_summary": held_out_summary_path.name,
            "source_representative_arrays": representative_path.name,
            "representative_morphology": representative_metadata[
                "morphology_name"
            ],
            "representative_fluence_mw_us": representative_metadata[
                "fluence_mw_us"
            ],
            "representative_realization_index": representative_metadata[
                "realization_index"
            ],
            "representative_seeds": representative.get("seeds", {}),
            "representative_supported_band_errors": representative_metadata[
                "supported_band_errors"
            ],
            "quality_curve_statistic": (
                "median over the finite held-out development ensemble; "
                "no uncertainty interval claimed"
            ),
            "representative": representative_metadata,
            "absolute_faraday_scale": summary.get("interpretation", {}).get(
                "faraday_absolute_scale", "not recorded"
            ),
            "output_files": [
                held_out_summary_path.name,
                *(path.name for path in rendered_outputs),
                representative_path.name,
            ],
        },
    )
    return {
        "quality_curve": quality_outputs[0],
        "representative_map": map_outputs[0],
        "metadata": metadata_path,
    }
