"""Render the DPFI initial-condition sweep from sealed disk artifacts only.

The renderer never simulates a camera observation or runs the inverse.  It
first verifies every manifest-recorded artifact and the deterministic raw-run
identity, then reads the frozen tables and arrays needed for the figures.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import textwrap
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from .initial_condition_suite import (
    CONDITION_SUMMARY_FILENAME,
    MANIFEST_FILENAME,
    MANIFEST_SCHEMA_VERSION,
    MAPS_FILENAME,
    RAW_ARTIFACT_ROLES,
    SOURCE_CONFIG_FILENAME,
    SOURCE_MANIFEST_FILENAME,
    STUDY_CONFIG_FILENAME,
    SUMMARY_FILENAME,
    TRIALS_FILENAME,
    _deterministic_run_id,
)
from .io import file_sha256, load_json, load_rows, write_json


FIGURE_METADATA_FILENAME = "figure_metadata.json"


@dataclass(frozen=True)
class VerifiedInitialConditionArtifactSet:
    """Paths and metadata proven to belong to one sealed suite run."""

    directory: Path
    manifest: dict[str, Any] = field(repr=False)
    summary: dict[str, Any] = field(repr=False)
    config: dict[str, Any] = field(repr=False)
    paths: dict[str, Path] = field(repr=False)
    trial_rows: tuple[dict[str, str], ...] = field(repr=False)


def _safe_artifact_path(directory: Path, filename: str) -> Path:
    candidate = Path(filename)
    if candidate.is_absolute() or candidate.name != filename:
        raise ValueError(f"artifact name is not a safe basename: {filename!r}")
    return directory / candidate


def _parse_bool(value: str) -> bool:
    normalised = str(value).strip().lower()
    if normalised in {"true", "1", "yes"}:
        return True
    if normalised in {"false", "0", "no"}:
        return False
    raise ValueError(f"invalid boolean value in sealed suite table: {value!r}")


def _suite_axes(
    summary: dict[str, Any],
    config: dict[str, Any],
) -> tuple[tuple[str, ...], tuple[float, ...], tuple[int, ...]]:
    suite = summary.get("suite")
    if not isinstance(suite, dict):
        raise ValueError("suite summary has no suite record")
    condition_ids = tuple(str(value) for value in suite.get("condition_ids", ()))
    fluences = tuple(float(value) for value in suite.get("fluence_mw_us", ()))
    realizations = tuple(int(value) for value in suite.get("realization_indices", ()))
    if not condition_ids or not fluences or not realizations:
        raise ValueError("suite summary declares an empty coordinate axis")
    if len(condition_ids) != len(set(condition_ids)):
        raise ValueError("suite summary condition ids are not unique")
    if len(fluences) != len(set(fluences)):
        raise ValueError("suite summary fluences are not unique")
    if len(realizations) != len(set(realizations)):
        raise ValueError("suite summary realization indices are not unique")

    configured_conditions = config.get("initial_conditions")
    if not isinstance(configured_conditions, list):
        raise ValueError("suite config has no initial_conditions list")
    configured_ids = tuple(str(item.get("id", "")) for item in configured_conditions)
    if configured_ids != condition_ids:
        raise ValueError("condition coordinates differ between config and summary")
    ensemble = config.get("ensemble")
    if not isinstance(ensemble, dict):
        raise ValueError("suite config has no ensemble record")
    if tuple(float(value) for value in ensemble.get("fluence_mw_us", ())) != fluences:
        raise ValueError("fluence coordinates differ between config and summary")
    expected_realizations = tuple(
        range(int(ensemble.get("realizations_per_condition", 0)))
    )
    if expected_realizations != realizations:
        raise ValueError("realization coordinates differ between config and summary")

    expected_trial_count = len(condition_ids) * len(fluences) * len(realizations)
    if int(suite.get("trial_count", -1)) != expected_trial_count:
        raise ValueError("suite summary trial count disagrees with its coordinate axes")
    if suite.get("readout") != "dual_port":
        raise ValueError("initial-condition suite is not the declared DPFI readout")
    usability = summary.get("observable_usability")
    if not isinstance(usability, dict):
        raise ValueError("suite summary has no observable-usability record")
    if usability.get("scores_are_reported_independently") is not True:
        raise ValueError("suite summary does not preserve independent observable scores")
    if usability.get("overall_image_usability_score", object()) is not None:
        raise ValueError("suite summary unexpectedly defines an overall usability score")
    return condition_ids, fluences, realizations


def _validate_trial_table(
    rows: list[dict[str, str]],
    condition_ids: tuple[str, ...],
    fluences: tuple[float, ...],
    realizations: tuple[int, ...],
) -> dict[tuple[float, str, int], dict[str, str]]:
    required_columns = {
        "condition_id",
        "fluence_mw_us",
        "realization_index",
        "seed",
        "fit_success",
        "c_A",
        "c_r",
        "c_w",
    }
    if not rows or not required_columns.issubset(rows[0]):
        missing = required_columns.difference(rows[0] if rows else {})
        raise ValueError(
            "condition trial table is missing columns: " + ", ".join(sorted(missing))
        )
    expected = {
        (fluence, condition_id, realization)
        for fluence in fluences
        for condition_id in condition_ids
        for realization in realizations
    }
    by_key: dict[tuple[float, str, int], dict[str, str]] = {}
    for row in rows:
        key = (
            float(row["fluence_mw_us"]),
            str(row["condition_id"]),
            int(row["realization_index"]),
        )
        if key in by_key:
            raise ValueError(f"duplicate condition trial coordinate: {key}")
        _parse_bool(row["fit_success"])
        int(row["seed"])
        for score_name in ("c_A", "c_r", "c_w"):
            raw_value = row.get(score_name, "")
            if raw_value == "":
                continue
            score = float(raw_value)
            if not np.isfinite(score) or score < 0.0:
                raise ValueError(f"invalid {score_name} value at {key}: {raw_value!r}")
        by_key[key] = row
    if set(by_key) != expected:
        missing = sorted(expected.difference(by_key))
        unexpected = sorted(set(by_key).difference(expected))
        raise ValueError(
            "condition trial axes disagree with the sealed suite; "
            f"missing={missing}, unexpected={unexpected}"
        )
    return by_key


def _validate_condition_summary(
    path: Path,
    condition_ids: tuple[str, ...],
    trial_rows: dict[tuple[float, str, int], dict[str, str]],
) -> None:
    rows = load_rows(path)
    if not rows or "condition_id" not in rows[0]:
        raise ValueError("condition summary table has no condition_id column")
    by_id: dict[str, dict[str, str]] = {}
    for row in rows:
        identifier = str(row["condition_id"])
        if identifier in by_id:
            raise ValueError(f"duplicate condition summary row: {identifier!r}")
        by_id[identifier] = row
    if set(by_id) != set(condition_ids):
        raise ValueError("condition summary ids disagree with suite coordinates")
    for identifier, row in by_id.items():
        matching = [
            trial for key, trial in trial_rows.items() if key[1] == identifier
        ]
        if int(row.get("trial_count", -1)) != len(matching):
            raise ValueError(f"condition summary trial count is wrong for {identifier!r}")
        successful = sum(_parse_bool(trial["fit_success"]) for trial in matching)
        if int(row.get("successful_fit_count", -1)) != successful:
            raise ValueError(
                f"condition summary successful-fit count is wrong for {identifier!r}"
            )


def _validate_maps_schema(
    path: Path,
    summary: dict[str, Any],
    condition_ids: tuple[str, ...],
    fluences: tuple[float, ...],
    realizations: tuple[int, ...],
    trial_rows: dict[tuple[float, str, int], dict[str, str]],
) -> None:
    required = {
        "y_grid_m",
        "z_grid_m",
        "cell_area_m2",
        "support_mask",
        "fluence_mw_us",
        "condition_ids",
        "realization_indices",
        "candidate_label",
        "seeds",
        "fit_success",
        "truth_column_density_m2",
        "reconstructed_column_density_m2",
        "fitted_coefficients",
        "dual_port_h_counts_e",
        "dual_port_v_counts_e",
        "source_run_id",
        "tensor_axis_order",
    }
    with np.load(path, allow_pickle=False) as arrays:
        missing = required.difference(arrays.files)
        if missing:
            raise ValueError(
                "condition maps NPZ is missing keys: " + ", ".join(sorted(missing))
            )
        grid_shape = tuple(
            int(value) for value in summary["integration_support"]["grid_shape"]
        )
        if len(grid_shape) != 2:
            raise ValueError("condition maps require a two-dimensional object grid")
        for name in ("y_grid_m", "z_grid_m", "cell_area_m2", "support_mask"):
            if arrays[name].shape != grid_shape:
                raise ValueError(f"condition maps {name} shape disagrees with summary")
        if int(np.count_nonzero(arrays["support_mask"])) != int(
            summary["integration_support"]["supported_cell_count"]
        ):
            raise ValueError("condition maps support mask disagrees with summary")

        leading_shape = (len(fluences), len(condition_ids), len(realizations))
        if arrays["truth_column_density_m2"].shape != (len(condition_ids), *grid_shape):
            raise ValueError("condition truth-map tensor shape disagrees with suite axes")
        if arrays["reconstructed_column_density_m2"].shape != (
            *leading_shape,
            *grid_shape,
        ):
            raise ValueError(
                "condition reconstruction-map tensor shape disagrees with suite axes"
            )
        for name in ("seeds", "fit_success"):
            if arrays[name].shape != leading_shape:
                raise ValueError(f"condition maps {name} tensor shape is invalid")
        parameter_count = int(summary["suite"]["parameter_count"])
        if arrays["fitted_coefficients"].shape != (*leading_shape, parameter_count):
            raise ValueError("condition coefficient tensor shape disagrees with summary")
        h_shape = arrays["dual_port_h_counts_e"].shape
        v_shape = arrays["dual_port_v_counts_e"].shape
        if h_shape != v_shape or h_shape[:3] != leading_shape or len(h_shape) != 5:
            raise ValueError("condition raw-channel tensors have invalid axes")

        if arrays["condition_ids"].tolist() != list(condition_ids):
            raise ValueError("condition coordinates disagree between maps and summary")
        if not np.array_equal(
            arrays["fluence_mw_us"], np.asarray(fluences, dtype=float)
        ):
            raise ValueError("fluence coordinates disagree between maps and summary")
        if arrays["realization_indices"].tolist() != list(realizations):
            raise ValueError("realization coordinates disagree between maps and summary")
        candidate = np.asarray(arrays["candidate_label"])
        if candidate.ndim != 0 or str(candidate.item()) != str(
            summary["suite"]["selected_candidate"]
        ):
            raise ValueError("condition maps candidate label disagrees with summary")
        source_run = np.asarray(arrays["source_run_id"])
        if source_run.ndim != 0 or str(source_run.item()) != str(
            summary["source_benchmark"]["run_id"]
        ):
            raise ValueError("condition maps source run_id disagrees with summary")
        axis_order = np.asarray(arrays["tensor_axis_order"])
        expected_axis_order = ",".join(summary["maps_schema"]["reconstruction_axis_order"])
        if axis_order.ndim != 0 or str(axis_order.item()) != expected_axis_order:
            raise ValueError("condition maps tensor-axis declaration disagrees with summary")

        for fluence_index, fluence in enumerate(fluences):
            for condition_index, identifier in enumerate(condition_ids):
                for realization_index, realization in enumerate(realizations):
                    row = trial_rows[(fluence, identifier, realization)]
                    array_index = (fluence_index, condition_index, realization_index)
                    if int(arrays["seeds"][array_index]) != int(row["seed"]):
                        raise ValueError("condition map seed disagrees with trial table")
                    if bool(arrays["fit_success"][array_index]) != _parse_bool(
                        row["fit_success"]
                    ):
                        raise ValueError(
                            "condition map fit-success flag disagrees with trial table"
                        )


def verify_initial_condition_artifact_set(
    output_directory: Path,
) -> VerifiedInitialConditionArtifactSet:
    """Verify a complete suite artifact set before parsing its raw payloads."""

    output_directory = output_directory.resolve()
    manifest_path = output_directory / MANIFEST_FILENAME
    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)
    manifest = load_json(manifest_path)
    if manifest.get("schema_version") != MANIFEST_SCHEMA_VERSION:
        raise ValueError("unsupported initial-condition artifact-manifest schema version")
    records = manifest.get("artifacts")
    if not isinstance(records, dict):
        raise ValueError("initial-condition artifact manifest has no artifact records")

    # Do not parse any raw JSON, CSV or NPZ until every recorded digest passes.
    paths: dict[str, Path] = {}
    for raw_filename, raw_record in records.items():
        filename = str(raw_filename)
        if not isinstance(raw_record, dict):
            raise ValueError(f"invalid initial-condition artifact record for {filename!r}")
        expected_hash = str(raw_record.get("sha256", ""))
        if len(expected_hash) != 64 or any(
            character not in "0123456789abcdef" for character in expected_hash
        ):
            raise ValueError(f"invalid initial-condition artifact hash for {filename!r}")
        path = _safe_artifact_path(output_directory, filename)
        if not path.is_file():
            raise FileNotFoundError(path)
        actual_hash = file_sha256(path)
        if actual_hash != expected_hash:
            raise ValueError(
                f"initial-condition artifact hash mismatch for {filename}: "
                f"expected {expected_hash}, got {actual_hash}"
            )
        paths[filename] = path

    required = set(RAW_ARTIFACT_ROLES) | {SUMMARY_FILENAME}
    missing = required.difference(paths)
    if missing:
        raise ValueError(
            "initial-condition manifest is missing records: "
            + ", ".join(sorted(missing))
        )
    for filename, expected_role in RAW_ARTIFACT_ROLES.items():
        if records[filename].get("role") != expected_role:
            raise ValueError(f"initial-condition artifact role is wrong for {filename}")
    if records[SUMMARY_FILENAME].get("role") != "run_summary":
        raise ValueError("initial-condition suite summary has the wrong artifact role")

    config_record = manifest.get("config")
    if not isinstance(config_record, dict):
        raise ValueError("initial-condition manifest has no config record")
    if config_record.get("artifact") != STUDY_CONFIG_FILENAME:
        raise ValueError("initial-condition manifest names an unexpected config")
    config_sha256 = str(config_record.get("sha256", ""))
    if config_sha256 != file_sha256(paths[STUDY_CONFIG_FILENAME]):
        raise ValueError("initial-condition manifest config hash disagrees with snapshot")
    source_run_id = str(manifest.get("source_run_id", ""))
    raw_hashes = {
        filename: file_sha256(paths[filename]) for filename in RAW_ARTIFACT_ROLES
    }
    expected_run_id = _deterministic_run_id(
        config_sha256,
        source_run_id,
        raw_hashes,
    )
    if manifest.get("run_id") != expected_run_id:
        raise ValueError("initial-condition run_id does not match sealed raw artifacts")

    summary = load_json(paths[SUMMARY_FILENAME])
    config = load_json(paths[STUDY_CONFIG_FILENAME])
    if summary.get("run_id") != expected_run_id:
        raise ValueError("initial-condition summary run_id disagrees with manifest")
    if not (
        summary.get("label") == manifest.get("label") == config.get("label")
    ):
        raise ValueError("initial-condition label differs across config, summary and manifest")
    if summary["source_benchmark"]["run_id"] != source_run_id:
        raise ValueError("initial-condition source run_id differs across summary and manifest")
    source_manifest = load_json(paths[SOURCE_MANIFEST_FILENAME])
    if source_manifest.get("run_id") != source_run_id:
        raise ValueError("source-manifest snapshot disagrees with suite lineage")
    if summary["source_benchmark"].get("manifest_sha256") != file_sha256(
        paths[SOURCE_MANIFEST_FILENAME]
    ):
        raise ValueError("source-manifest snapshot hash disagrees with suite summary")
    if summary["source_benchmark"].get("config_sha256") != file_sha256(
        paths[SOURCE_CONFIG_FILENAME]
    ):
        raise ValueError("source-config snapshot hash disagrees with suite summary")

    condition_ids, fluences, realizations = _suite_axes(summary, config)
    rows = load_rows(paths[TRIALS_FILENAME])
    rows_by_key = _validate_trial_table(rows, condition_ids, fluences, realizations)
    successful_count = sum(
        _parse_bool(row["fit_success"]) for row in rows_by_key.values()
    )
    if int(summary["suite"].get("successful_fit_count", -1)) != successful_count:
        raise ValueError("suite successful-fit count disagrees with trial table")
    _validate_condition_summary(
        paths[CONDITION_SUMMARY_FILENAME],
        condition_ids,
        rows_by_key,
    )
    _validate_maps_schema(
        paths[MAPS_FILENAME],
        summary,
        condition_ids,
        fluences,
        realizations,
        rows_by_key,
    )
    return VerifiedInitialConditionArtifactSet(
        directory=output_directory,
        manifest=manifest,
        summary=summary,
        config=config,
        paths=paths,
        trial_rows=tuple(rows),
    )


def _condition_label(identifier: str) -> str:
    return "\n".join(textwrap.wrap(identifier.replace("_", " "), width=28))


def _shared_density_limit(
    truth_column_density_m2: np.ndarray,
    reconstructed_column_density_m2: np.ndarray,
) -> float:
    finite_values = np.concatenate(
        (
            np.asarray(truth_column_density_m2, dtype=float).ravel(),
            np.asarray(reconstructed_column_density_m2, dtype=float).ravel(),
        )
    )
    finite_values = finite_values[np.isfinite(finite_values)]
    if finite_values.size == 0:
        return 1.0
    upper = float(np.max(finite_values))
    return upper if upper > 0.0 else 1.0


def _save_figure(
    figure: plt.Figure,
    directory: Path,
    basename: str,
    formats: list[str],
) -> list[Path]:
    outputs: list[Path] = []
    for extension in formats:
        path = directory / f"{basename}.{extension}"
        if extension == "png":
            figure.savefig(path, dpi=300)
        else:
            figure.savefig(path)
        outputs.append(path)
    plt.close(figure)
    return outputs


def _plot_condition_overview(
    verified: VerifiedInitialConditionArtifactSet,
) -> list[Path]:
    figure_config = verified.config["figure"]
    representative_fluence = float(figure_config["representative_fluence_mw_us"])
    condition_ids, fluences, realizations = _suite_axes(
        verified.summary,
        verified.config,
    )
    matching = [
        index
        for index, value in enumerate(fluences)
        if np.isclose(value, representative_fluence, rtol=0.0, atol=1e-12)
    ]
    if len(matching) != 1:
        raise ValueError("representative fluence is not unique on the sealed suite axis")
    if len(realizations) != 1:
        raise ValueError("version-1 overview requires one realization per condition")
    fluence_index = matching[0]

    with np.load(verified.paths[MAPS_FILENAME], allow_pickle=False) as arrays:
        support = np.asarray(arrays["support_mask"], dtype=bool)
        supported_rows, supported_columns = np.nonzero(support)
        if supported_rows.size == 0:
            raise ValueError("condition overview support is empty")
        row_slice = slice(int(supported_rows.min()), int(supported_rows.max()) + 1)
        column_slice = slice(
            int(supported_columns.min()), int(supported_columns.max()) + 1
        )
        cropped_support = support[row_slice, column_slice]
        y_grid_um = np.asarray(arrays["y_grid_m"], dtype=float)[
            row_slice, column_slice
        ] * 1e6
        z_grid_um = np.asarray(arrays["z_grid_m"], dtype=float)[
            row_slice, column_slice
        ] * 1e6
        extent = (
            float(np.min(y_grid_um)),
            float(np.max(y_grid_um)),
            float(np.min(z_grid_um)),
            float(np.max(z_grid_um)),
        )
        truth_maps = np.asarray(arrays["truth_column_density_m2"], dtype=float)
        reconstruction_maps = np.asarray(
            arrays["reconstructed_column_density_m2"], dtype=float
        )[fluence_index, :, 0]
        fit_success = np.asarray(arrays["fit_success"], dtype=bool)[
            fluence_index, :, 0
        ]

    pair_columns = 2
    row_count = int(np.ceil(len(condition_ids) / pair_columns))
    figure, axes = plt.subplots(
        row_count,
        2 * pair_columns,
        figsize=(14.0, max(9.0, 2.45 * row_count)),
        constrained_layout=True,
        squeeze=False,
    )
    colour_map = plt.get_cmap("viridis").with_extremes(bad="#e6e6e6")
    for condition_index, identifier in enumerate(condition_ids):
        row_index = condition_index // pair_columns
        pair_index = condition_index % pair_columns
        truth_axis = axes[row_index, 2 * pair_index]
        reconstruction_axis = axes[row_index, 2 * pair_index + 1]
        truth = np.where(
            cropped_support,
            truth_maps[condition_index, row_slice, column_slice] / 1e14,
            np.nan,
        )
        reconstruction = np.where(
            cropped_support,
            reconstruction_maps[condition_index, row_slice, column_slice] / 1e14,
            np.nan,
        )
        shared_upper = _shared_density_limit(truth, reconstruction)
        truth_image = truth_axis.imshow(
            np.ma.masked_invalid(truth),
            origin="lower",
            extent=extent,
            cmap=colour_map,
            vmin=0.0,
            vmax=shared_upper,
            aspect="equal",
            rasterized=True,
        )
        reconstruction_axis.imshow(
            np.ma.masked_invalid(reconstruction),
            origin="lower",
            extent=extent,
            cmap=colour_map,
            vmin=0.0,
            vmax=shared_upper,
            aspect="equal",
            rasterized=True,
        )
        truth_axis.set_title(f"{_condition_label(identifier)}\nTruth", fontsize=8.4)
        reconstruction_axis.set_title("Reconstruction", fontsize=8.4)
        if not bool(fit_success[condition_index]):
            reconstruction_axis.text(
                0.5,
                0.5,
                "fit failed",
                transform=reconstruction_axis.transAxes,
                ha="center",
                va="center",
                color="#B22222",
                fontsize=9,
                fontweight="bold",
            )
        for axis in (truth_axis, reconstruction_axis):
            axis.set_xlabel(r"$y$ ($\mu$m)", fontsize=8)
            axis.tick_params(labelsize=7)
        truth_axis.set_ylabel(r"$z$ ($\mu$m)", fontsize=8)
        colour_bar = figure.colorbar(
            truth_image,
            ax=(truth_axis, reconstruction_axis),
            fraction=0.035,
            pad=0.02,
        )
        colour_bar.set_label(r"$n_{\rm col}$ ($10^{14}$ m$^{-2}$)", fontsize=7.5)
        colour_bar.ax.tick_params(labelsize=7)

    for unused_index in range(len(condition_ids), row_count * pair_columns):
        row_index = unused_index // pair_columns
        pair_index = unused_index % pair_columns
        axes[row_index, 2 * pair_index].set_visible(False)
        axes[row_index, 2 * pair_index + 1].set_visible(False)
    figure.suptitle(
        "DPFI truth and independent reconstruction at "
        rf"$F={representative_fluence:g}$ mW $\mu$s"
        "\n"
        "fixed observable support; one shared colour scale per condition",
        fontsize=12,
    )
    return _save_figure(
        figure,
        verified.directory,
        str(figure_config["overview_basename"]),
        [str(value) for value in figure_config["formats"]],
    )


_SCORE_PANELS = (
    ("c_A", r"Integrated response $c_A$"),
    ("c_r", r"Centroid position $c_r$"),
    ("c_w", r"Major rms width $c_w$"),
)


def _plot_metric_scores(
    verified: VerifiedInitialConditionArtifactSet,
) -> list[Path]:
    condition_ids, fluences, realizations = _suite_axes(
        verified.summary,
        verified.config,
    )
    if len(realizations) != 1:
        raise ValueError("version-1 metric figure requires one realization per condition")
    rows = {
        (
            float(row["fluence_mw_us"]),
            str(row["condition_id"]),
            int(row["realization_index"]),
        ): row
        for row in verified.trial_rows
    }
    figure, axes = plt.subplots(
        len(_SCORE_PANELS),
        1,
        figsize=(14.0, 8.5),
        sharex=True,
        constrained_layout=True,
    )
    x_positions = np.arange(len(condition_ids), dtype=float)
    offsets = np.linspace(-0.22, 0.22, len(fluences))
    colours = plt.get_cmap("viridis")(np.linspace(0.15, 0.85, len(fluences)))
    for axis, (score_name, ylabel) in zip(axes, _SCORE_PANELS, strict=True):
        for offset, colour, fluence in zip(offsets, colours, fluences, strict=True):
            values = np.asarray(
                [
                    float(rows[(fluence, identifier, realizations[0])][score_name])
                    if rows[(fluence, identifier, realizations[0])][score_name] != ""
                    else np.nan
                    for identifier in condition_ids
                ],
                dtype=float,
            )
            finite = np.isfinite(values)
            axis.scatter(
                x_positions[finite] + offset,
                values[finite],
                color=colour,
                edgecolor="white",
                linewidth=0.5,
                s=34,
                label=rf"$F={fluence:g}$ mW $\mu$s",
                zorder=3,
            )
        axis.axhline(
            1.0,
            color="#B22222",
            linestyle="--",
            linewidth=1.2,
            label="usability threshold" if score_name == "c_A" else None,
            zorder=2,
        )
        axis.set_ylabel(ylabel)
        axis.set_ylim(bottom=0.0)
        axis.grid(axis="y", alpha=0.22, linewidth=0.7)
    axes[0].legend(frameon=False, ncol=len(fluences) + 1, loc="upper left")
    axes[-1].set_xticks(x_positions)
    axes[-1].set_xticklabels(
        [_condition_label(identifier) for identifier in condition_ids],
        rotation=28,
        ha="right",
        fontsize=8,
    )
    axes[-1].set_xlabel("Initial condition")
    figure.suptitle(
        "Observable-specific DPFI reconstruction usability\n"
        r"each coefficient is judged independently against $c=1$",
        fontsize=12,
    )
    figure_config = verified.config["figure"]
    return _save_figure(
        figure,
        verified.directory,
        str(figure_config["metric_basename"]),
        [str(value) for value in figure_config["formats"]],
    )


def _preferred_path(paths: list[Path]) -> Path:
    return next((path for path in paths if path.suffix.lower() == ".png"), paths[0])


def generate_initial_condition_suite_figures(
    output_directory: Path,
) -> dict[str, Path]:
    """Verify the sealed suite, render its two declared figures and hash them."""

    verified = verify_initial_condition_artifact_set(output_directory)
    overview_paths = _plot_condition_overview(verified)
    metric_paths = _plot_metric_scores(verified)
    figure_paths = [*overview_paths, *metric_paths]
    metadata_path = verified.directory / FIGURE_METADATA_FILENAME
    figure_config = verified.config["figure"]
    write_json(
        metadata_path,
        {
            "source_run_id": verified.manifest["run_id"],
            "source_manifest": MANIFEST_FILENAME,
            "source_manifest_sha256": file_sha256(
                verified.directory / MANIFEST_FILENAME
            ),
            "source_artifact_sha256": {
                filename: record["sha256"]
                for filename, record in verified.manifest["artifacts"].items()
            },
            "source_benchmark": verified.summary["source_benchmark"],
            "source_trial_table": TRIALS_FILENAME,
            "source_maps": MAPS_FILENAME,
            "representative_fluence_mw_us": float(
                figure_config["representative_fluence_mw_us"]
            ),
            "overview_colour_scale": "shared by truth and reconstruction within each condition",
            "metric_coefficients": [name for name, _ in _SCORE_PANELS],
            "metric_threshold": 1.0,
            "scores_are_reported_independently": True,
            "renderer": Path(__file__).name,
            "renderer_sha256": file_sha256(Path(__file__)),
            "claims_boundary": verified.summary["claims_boundary"],
            "output_files": [path.name for path in figure_paths],
            "rendered_output_sha256": {
                path.name: file_sha256(path) for path in figure_paths
            },
        },
    )
    return {
        "overview": _preferred_path(overview_paths),
        "metrics": _preferred_path(metric_paths),
        "metadata": metadata_path,
    }


__all__ = [
    "FIGURE_METADATA_FILENAME",
    "VerifiedInitialConditionArtifactSet",
    "generate_initial_condition_suite_figures",
    "verify_initial_condition_artifact_set",
]
