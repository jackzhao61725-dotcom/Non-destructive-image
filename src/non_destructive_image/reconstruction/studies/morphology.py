"""Feature-balanced calibration, freezing and held-out reconstruction study."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import asdict, dataclass, field
from time import perf_counter
from typing import Any, Literal

import numpy as np
from numpy.typing import NDArray

from ..contracts import DetectorContract, FaradayResponseContract, ReconstructionGrid
from ..density_fit import DensityFitOptions, DensityFitResult
from ..ensemble import (
    CandidateEnsembleAssessment,
    FrozenReconstructionChoice,
    ReconstructionCandidate,
    ReconstructionTrial,
    SyntheticNoisyObservation,
    assess_reconstruction_candidate,
    evaluate_frozen_candidate_on_held_out,
    generate_noisy_observation_ensemble,
    make_dark_field_candidate_initialiser,
    make_linear_candidate_initialiser,
    select_and_freeze_candidate,
)
from ..measurements import DarkFieldFaradayMeasurement, DualPortFaradayMeasurement
from ..object_models import NonnegativeBilinearDensityModel
from ..parameters import SmoothTFBounds
from ..regularisation import CurvatureAxisWeights, build_curvature_regularisation
from ..resolution import (
    CameraAlignedReduction,
    PhysicalCameraAlignedReduction,
    build_camera_aligned_reduced_grid,
    build_physical_camera_aligned_reduced_grid,
    build_uniform_reconstruction_grid,
    build_uniform_physical_camera_grid,
    faraday_grid_agreement,
)
from ..synthetic_morphologies import (
    SyntheticMorphologySplit,
    build_calibration_held_out_morphology_suite,
)
from .summaries import aggregate_held_out_trials


ReadoutName = Literal["dual_port", "dark_field"]
ProgressCallback = Callable[[str], None]


@dataclass(frozen=True)
class MorphologyStudyContext:
    """Resolved physical, numerical and synthetic inputs for one study."""

    config: dict[str, Any] = field(repr=False)
    canonical_grid: ReconstructionGrid = field(repr=False)
    grid: ReconstructionGrid = field(repr=False)
    reduction: CameraAlignedReduction | PhysicalCameraAlignedReduction
    response: FaradayResponseContract
    morphology_split: SyntheticMorphologySplit = field(repr=False)
    candidates: tuple[ReconstructionCandidate, ...] = field(repr=False)
    smooth_bounds: SmoothTFBounds = field(repr=False)
    reference_fluence_mw_us: float
    reference_count_scale_e: float
    read_noise_e: float
    jacobian_batch_size: int


@dataclass(frozen=True)
class RepresentativeReconstructionArtifact:
    """One declared held-out truth and frozen reconstruction per readout."""

    morphology_name: str
    fluence_mw_us: float
    realization_index: int
    seeds: dict[str, int]
    supported_band_errors: dict[str, float]
    truth_column_density_m2: NDArray[np.floating] = field(repr=False)
    dual_port_column_density_m2: NDArray[np.floating] = field(repr=False)
    dark_field_column_density_m2: NDArray[np.floating] = field(repr=False)
    dual_port_h_counts_e: NDArray[np.floating] = field(repr=False)
    dual_port_v_counts_e: NDArray[np.floating] = field(repr=False)
    dark_field_counts_e: NDArray[np.floating] = field(repr=False)


@dataclass(frozen=True)
class MorphologyBenchmarkRun:
    """Complete in-memory result before artifact serialization."""

    context: MorphologyStudyContext = field(repr=False)
    convergence_rows: tuple[dict[str, Any], ...]
    candidate_rows: tuple[dict[str, Any], ...]
    held_out_rows: tuple[dict[str, Any], ...]
    held_out_summary_rows: tuple[dict[str, Any], ...]
    selected: dict[str, FrozenReconstructionChoice] = field(repr=False)
    held_out_summaries: dict[str, dict[str, Any]]
    representative: RepresentativeReconstructionArtifact = field(repr=False)
    elapsed_seconds: float


@dataclass(frozen=True)
class _RepresentativeCapture:
    """Exact held-out observation and fit selected for the figure artifact."""

    seed: int
    supported_band_error: float
    truth_column_density_m2: NDArray[np.floating] = field(repr=False)
    recovered_column_density_m2: NDArray[np.floating] = field(repr=False)
    observed_channels_e: tuple[NDArray[np.floating], ...] = field(repr=False)


def _progress(callback: ProgressCallback | None, message: str) -> None:
    if callback is not None:
        callback(message)


def make_study_measurement(
    context: MorphologyStudyContext,
    readout: ReadoutName,
    *,
    fluence_mw_us: float | None = None,
    grid: ReconstructionGrid | None = None,
    count_scale_e: float | None = None,
    read_noise_e: float | None = None,
    jacobian_batch_size: int | None = None,
) -> DualPortFaradayMeasurement | DarkFieldFaradayMeasurement:
    """Build one raw-channel measurement from the resolved study contract."""

    if fluence_mw_us is not None and count_scale_e is not None:
        raise ValueError("specify fluence or count scale, not both")
    if count_scale_e is None:
        fluence = (
            context.reference_fluence_mw_us
            if fluence_mw_us is None
            else float(fluence_mw_us)
        )
        if not np.isfinite(fluence) or fluence <= 0:
            raise ValueError("fluence must be finite and positive")
        count_scale_e = (
            context.reference_count_scale_e
            * fluence
            / context.reference_fluence_mw_us
        )
    measurement_type = {
        "dual_port": DualPortFaradayMeasurement,
        "dark_field": DarkFieldFaradayMeasurement,
    }[readout]
    return measurement_type(
        grid=context.grid if grid is None else grid,
        detector=DetectorContract(
            photoelectrons_per_i0_pixel=float(count_scale_e),
            read_noise_electrons_per_pixel_per_readout=(
                context.read_noise_e if read_noise_e is None else float(read_noise_e)
            ),
        ),
        response=context.response,
        jacobian_batch_size=(
            context.jacobian_batch_size
            if jacobian_batch_size is None
            else int(jacobian_batch_size)
        ),
    )


def build_reconstruction_candidates(
    config: dict[str, Any],
    grid: ReconstructionGrid,
) -> tuple[ReconstructionCandidate, ...]:
    """Resolve every declared basis and regularisation candidate."""

    density = config["density_basis"]
    regularisation_config = density["regularisation"]
    fit = config["fit"]
    support = (
        (np.abs(grid.y_grid_m * 1e6) <= float(density["support_half_width_y_um"]))
        & (np.abs(grid.z_grid_m * 1e6) <= float(density["support_half_width_z_um"]))
    )
    fit_options = DensityFitOptions(
        irls_iterations=int(fit["irls_iterations"]),
        max_nfev=int(fit["maximum_function_evaluations"]),
        xtol=float(fit["xtol"]),
        ftol=float(fit["ftol"]),
        gtol=float(fit["gtol"]),
    )
    candidates: list[ReconstructionCandidate] = []
    for basis in density["candidates"]:
        model = NonnegativeBilinearDensityModel.from_grid(
            y_grid_m=grid.y_grid_m,
            z_grid_m=grid.z_grid_m,
            knot_y_um=basis["knot_y_um"],
            knot_z_um=basis["knot_z_um"],
            coefficient_scale_m2=float(density["coefficient_scale_m2"]),
            support_mask=support,
        )
        axis_weights = CurvatureAxisWeights(
            *regularisation_config["axis_weights"]
        )
        for weight_um2 in basis["curvature_weights_um2"]:
            weight = float(weight_um2)
            regularisation = (
                None
                if weight == 0.0
                else build_curvature_regularisation(
                    basis["knot_y_um"],
                    basis["knot_z_um"],
                    density_scale_m2=float(
                        regularisation_config["density_scale_m2"]
                    ),
                    weight_um2=weight,
                    boundary_policy=regularisation_config["boundary_policy"],
                    axis_weights=axis_weights,
                )
            )
            candidates.append(
                ReconstructionCandidate(
                    label=f"{basis['basis_label']}__curvature_{weight:g}_um2",
                    model=model,
                    coefficient_upper=float(density["coefficient_upper"]),
                    regularisation=regularisation,
                    fit_options=fit_options,
                )
            )
    labels = [candidate.label for candidate in candidates]
    if len(labels) != len(set(labels)):
        raise ValueError("reconstruction candidate labels must be unique")
    return tuple(candidates)


def build_morphology_study_context(config: dict[str, Any]) -> MorphologyStudyContext:
    """Resolve and validate a benchmark configuration without running fits."""

    grid_config = config["grid"]
    physics = config["physics"]
    detector = config["detector"]
    fit = config["fit"]
    morphology = config["synthetic_morphologies"]
    response = FaradayResponseContract(
        phase_per_column_density_rad_m2=float(
            physics["phase_per_column_density_rad_m2"]
        ),
        kappa_f=float(physics["kappa_F"]),
    )
    sampling_mode = str(grid_config.get("sampling_mode", "integer_block"))
    common_grid_arguments = {
        "canonical_ngrid": int(grid_config["canonical_ngrid"]),
        "canonical_field_of_view_m": float(
            grid_config["canonical_field_of_view_m"]
        ),
        "numerical_aperture": float(physics["numerical_aperture"]),
        "wavelength_m": float(physics["wavelength_m"]),
        "roi_half_width_y_um": float(grid_config["roi_half_width_y_um"]),
        "roi_half_width_z_um": float(grid_config["roi_half_width_z_um"]),
    }
    if sampling_mode == "integer_block":
        canonical_grid = build_uniform_reconstruction_grid(
            ngrid=common_grid_arguments["canonical_ngrid"],
            field_of_view_m=common_grid_arguments["canonical_field_of_view_m"],
            bin_size=int(grid_config["canonical_bin_size"]),
            numerical_aperture=common_grid_arguments["numerical_aperture"],
            wavelength_m=common_grid_arguments["wavelength_m"],
        )
        grid, reduction = build_camera_aligned_reduced_grid(
            canonical_bin_size=int(grid_config["canonical_bin_size"]),
            reduced_bin_size=int(grid_config["reduced_bin_size"]),
            **common_grid_arguments,
        )
    elif sampling_mode == "physical_pixel":
        camera_output_shape = tuple(int(value) for value in grid_config["camera_output_shape"])
        if len(camera_output_shape) != 2:
            raise ValueError("camera_output_shape must contain two dimensions")
        camera_pixel_size_m = float(grid_config["camera_pixel_size_m"])
        canonical_grid = build_uniform_physical_camera_grid(
            ngrid=common_grid_arguments["canonical_ngrid"],
            field_of_view_m=common_grid_arguments["canonical_field_of_view_m"],
            camera_pixel_size_m=camera_pixel_size_m,
            camera_output_shape=(camera_output_shape[0], camera_output_shape[1]),
            numerical_aperture=common_grid_arguments["numerical_aperture"],
            wavelength_m=common_grid_arguments["wavelength_m"],
        )
        grid, reduction = build_physical_camera_aligned_reduced_grid(
            camera_pixel_size_m=camera_pixel_size_m,
            camera_output_shape=(camera_output_shape[0], camera_output_shape[1]),
            reduced_ngrid=int(grid_config["reduced_ngrid"]),
            **common_grid_arguments,
        )
    else:
        raise ValueError(f"unsupported reconstruction sampling mode: {sampling_mode}")
    split = build_calibration_held_out_morphology_suite(
        grid.y_grid_m,
        grid.z_grid_m,
        peak_column_density_m2=float(morphology["peak_column_density_m2"]),
        centre_y_um=float(morphology["centre_y_um"]),
        centre_z_um=float(morphology["centre_z_um"]),
        radius_y_um=float(morphology["radius_y_um"]),
        radius_z_um=float(morphology["radius_z_um"]),
    )
    if tuple(case.name for case in split.calibration) != tuple(
        morphology["calibration_names"]
    ):
        raise RuntimeError("configured calibration names do not match the generated split")
    if tuple(case.name for case in split.held_out) != tuple(
        morphology["held_out_names"]
    ):
        raise RuntimeError("configured held-out names do not match the generated split")
    candidates = build_reconstruction_candidates(config, grid)
    return MorphologyStudyContext(
        config=config,
        canonical_grid=canonical_grid,
        grid=grid,
        reduction=reduction,
        response=response,
        morphology_split=split,
        candidates=candidates,
        smooth_bounds=SmoothTFBounds.from_mapping(fit["smooth_seed_bounds"]),
        reference_fluence_mw_us=float(detector["reference_fluence_mw_us"]),
        reference_count_scale_e=float(
            detector["photoelectrons_per_i0_pixel_at_reference"]
        ),
        read_noise_e=float(detector["read_noise_electrons_per_pixel_per_readout"]),
        jacobian_batch_size=int(fit["jacobian_batch_size"]),
    )


def _grid_convergence(context: MorphologyStudyContext) -> tuple[dict[str, Any], ...]:
    config = context.config
    morphology = config["synthetic_morphologies"]
    arguments = {
        "peak_column_density_m2": float(morphology["peak_column_density_m2"]),
        "centre_y_um": float(morphology["centre_y_um"]),
        "centre_z_um": float(morphology["centre_z_um"]),
        "radius_y_um": float(morphology["radius_y_um"]),
        "radius_z_um": float(morphology["radius_z_um"]),
    }
    canonical_split = build_calibration_held_out_morphology_suite(
        context.canonical_grid.y_grid_m,
        context.canonical_grid.z_grid_m,
        **arguments,
    )
    reduced_split = build_calibration_held_out_morphology_suite(
        context.grid.y_grid_m,
        context.grid.z_grid_m,
        **arguments,
    )
    canonical_suite = (*canonical_split.calibration, *canonical_split.held_out)
    reduced_suite = (*reduced_split.calibration, *reduced_split.held_out)
    canonical_dual = make_study_measurement(
        context,
        "dual_port",
        grid=context.canonical_grid,
        count_scale_e=1.0,
        read_noise_e=0.0,
        jacobian_batch_size=1,
    )
    reduced_dual = make_study_measurement(
        context,
        "dual_port",
        count_scale_e=1.0,
        read_noise_e=0.0,
        jacobian_batch_size=1,
    )
    canonical_dark = make_study_measurement(
        context,
        "dark_field",
        grid=context.canonical_grid,
        count_scale_e=1.0,
        read_noise_e=0.0,
        jacobian_batch_size=1,
    )
    reduced_dark = make_study_measurement(
        context,
        "dark_field",
        count_scale_e=1.0,
        read_noise_e=0.0,
        jacobian_batch_size=1,
    )
    rows: list[dict[str, Any]] = []
    for canonical_truth, reduced_truth in zip(canonical_suite, reduced_suite, strict=True):
        canonical_h, canonical_v = canonical_dual.expected_channels_from_density(
            canonical_truth.column_density_m2
        )
        reduced_h, reduced_v = reduced_dual.expected_channels_from_density(
            reduced_truth.column_density_m2
        )
        agreement = faraday_grid_agreement(
            canonical_h=canonical_h,
            canonical_v=canonical_v,
            reduced_h=reduced_h,
            reduced_v=reduced_v,
            canonical_dark=canonical_dark.expected_channels_from_density(
                canonical_truth.column_density_m2
            )[0],
            reduced_dark=reduced_dark.expected_channels_from_density(
                reduced_truth.column_density_m2
            )[0],
        )
        rows.append({"morphology": canonical_truth.name, **asdict(agreement)})
    gate = config["grid"]["convergence_gate"]
    signal_limit = float(gate["maximum_signal_relative_l2_error"])
    peak_limit = float(gate["maximum_peak_relative_error"])
    for row in rows:
        signal_errors = (
            row["dual_port_signal_relative_l2_error"],
            row["dual_port_atom_dependent_channels_relative_l2_error"],
            row["dark_field_relative_l2_error"],
        )
        peak_errors = (
            abs(row["dual_port_peak_relative_error"]),
            abs(row["dark_field_peak_relative_error"]),
        )
        if max(signal_errors) > signal_limit or max(peak_errors) > peak_limit:
            raise RuntimeError(
                f"reduced-grid convergence gate failed for {row['morphology']}: {row}"
            )
    return tuple(rows)


def _candidate_initialiser(
    context: MorphologyStudyContext,
    readout: ReadoutName,
    measurement: DualPortFaradayMeasurement | DarkFieldFaradayMeasurement,
) -> Any:
    fit = context.config["fit"]
    if readout == "dual_port":
        return make_linear_candidate_initialiser(
            measurement,
            ridge_strength=float(fit["dual_port_initialisation_ridge_strength"]),
        )
    if not isinstance(measurement, DarkFieldFaradayMeasurement):
        raise TypeError("dark-field initialisation requires a dark-field measurement")
    return make_dark_field_candidate_initialiser(
        measurement,
        smooth_bounds=context.smooth_bounds,
        projection_ridge_strength=float(fit["dark_field_projection_ridge_strength"]),
    )


def _candidate_row(
    readout: ReadoutName,
    assessment: CandidateEnsembleAssessment,
) -> dict[str, Any]:
    summary = assessment.summary
    regularisation = assessment.candidate.regularisation
    return {
        "readout": readout,
        "candidate": summary.candidate_label,
        "parameter_count": summary.parameter_count,
        "curvature_weight_um2": (
            0.0 if regularisation is None else regularisation.weight_um2
        ),
        "regularisation_density_scale_m2": (
            "" if regularisation is None else regularisation.density_scale_m2
        ),
        "regularisation_boundary_policy": (
            "none" if regularisation is None else regularisation.boundary_policy
        ),
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
    fluence_mw_us: float,
    assessment: CandidateEnsembleAssessment,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for trial in assessment.trials:
        metrics = trial.metrics
        rows.append(
            {
                "readout": readout,
                "fluence_mw_us": fluence_mw_us,
                "candidate": trial.candidate_label,
                "morphology": trial.morphology_name,
                "realization_index": trial.realization_index,
                "seed": trial.observation_seed,
                "success": trial.success,
                "data_jacobian_rank": trial.data_jacobian_rank,
                "data_jacobian_condition": trial.data_jacobian_condition,
                "parameter_count": trial.parameter_count,
                "full_map_relative_l2_error": (
                    metrics.full_map_relative_l2_error if metrics is not None else ""
                ),
                "supported_band_relative_l2_error": (
                    metrics.supported_band_relative_l2_error if metrics is not None else ""
                ),
                "integrated_density_relative_error": (
                    metrics.integrated_density_relative_error if metrics is not None else ""
                ),
                "centroid_y_error_um": metrics.centroid_y_error_um if metrics is not None else "",
                "centroid_z_error_um": metrics.centroid_z_error_um if metrics is not None else "",
                "rms_y_relative_error": metrics.rms_y_relative_error if metrics is not None else "",
                "rms_z_relative_error": metrics.rms_z_relative_error if metrics is not None else "",
                "message": trial.message,
            }
        )
    return rows


def _representative_specification(
    context: MorphologyStudyContext,
) -> tuple[str, float, int]:
    representative = context.config["representative"]
    morphology_name = str(representative["morphology_name"])
    realization_index = int(representative["realization_index"])
    representative_fluence = float(representative["fluence_mw_us"])
    held_out_names = {case.name for case in context.morphology_split.held_out}
    if morphology_name not in held_out_names:
        raise ValueError("representative morphology must belong to the held-out split")
    scanned_fluences = {
        float(value) for value in context.config["ensemble"]["held_out_fluence_mw_us"]
    }
    if representative_fluence not in scanned_fluences:
        raise ValueError("representative fluence must belong to the held-out scan")
    realization_count = int(
        context.config["ensemble"]["held_out_realizations_per_morphology"]
    )
    if not 0 <= realization_index < realization_count:
        raise ValueError(
            "representative realization_index must be within the held-out "
            f"ensemble range [0, {realization_count})"
        )
    return morphology_name, representative_fluence, realization_index


def _representative_callback(
    *,
    readout: ReadoutName,
    morphology_name: str,
    realization_index: int,
    captures: dict[str, _RepresentativeCapture],
) -> Callable[[ReconstructionTrial, SyntheticNoisyObservation, DensityFitResult], None]:
    def capture(
        trial: ReconstructionTrial,
        observation: SyntheticNoisyObservation,
        fit: DensityFitResult,
    ) -> None:
        if (
            trial.morphology_name != morphology_name
            or trial.realization_index != realization_index
        ):
            return
        if readout in captures:
            raise RuntimeError(f"duplicate representative capture for {readout}")
        if trial.metrics is None:
            raise RuntimeError("successful representative trial has no recovery metrics")
        captures[readout] = _RepresentativeCapture(
            seed=trial.observation_seed,
            supported_band_error=float(
                trial.metrics.supported_band_relative_l2_error
            ),
            truth_column_density_m2=observation.morphology.column_density_m2.copy(),
            recovered_column_density_m2=fit.column_density_m2.copy(),
            observed_channels_e=tuple(channel.copy() for channel in observation.channels),
        )

    return capture


def _build_representative_artifact(
    *,
    morphology_name: str,
    fluence_mw_us: float,
    realization_index: int,
    captures: dict[str, _RepresentativeCapture],
) -> RepresentativeReconstructionArtifact:
    expected = {"dual_port", "dark_field"}
    if set(captures) != expected:
        missing = expected.difference(captures)
        raise RuntimeError(
            "representative held-out fit was not captured for "
            + ", ".join(sorted(missing))
        )
    dual = captures["dual_port"]
    dark = captures["dark_field"]
    if not np.array_equal(dual.truth_column_density_m2, dark.truth_column_density_m2):
        raise RuntimeError("representative readouts do not share an identical truth map")
    if len(dual.observed_channels_e) != 2 or len(dark.observed_channels_e) != 1:
        raise RuntimeError("representative raw-channel count does not match the readout")
    return RepresentativeReconstructionArtifact(
        morphology_name=morphology_name,
        fluence_mw_us=fluence_mw_us,
        realization_index=realization_index,
        seeds={"dual_port": dual.seed, "dark_field": dark.seed},
        supported_band_errors={
            "dual_port": dual.supported_band_error,
            "dark_field": dark.supported_band_error,
        },
        truth_column_density_m2=dual.truth_column_density_m2,
        dual_port_column_density_m2=dual.recovered_column_density_m2,
        dark_field_column_density_m2=dark.recovered_column_density_m2,
        dual_port_h_counts_e=dual.observed_channels_e[0],
        dual_port_v_counts_e=dual.observed_channels_e[1],
        dark_field_counts_e=dark.observed_channels_e[0],
    )


def run_morphology_benchmark_study(
    config: dict[str, Any],
    *,
    progress: ProgressCallback | None = None,
) -> MorphologyBenchmarkRun:
    """Run the complete grid-gated calibration/freeze/held-out study."""

    started = perf_counter()
    context = build_morphology_study_context(config)
    representative_name, representative_fluence, representative_index = (
        _representative_specification(context)
    )
    _progress(
        progress,
        f"grid gate: canonical {context.canonical_grid.y_grid_m.shape[0]} -> "
        f"reduced {context.grid.y_grid_m.shape[0]}, camera {context.grid.camera_shape}",
    )
    convergence_rows = _grid_convergence(context)
    _progress(progress, "grid gate passed for all morphology stress tests")
    ensemble = config["ensemble"]
    candidate_rows: list[dict[str, Any]] = []
    held_out_rows: list[dict[str, Any]] = []
    selected: dict[str, FrozenReconstructionChoice] = {}
    held_out_summaries: dict[str, dict[str, Any]] = {}
    representative_captures: dict[str, _RepresentativeCapture] = {}
    for mode_index, readout in enumerate(("dual_port", "dark_field")):
        _progress(progress, f"{readout}: generating shared calibration ensemble")
        measurement = make_study_measurement(context, readout)
        calibration_observations = generate_noisy_observation_ensemble(
            measurement,
            context.morphology_split.calibration,
            realizations_per_morphology=int(
                ensemble["calibration_realizations_per_morphology"]
            ),
            base_seed=int(ensemble["calibration_seed"]) + 10000 * mode_index,
        )
        initialiser = _candidate_initialiser(context, readout, measurement)
        assessments: list[CandidateEnsembleAssessment] = []
        for candidate_index, candidate in enumerate(context.candidates, start=1):
            _progress(
                progress,
                f"{readout}: candidate {candidate_index}/{len(context.candidates)} "
                f"{candidate.label}",
            )
            assessment = assess_reconstruction_candidate(
                measurement,
                candidate,
                calibration_observations,
                initialise=initialiser,
            )
            assessments.append(assessment)
            candidate_rows.append(_candidate_row(readout, assessment))
        choice = select_and_freeze_candidate(
            assessments,
            minimum_success_fraction=float(ensemble["minimum_success_fraction"]),
            relative_error_tolerance=float(
                ensemble["near_equivalent_relative_tolerance"]
            ),
        )
        selected[readout] = choice
        _progress(progress, f"{readout}: frozen candidate {choice.candidate.label}")
        held_out_summaries[readout] = {}
        for fluence_index, fluence in enumerate(ensemble["held_out_fluence_mw_us"]):
            fluence_value = float(fluence)
            held_measurement = make_study_measurement(
                context,
                readout,
                fluence_mw_us=fluence_value,
            )
            observations = generate_noisy_observation_ensemble(
                held_measurement,
                context.morphology_split.held_out,
                realizations_per_morphology=int(
                    ensemble["held_out_realizations_per_morphology"]
                ),
                base_seed=(
                    int(ensemble["held_out_seed"])
                    + 10000 * mode_index
                    + 100 * fluence_index
                ),
            )
            _progress(
                progress,
                f"{readout}: held-out assessment at F={fluence_value:g} mW us",
            )
            held = evaluate_frozen_candidate_on_held_out(
                held_measurement,
                choice,
                observations,
                initialise=_candidate_initialiser(context, readout, held_measurement),
                on_successful_fit=(
                    _representative_callback(
                        readout=readout,
                        morphology_name=representative_name,
                        realization_index=representative_index,
                        captures=representative_captures,
                    )
                    if fluence_value == representative_fluence
                    else None
                ),
            )
            held_out_rows.extend(_trial_rows(readout, fluence_value, held.assessment))
            held_out_summaries[readout][f"{fluence_value:g}"] = asdict(
                held.assessment.summary
            )
    representative = _build_representative_artifact(
        morphology_name=representative_name,
        fluence_mw_us=representative_fluence,
        realization_index=representative_index,
        captures=representative_captures,
    )
    held_tuple = tuple(held_out_rows)
    return MorphologyBenchmarkRun(
        context=context,
        convergence_rows=convergence_rows,
        candidate_rows=tuple(candidate_rows),
        held_out_rows=held_tuple,
        held_out_summary_rows=aggregate_held_out_trials(held_tuple),
        selected=selected,
        held_out_summaries=held_out_summaries,
        representative=representative,
        elapsed_seconds=perf_counter() - started,
    )
