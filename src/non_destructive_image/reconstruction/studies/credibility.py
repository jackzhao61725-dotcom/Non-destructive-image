"""Executable credibility study for the frozen morphology reconstruction."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from time import perf_counter
from typing import Any, Literal

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

from ..credibility import (
    LocalCredibilityDiagnostics,
    ParametricBootstrapResult,
    ReconstructionStabilitySummary,
    analyse_local_credibility,
    parametric_bootstrap_reconstruction,
    summarise_density_features,
    summarise_reconstruction_stability,
)
from ..density_fit import DensityFitResult, fit_nonnegative_basis_density
from ..density_initialise import (
    dark_field_sqrt_moment_initialisation,
    linearised_nonnegative_initialisation,
)
from ..ensemble import ReconstructionCandidate
from ..evidence import NullReferenceDistribution, compare_to_zero_density
from ..measurements import (
    DarkFieldFaradayMeasurement,
    DifferentiableDensityMeasurement,
    DualPortFaradayMeasurement,
)
from ..noise import simulate_poisson_gaussian_counts
from ..regularisation import CurvatureAxisWeights, build_curvature_regularisation
from .io import file_sha256, write_json, write_rows
from .morphology import (
    MorphologyStudyContext,
    build_morphology_study_context,
    make_study_measurement,
)
from .provenance import capture_reconstruction_provenance


ReadoutName = Literal["dual_port", "dark_field"]


def _mapping_fingerprint(payload: dict[str, Any]) -> str:
    serialised = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(serialised).hexdigest()


@dataclass(frozen=True)
class CredibilityReadoutResult:
    """All credibility outputs for one raw-count readout."""

    readout: ReadoutName
    measurement: DifferentiableDensityMeasurement = field(repr=False)
    candidate: ReconstructionCandidate = field(repr=False)
    observed_channels: tuple[NDArray[np.floating], ...] = field(repr=False)
    fit: DensityFitResult = field(repr=False)
    local: LocalCredibilityDiagnostics = field(repr=False)
    bootstrap: ParametricBootstrapResult = field(repr=False)
    stability: ReconstructionStabilitySummary = field(repr=False)
    prior_variant_fits: dict[str, DensityFitResult] = field(repr=False)
    blank_rows: tuple[dict[str, Any], ...]
    pipeline_fingerprint: str
    condition_fingerprint: str


@dataclass(frozen=True)
class CredibilityStudyRun:
    """Complete in-memory credibility study before serialization."""

    context: MorphologyStudyContext = field(repr=False)
    study_config: dict[str, Any] = field(repr=False)
    truth_column_density_m2: NDArray[np.floating] = field(repr=False)
    results: dict[str, CredibilityReadoutResult] = field(repr=False)
    elapsed_seconds: float


def _candidate(
    context: MorphologyStudyContext,
    label: str,
) -> ReconstructionCandidate:
    matches = [candidate for candidate in context.candidates if candidate.label == label]
    if len(matches) != 1:
        raise ValueError(f"expected one reconstruction candidate labelled {label!r}")
    return matches[0]


def _truth_map(
    context: MorphologyStudyContext,
    name: str,
) -> NDArray[np.floating]:
    matches = [
        morphology.column_density_m2
        for morphology in context.morphology_split.held_out
        if morphology.name == name
    ]
    if len(matches) != 1:
        raise ValueError(f"expected one held-out morphology named {name!r}")
    return np.asarray(matches[0], dtype=float)


def _blind_initialisation(
    context: MorphologyStudyContext,
    readout: ReadoutName,
    measurement: DifferentiableDensityMeasurement,
    candidate: ReconstructionCandidate,
    channels: tuple[NDArray[np.floating], ...],
) -> NDArray[np.floating]:
    fit_config = context.config["fit"]
    if readout == "dual_port":
        return linearised_nonnegative_initialisation(
            measurement,
            candidate.model,
            channels,
            coefficient_upper=candidate.coefficient_upper,
            ridge_strength=float(
                fit_config["dual_port_initialisation_ridge_strength"]
            ),
        ).coefficients
    if not isinstance(measurement, DarkFieldFaradayMeasurement):
        raise TypeError("dark-field initialisation requires its measurement type")
    if len(channels) != 1:
        raise ValueError("dark-field initialisation requires one raw channel")
    return dark_field_sqrt_moment_initialisation(
        measurement,
        candidate.model,
        channels[0],
        smooth_bounds=context.smooth_bounds,
        coefficient_upper=candidate.coefficient_upper,
        projection_ridge_strength=float(
            fit_config["dark_field_projection_ridge_strength"]
        ),
    ).coefficients


def _fit_candidate(
    context: MorphologyStudyContext,
    readout: ReadoutName,
    measurement: DifferentiableDensityMeasurement,
    candidate: ReconstructionCandidate,
    channels: tuple[NDArray[np.floating], ...],
) -> DensityFitResult:
    initial = _blind_initialisation(
        context,
        readout,
        measurement,
        candidate,
        channels,
    )
    return fit_nonnegative_basis_density(
        measurement,
        candidate.model,
        channels,
        initial_coefficients=initial,
        coefficient_upper=candidate.coefficient_upper,
        regularisation=candidate.regularisation,
        options=candidate.fit_options,
    )


def _regularisation_variant(
    context: MorphologyStudyContext,
    candidate: ReconstructionCandidate,
    weight_um2: float,
) -> ReconstructionCandidate:
    regularisation_config = context.config["density_basis"]["regularisation"]
    regularisation = (
        None
        if weight_um2 == 0.0
        else build_curvature_regularisation(
            candidate.model.knot_y_um,
            candidate.model.knot_z_um,
            density_scale_m2=float(regularisation_config["density_scale_m2"]),
            weight_um2=float(weight_um2),
            boundary_policy=regularisation_config["boundary_policy"],
            axis_weights=CurvatureAxisWeights(
                *regularisation_config["axis_weights"]
            ),
        )
    )
    return ReconstructionCandidate(
        label=f"{candidate.label.split('__curvature_')[0]}__credibility_{weight_um2:g}_um2",
        model=candidate.model,
        coefficient_upper=candidate.coefficient_upper,
        regularisation=regularisation,
        fit_options=candidate.fit_options,
    )


def _prior_sensitivity(
    context: MorphologyStudyContext,
    readout: ReadoutName,
    measurement: DifferentiableDensityMeasurement,
    candidate: ReconstructionCandidate,
    channels: tuple[NDArray[np.floating], ...],
    weights: list[float],
) -> tuple[dict[str, DensityFitResult], ReconstructionStabilitySummary]:
    fits: dict[str, DensityFitResult] = {}
    for weight in weights:
        variant = _regularisation_variant(context, candidate, float(weight))
        fit = _fit_candidate(context, readout, measurement, variant, channels)
        if not fit.diagnostics.success:
            raise RuntimeError(
                f"{readout} prior-sensitivity fit failed at {weight:g} um2"
            )
        fits[f"curvature_{weight:g}_um2"] = fit
    return fits, summarise_reconstruction_stability(fits, measurement)


def _blank_ensemble(
    context: MorphologyStudyContext,
    readout: ReadoutName,
    measurement: DifferentiableDensityMeasurement,
    candidate: ReconstructionCandidate,
    *,
    draws: int,
    base_seed: int,
    reference_integrated_density: float,
    pipeline_fingerprint: str,
    condition_fingerprint: str,
) -> tuple[dict[str, Any], ...]:
    if draws <= 0:
        raise ValueError("blank draws must be positive")
    zero = np.zeros_like(context.grid.y_grid_m)
    child_sequences = np.random.SeedSequence(base_seed).spawn(draws)
    rows: list[dict[str, Any]] = []
    for index, child in enumerate(child_sequences):
        seed = int(child.generate_state(1, dtype=np.uint32)[0])
        rng = np.random.default_rng(seed)
        expected = measurement.expected_channels_from_density(zero)
        observed = tuple(
            simulate_poisson_gaussian_counts(
                channel,
                read_noise_electrons=measurement.read_noise_electrons,
                rng=rng,
            )
            for channel in expected
        )
        try:
            fit = _fit_candidate(
                context,
                readout,
                measurement,
                candidate,
                observed,
            )
            features = summarise_density_features(
                fit.column_density_m2,
                measurement,
            )
            zero_evidence = compare_to_zero_density(
                measurement,
                observed,
                fit,
                pipeline_fingerprint=pipeline_fingerprint,
                condition_fingerprint=condition_fingerprint,
                target_origin="synthetic_development",
            )
            rows.append(
                {
                    "readout": readout,
                    "draw_index": index,
                    "seed": seed,
                    "success": bool(fit.diagnostics.success),
                    "false_integrated_density": features.integrated_column_density,
                    "false_integrated_density_fraction_of_reference": (
                        features.integrated_column_density
                        / reference_integrated_density
                    ),
                    "false_peak_column_density_m2": features.peak_column_density_m2,
                    "false_peak_rotation_rad": (
                        measurement.response.rotation_per_column_density_rad_m2
                        * features.peak_column_density_m2
                    ),
                    "weighted_chi_square": fit.diagnostics.weighted_chi_square,
                    "null_quasi_deviance": zero_evidence.null_quasi_deviance,
                    "fitted_quasi_deviance": zero_evidence.fitted_quasi_deviance,
                    "quasi_deviance_improvement_over_blank": (
                        zero_evidence.delta_quasi_deviance
                    ),
                    "message": fit.diagnostics.message,
                }
            )
        except (RuntimeError, ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
            rows.append(
                {
                    "readout": readout,
                    "draw_index": index,
                    "seed": seed,
                    "success": False,
                    "false_integrated_density": "",
                    "false_integrated_density_fraction_of_reference": "",
                    "false_peak_column_density_m2": "",
                    "false_peak_rotation_rad": "",
                    "weighted_chi_square": "",
                    "null_quasi_deviance": "",
                    "fitted_quasi_deviance": "",
                    "quasi_deviance_improvement_over_blank": "",
                    "message": f"{type(error).__name__}: {error}",
                }
            )
    return tuple(rows)


def run_credibility_study(
    source_benchmark_config: dict[str, Any],
    study_config: dict[str, Any],
    *,
    progress: Any | None = None,
) -> CredibilityStudyRun:
    """Execute the declared credibility study on one held-out observation."""

    started = perf_counter()
    context = build_morphology_study_context(source_benchmark_config)
    representative = study_config["representative"]
    truth = _truth_map(context, str(representative["morphology_name"]))
    reference_measurement = make_study_measurement(
        context,
        "dual_port",
        fluence_mw_us=float(representative["fluence_mw_us"]),
    )
    reference_integrated = summarise_density_features(
        truth,
        reference_measurement,
    ).integrated_column_density
    results: dict[str, CredibilityReadoutResult] = {}
    for readout in ("dual_port", "dark_field"):
        if progress is not None:
            progress(f"{readout}: representative fit")
        measurement = make_study_measurement(
            context,
            readout,
            fluence_mw_us=float(representative["fluence_mw_us"]),
        )
        candidate = _candidate(
            context,
            str(study_config["frozen_candidates"][readout]),
        )
        seed = int(representative["seeds"][readout])
        observed = measurement.simulate_channels_from_density(
            truth,
            np.random.default_rng(seed),
        )
        fit = _fit_candidate(
            context,
            readout,
            measurement,
            candidate,
            observed,
        )
        if not fit.diagnostics.success:
            raise RuntimeError(f"{readout} representative fit failed")
        local = analyse_local_credibility(
            measurement,
            candidate.model,
            observed,
            fit,
            coefficient_upper=candidate.coefficient_upper,
            regularisation=candidate.regularisation,
        )
        if progress is not None:
            progress(f"{readout}: conditional bootstrap")
        bootstrap_config = study_config["conditional_bootstrap"]
        bootstrap = parametric_bootstrap_reconstruction(
            measurement,
            candidate.model,
            fit,
            coefficient_upper=candidate.coefficient_upper,
            regularisation=candidate.regularisation,
            draws=int(bootstrap_config["draws_per_readout"]),
            seed=int(bootstrap_config["seeds"][readout]),
            confidence_level=float(bootstrap_config["confidence_level"]),
            options=candidate.fit_options,
        )
        if progress is not None:
            progress(f"{readout}: prior sensitivity")
        prior_fits, stability = _prior_sensitivity(
            context,
            readout,
            measurement,
            candidate,
            observed,
            [
                float(value)
                for value in study_config["prior_sensitivity"][
                    "curvature_weights_um2"
                ]
            ],
        )
        if progress is not None:
            progress(f"{readout}: blank false-positive ensemble")
        blank_config = study_config["blank_ensemble"]
        pipeline_fingerprint = _mapping_fingerprint(
            {
                "source_benchmark_config": source_benchmark_config,
                "credibility_study_config": study_config,
                "candidate": candidate.label,
                "readout": readout,
            }
        )
        condition_fingerprint = _mapping_fingerprint(
            {
                "source_benchmark_config": source_benchmark_config,
                "readout": readout,
            }
        )
        blank_rows = _blank_ensemble(
            context,
            readout,
            measurement,
            candidate,
            draws=int(blank_config["draws_per_readout"]),
            base_seed=int(blank_config["seeds"][readout]),
            reference_integrated_density=reference_integrated,
            pipeline_fingerprint=pipeline_fingerprint,
            condition_fingerprint=condition_fingerprint,
        )
        results[readout] = CredibilityReadoutResult(
            readout=readout,
            measurement=measurement,
            candidate=candidate,
            observed_channels=observed,
            fit=fit,
            local=local,
            bootstrap=bootstrap,
            stability=stability,
            prior_variant_fits=prior_fits,
            blank_rows=blank_rows,
            pipeline_fingerprint=pipeline_fingerprint,
            condition_fingerprint=condition_fingerprint,
        )
    return CredibilityStudyRun(
        context=context,
        study_config=study_config,
        truth_column_density_m2=truth,
        results=results,
        elapsed_seconds=perf_counter() - started,
    )


def _feature_rows(run: CredibilityStudyRun) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for readout, result in run.results.items():
        for feature, interval in result.bootstrap.feature_intervals.items():
            rows.append(
                {
                    "readout": readout,
                    "feature": feature,
                    "estimate": interval.estimate,
                    "bootstrap_mean": interval.bootstrap_mean,
                    "bootstrap_standard_deviation": (
                        interval.bootstrap_standard_deviation
                    ),
                    "bootstrap_bias": interval.bootstrap_bias,
                    "conditional_interval_lower": interval.lower,
                    "conditional_interval_upper": interval.upper,
                    "confidence_level": result.bootstrap.confidence_level,
                    "successful_draws": result.bootstrap.successful_draws,
                    "requested_draws": result.bootstrap.requested_draws,
                }
            )
    return rows


def _mode_rows(run: CredibilityStudyRun) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for readout, result in run.results.items():
        singular = result.local.relative_data_singular_values
        fractions = result.local.generalised_data_mode_fractions
        count = max(singular.size, fractions.size)
        for index in range(count):
            rows.append(
                {
                    "readout": readout,
                    "sorted_index": index,
                    "relative_data_singular_value": (
                        float(singular[index]) if index < singular.size else ""
                    ),
                    "sorted_generalised_data_mode_fraction": (
                        float(fractions[index]) if index < fractions.size else ""
                    ),
                }
            )
    return rows


def _residual_rows(run: CredibilityStudyRun) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for readout, result in run.results.items():
        for channel in result.local.residual_channels:
            rows.append(
                {
                    "readout": readout,
                    "channel": channel.channel_name,
                    "residual_mean": channel.residual_mean,
                    "residual_rms": channel.residual_rms,
                    "residual_standard_deviation": channel.residual_standard_deviation,
                    "lag_one_correlation_y": channel.lag_one_correlation_y,
                    "lag_one_correlation_z": channel.lag_one_correlation_z,
                    "dual_port_cross_correlation": (
                        result.local.residual_port_cross_correlation
                        if result.local.residual_port_cross_correlation is not None
                        else ""
                    ),
                }
            )
    return rows


def _prior_rows(run: CredibilityStudyRun) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for readout, result in run.results.items():
        for label, fit in result.prior_variant_fits.items():
            features = asdict(
                summarise_density_features(fit.column_density_m2, result.measurement)
            )
            rows.append(
                {
                    "readout": readout,
                    "variant": label,
                    "weighted_chi_square": fit.diagnostics.weighted_chi_square,
                    "regularisation_objective": fit.diagnostics.regularisation_objective,
                    **features,
                }
            )
    return rows


def _summary(run: CredibilityStudyRun) -> dict[str, Any]:
    output: dict[str, Any] = {
        "label": run.study_config["label"],
        "purpose": run.study_config["purpose"],
        "elapsed_seconds": run.elapsed_seconds,
        "claims_boundary": run.study_config["claims_boundary"],
        "readouts": {},
    }
    for readout, result in run.results.items():
        blank_success = [row for row in result.blank_rows if row["success"]]
        blank_fractions = np.asarray(
            [
                float(row["false_integrated_density_fraction_of_reference"])
                for row in blank_success
            ],
            dtype=float,
        )
        blank_improvements = np.asarray(
            [
                float(row["quasi_deviance_improvement_over_blank"])
                for row in blank_success
            ],
            dtype=float,
        )
        blank_reference = (
            NullReferenceDistribution(
                delta_quasi_deviance=blank_improvements,
                origin="synthetic_development",
                acquisition_ids=tuple(
                    f"{readout}:blank:{row['draw_index']}" for row in blank_success
                ),
                attempted_count=len(result.blank_rows),
                failed_acquisition_ids=tuple(
                    f"{readout}:blank:{row['draw_index']}"
                    for row in result.blank_rows
                    if not row["success"]
                ),
                pipeline_fingerprint=result.pipeline_fingerprint,
                condition_fingerprint=result.condition_fingerprint,
                independent_of_target=True,
                pipeline_frozen_before_target=True,
            )
            if blank_improvements.size
            else None
        )
        representative_evidence = compare_to_zero_density(
            result.measurement,
            result.observed_channels,
            result.fit,
            pipeline_fingerprint=result.pipeline_fingerprint,
            condition_fingerprint=result.condition_fingerprint,
            target_origin="synthetic_development",
            reference=blank_reference,
        )
        representative_improvement = representative_evidence.delta_quasi_deviance
        output["readouts"][readout] = {
            "candidate": result.candidate.label,
            "free_coefficient_count": int(
                np.count_nonzero(result.local.free_coefficient_mask)
            ),
            "combined_constrained_rank": result.local.combined_constrained_rank,
            "effective_data_degrees_of_freedom": (
                result.local.effective_data_degrees_of_freedom
            ),
            "effective_prior_degrees_of_freedom": (
                result.local.effective_prior_degrees_of_freedom
            ),
            "locally_unconstrained_degrees_of_freedom": (
                result.local.locally_unconstrained_degrees_of_freedom
            ),
            "minimum_generalised_data_mode_fraction": float(
                np.min(result.local.generalised_data_mode_fractions)
            ),
            "median_generalised_data_mode_fraction": float(
                np.median(result.local.generalised_data_mode_fractions)
            ),
            "bootstrap_successful_draws": result.bootstrap.successful_draws,
            "bootstrap_requested_draws": result.bootstrap.requested_draws,
            "blank_successful_draws": len(blank_success),
            "blank_requested_draws": len(result.blank_rows),
            "blank_median_false_integrated_density_fraction_of_reference": (
                float(np.median(blank_fractions)) if blank_fractions.size else ""
            ),
            "blank_upper_95_false_integrated_density_fraction_of_reference": (
                float(np.quantile(blank_fractions, 0.95))
                if blank_fractions.size
                else ""
            ),
            "representative_quasi_deviance_improvement_over_blank": (
                representative_improvement
            ),
            "representative_evidence_level": representative_evidence.evidence_level,
            "representative_development_rank_exceedance_count": (
                representative_evidence.exceedance_count
            ),
            "representative_development_rank_reference_count": (
                representative_evidence.reference_count
            ),
            "blank_maximum_quasi_deviance_improvement": (
                float(np.max(blank_improvements))
                if blank_improvements.size
                else ""
            ),
            "development_blank_boundary": (
                "Synthetic blank ranks are a development diagnostic only; they do "
                "not define an experimental threshold or p value."
            ),
            "local_assumptions": result.local.assumptions,
            "bootstrap_assumptions": result.bootstrap.assumptions,
        }
    return output


def _plot_credibility_maps(
    run: CredibilityStudyRun,
    output_directory: Path,
) -> tuple[Path, Path]:
    extent = [
        float(run.context.grid.y_grid_m.min() * 1e6),
        float(run.context.grid.y_grid_m.max() * 1e6),
        float(run.context.grid.z_grid_m.min() * 1e6),
        float(run.context.grid.z_grid_m.max() * 1e6),
    ]
    camera_extent = [
        float(run.context.grid.camera_y_um.min()),
        float(run.context.grid.camera_y_um.max()),
        float(run.context.grid.camera_z_um.min()),
        float(run.context.grid.camera_z_um.max()),
    ]
    peak = max(
        float(np.max(result.fit.column_density_m2))
        for result in run.results.values()
    ) / 1e14
    uncertainty_peak = max(
        float(np.max(result.bootstrap.density_standard_uncertainty_m2))
        for result in run.results.values()
    ) / 1e14
    uncertainty_cmap = plt.get_cmap("magma").copy()
    uncertainty_cmap.set_bad("#d9d9d9")
    bootstrap_draws = int(
        run.study_config["conditional_bootstrap"]["draws_per_readout"]
    )
    figure, axes = plt.subplots(2, 4, figsize=(12.2, 5.8), constrained_layout=True)
    for row, readout in enumerate(("dual_port", "dark_field")):
        result = run.results[readout]
        if readout == "dual_port":
            h, v = result.observed_channels
            observable = np.divide(
                h - v,
                h + v,
                out=np.zeros_like(h),
                where=np.abs(h + v) > np.finfo(float).eps,
            )
            observable_label = "Observed $S$\n(fit uses $H,V$ counts)"
            limit = max(float(np.max(np.abs(observable))), 1e-6)
            observable_image = axes[row, 0].imshow(
                observable,
                origin="lower",
                cmap="RdBu_r",
                vmin=-limit,
                vmax=limit,
                extent=camera_extent,
                aspect="auto",
            )
            residual = (
                result.local.residual_channels[0].standardised_residual_map
                - result.local.residual_channels[1].standardised_residual_map
            ) / np.sqrt(2.0)
        else:
            observable = result.observed_channels[0]
            observable_label = "Calibrated readout ($e^-$)"
            observable_image = axes[row, 0].imshow(
                observable,
                origin="lower",
                cmap="viridis",
                extent=camera_extent,
                aspect="auto",
            )
            residual = result.local.residual_channels[0].standardised_residual_map
        axes[row, 0].set_title(observable_label)
        figure.colorbar(observable_image, ax=axes[row, 0], shrink=0.82)

        density_image = axes[row, 1].imshow(
            result.fit.column_density_m2 / 1e14,
            origin="lower",
            extent=extent,
            cmap="viridis",
            vmin=0.0,
            vmax=peak,
            aspect="auto",
        )
        axes[row, 1].set_title("Reconstructed density")
        figure.colorbar(
            density_image,
            ax=axes[row, 1],
            shrink=0.82,
            label=r"$10^{14}\,\mathrm{m}^{-2}$",
        )
        uncertainty = result.bootstrap.density_standard_uncertainty_m2 / 1e14
        uncertainty_image = axes[row, 2].imshow(
            np.ma.masked_where(uncertainty <= 0.0, uncertainty),
            origin="lower",
            extent=extent,
            cmap=uncertainty_cmap,
            vmin=0.0,
            vmax=uncertainty_peak,
            aspect="auto",
        )
        axes[row, 2].set_title(
            "Conditional bootstrap SD\n"
            f"({bootstrap_draws} draws)"
        )
        figure.colorbar(
            uncertainty_image,
            ax=axes[row, 2],
            shrink=0.82,
            label=r"$10^{14}\,\mathrm{m}^{-2}$",
        )
        residual_image = axes[row, 3].imshow(
            residual,
            origin="lower",
            cmap="RdBu_r",
            vmin=-4.0,
            vmax=4.0,
            extent=camera_extent,
            aspect="auto",
        )
        axes[row, 3].set_title("Standardised residual")
        figure.colorbar(residual_image, ax=axes[row, 3], shrink=0.82)
        readout_label = "Dual-port" if readout == "dual_port" else "Dark-field"
        axes[row, 0].set_ylabel(readout_label + "\n" + r"$z$ ($\mu$m)")
        for axis in axes[row, :]:
            axis.set_xlabel(r"$y$ ($\mu$m)")
        for axis in axes[row, 1:]:
            axis.set_ylabel(r"$z$ ($\mu$m)")
    figure.suptitle(
        r"Reconstruction credibility at $F=90$ mW $\mu$s" "\n"
        "uncertainty is conditional on the declared forward model and regulariser",
        fontsize=12,
    )
    png = output_directory / "representative_credibility_F90.png"
    pdf = output_directory / "representative_credibility_F90.pdf"
    figure.savefig(png, dpi=220)
    figure.savefig(pdf)
    plt.close(figure)
    return png, pdf


def _plot_mode_support(
    run: CredibilityStudyRun,
    output_directory: Path,
) -> tuple[Path, Path]:
    figure, axis = plt.subplots(figsize=(6.3, 4.1), constrained_layout=True)
    for readout, colour in (("dual_port", "#006BA4"), ("dark_field", "#A23B72")):
        result = run.results[readout]
        fractions = result.local.generalised_data_mode_fractions
        label = (
            f"{readout.replace('_', ' ')} "
            f"({result.local.effective_data_degrees_of_freedom:.1f} of "
            f"{np.count_nonzero(result.local.free_coefficient_mask)} free modes)"
        )
        axis.plot(
            np.arange(1, fractions.size + 1),
            fractions,
            marker="o",
            ms=2.5,
            label=label,
            color=colour,
        )
    axis.set_xlabel("Within-readout mode rank, sorted by data support")
    axis.set_ylabel("Fraction supplied by camera data (local)")
    axis.set_ylim(-0.03, 1.03)
    axis.grid(alpha=0.22)
    axis.legend(frameon=False)
    axis.set_title("Local data support of the regularised inverse")
    png = output_directory / "data_prior_mode_support_F90.png"
    pdf = output_directory / "data_prior_mode_support_F90.pdf"
    figure.savefig(png, dpi=220)
    figure.savefig(pdf)
    plt.close(figure)
    return png, pdf


def write_credibility_study_run(
    run: CredibilityStudyRun,
    study_config_path: Path,
    source_config_path: Path,
    repository_root: Path,
    *,
    provenance: dict[str, Any],
) -> dict[str, Path]:
    """Serialize one reproducible credibility run and its evidence figures."""

    output = repository_root / str(run.study_config["output_directory"])
    output.mkdir(parents=True, exist_ok=True)
    study_snapshot = output / "study_config.json"
    source_snapshot = output / "source_benchmark_config.json"
    write_json(study_snapshot, run.study_config)
    write_json(source_snapshot, run.context.config)
    summary_path = output / "credibility_summary.json"
    feature_path = output / "bootstrap_feature_intervals.csv"
    mode_path = output / "local_mode_support.csv"
    residual_path = output / "residual_diagnostics.csv"
    prior_path = output / "prior_sensitivity.csv"
    blank_path = output / "blank_false_positive_trials.csv"
    write_json(summary_path, _summary(run))
    write_rows(feature_path, _feature_rows(run))
    write_rows(mode_path, _mode_rows(run))
    write_rows(residual_path, _residual_rows(run))
    write_rows(prior_path, _prior_rows(run))
    write_rows(
        blank_path,
        [row for result in run.results.values() for row in result.blank_rows],
    )
    arrays_path = output / "credibility_arrays.npz"
    np.savez_compressed(
        arrays_path,
        truth_column_density_m2=run.truth_column_density_m2,
        dual_port_h_counts_e=run.results["dual_port"].observed_channels[0],
        dual_port_v_counts_e=run.results["dual_port"].observed_channels[1],
        dark_field_counts_e=run.results["dark_field"].observed_channels[0],
        dual_port_density_m2=run.results["dual_port"].fit.column_density_m2,
        dark_field_density_m2=run.results["dark_field"].fit.column_density_m2,
        dual_port_local_sigma_m2=run.results["dual_port"].local.density_standard_uncertainty_m2,
        dark_field_local_sigma_m2=run.results["dark_field"].local.density_standard_uncertainty_m2,
        dual_port_bootstrap_sigma_m2=run.results["dual_port"].bootstrap.density_standard_uncertainty_m2,
        dark_field_bootstrap_sigma_m2=run.results["dark_field"].bootstrap.density_standard_uncertainty_m2,
        dual_port_prior_spread_m2=run.results["dual_port"].stability.density_standard_deviation_m2,
        dark_field_prior_spread_m2=run.results["dark_field"].stability.density_standard_deviation_m2,
    )
    credibility_figures = _plot_credibility_maps(run, output)
    support_figures = _plot_mode_support(run, output)
    artifacts = (
        study_snapshot,
        source_snapshot,
        summary_path,
        feature_path,
        mode_path,
        residual_path,
        prior_path,
        blank_path,
        arrays_path,
        *credibility_figures,
        *support_figures,
    )
    metadata_path = output / "metadata.json"
    write_json(
        metadata_path,
        {
            "study_config_source": study_config_path.relative_to(repository_root).as_posix(),
            "study_config_source_sha256": file_sha256(study_config_path),
            "source_benchmark_config": source_config_path.relative_to(repository_root).as_posix(),
            "source_benchmark_config_sha256": file_sha256(source_config_path),
            "provenance": provenance,
            "artifacts_sha256": {
                path.name: file_sha256(path) for path in artifacts
            },
        },
    )
    return {
        "summary": summary_path,
        "features": feature_path,
        "modes": mode_path,
        "residuals": residual_path,
        "prior_sensitivity": prior_path,
        "blank_trials": blank_path,
        "arrays": arrays_path,
        "credibility_figure": credibility_figures[0],
        "mode_support_figure": support_figures[0],
        "metadata": metadata_path,
    }


def run_and_write_credibility_study(
    study_config: dict[str, Any],
    study_config_path: Path,
    source_config: dict[str, Any],
    source_config_path: Path,
    repository_root: Path,
    *,
    progress: Any | None = None,
) -> dict[str, Path]:
    provenance = capture_reconstruction_provenance(
        repository_root,
        entry_points=(Path("scripts/run_reconstruction_credibility_study.py"),),
    )
    run = run_credibility_study(
        source_config,
        study_config,
        progress=progress,
    )
    return write_credibility_study_run(
        run,
        study_config_path,
        source_config_path,
        repository_root,
        provenance=provenance,
    )
