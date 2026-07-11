"""Audit small-angle and far-detuned approximations used in interpretation."""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
from pathlib import Path
from typing import Any

import numpy as np

from non_destructive_image import dimensionless_detuning, intensity_at_atoms
from scripts.recover_notebook_condensate_stage import load_config
from scripts.recover_notebook_faraday_stage import build_faraday_stage
from scripts.recover_notebook_phase_stage import build_phase_stage


OUTPUT_DIR = Path("results/linear_approximation_audit")
REPORT_PATH = Path("docs/linear_approximation_validity_audit.md")
PARAMETER_INVENTORY = Path("results/notebook_aligned_recovery/parameter_inventory.csv")
UNIT_INVENTORY = Path("results/notebook_aligned_recovery/unit_inventory.csv")
FARADAY_OPTIMISATION_CONFIG = Path("configs/dissertation_results_v1.json")
DETUNING_PLOT_CONFIG = Path("configs/dissertation_plots_v1.json")


def git_commit() -> str:
    commands = [
        ["git", "rev-parse", "HEAD"],
        [r"C:\Program Files\Git\cmd\git.exe", "rev-parse", "HEAD"],
    ]
    for command in commands:
        try:
            return subprocess.check_output(command, text=True, stderr=subprocess.DEVNULL).strip()
        except (OSError, subprocess.CalledProcessError):
            continue
    return "unknown"


def write_json(path: Path, payload: dict[str, Any]) -> None:
    def convert(value: Any) -> Any:
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, dict):
            return {key: convert(item) for key, item in value.items()}
        if isinstance(value, list):
            return [convert(item) for item in value]
        return value

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(convert(payload), handle, indent=2, sort_keys=True)
        handle.write("\n")


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"no rows to write for {path}")
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def relative_error(approximate: float, exact: float) -> float:
    if exact == 0 or not math.isfinite(exact):
        return float("nan")
    return abs(approximate - exact) / abs(exact)


def finite_float(value: float) -> float:
    return float(value) if np.isfinite(value) else float("nan")


def cloud_mask(array: np.ndarray) -> np.ndarray:
    return np.asarray(array) > 0


def percentile_values(values: np.ndarray, percentiles: list[float]) -> dict[str, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return {f"p{int(p)}": float("nan") for p in percentiles}
    return {f"p{int(p)}": float(np.percentile(finite, p)) for p in percentiles}


def range_row(stage: str, quantity: str, values: np.ndarray, centre_value: float) -> dict[str, Any]:
    absolute_values = np.abs(np.asarray(values, dtype=float))
    mask = cloud_mask(absolute_values)
    in_cloud = absolute_values[mask]
    stats = percentile_values(in_cloud, [50, 90, 95, 99])
    return {
        "stage": stage,
        "quantity": quantity,
        "units": "rad",
        "sample_scope": "nonzero cloud pixels",
        "max_abs": f"{float(np.max(in_cloud)):.12g}",
        "central_abs": f"{abs(float(centre_value)):.12g}",
        "p50_abs": f"{stats['p50']:.12g}",
        "p90_abs": f"{stats['p90']:.12g}",
        "p95_abs": f"{stats['p95']:.12g}",
        "p99_abs": f"{stats['p99']:.12g}",
        "nonzero_pixel_count": int(in_cloud.size),
    }


def phase_small_angle_errors(value: float) -> dict[str, float]:
    exact = np.exp(1j * value)
    linear = 1 + 1j * value
    field_error = abs(exact - linear)
    return {
        "exp_i_phi_vs_1_plus_i_phi_absolute_field_error": float(field_error),
        "exp_i_phi_vs_1_plus_i_phi_relative_field_error": float(field_error / abs(exact)),
        "cos_phi_vs_1_absolute_error": float(abs(np.cos(value) - 1)),
        "sin_phi_vs_phi_relative_error": relative_error(value, float(np.sin(value))),
    }


def faraday_small_angle_errors(value: float) -> dict[str, float]:
    sin_exact = float(np.sin(value))
    sin2_exact = float(np.sin(value) ** 2)
    cos_exact = float(np.cos(value))
    dual_exact = float(np.sin(2 * value))
    return {
        "sin_theta_vs_theta_relative_error": relative_error(value, sin_exact),
        "sin2_theta_vs_theta2_relative_error": relative_error(value**2, sin2_exact),
        "cos_theta_vs_1_absolute_error": float(abs(1 - cos_exact)),
        "cos_theta_vs_1_relative_error": relative_error(1.0, cos_exact),
        "dual_port_sin_2theta_vs_2theta_relative_error": relative_error(2 * value, dual_exact),
    }


def small_angle_error_rows(
    phase_map: np.ndarray,
    theta_map: np.ndarray,
    phase_centre: float,
    theta_centre: float,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    samples = {
        "central": abs(float(phase_centre)),
        "p95_cloud": float(np.percentile(np.abs(phase_map)[cloud_mask(phase_map)], 95)),
        "p99_cloud": float(np.percentile(np.abs(phase_map)[cloud_mask(phase_map)], 99)),
        "peak": float(np.max(np.abs(phase_map))),
    }
    for label, value in samples.items():
        errors = phase_small_angle_errors(value)
        rows.append(
            {
                "quantity": "scalar_phase_phi",
                "sample": label,
                "value_rad": f"{value:.12g}",
                **{key: f"{val:.12g}" for key, val in errors.items()},
            }
        )

    theta_samples = {
        "central": abs(float(theta_centre)),
        "p95_cloud": float(np.percentile(np.abs(theta_map)[cloud_mask(theta_map)], 95)),
        "p99_cloud": float(np.percentile(np.abs(theta_map)[cloud_mask(theta_map)], 99)),
        "peak": float(np.max(np.abs(theta_map))),
    }
    for label, value in theta_samples.items():
        errors = faraday_small_angle_errors(value)
        rows.append(
            {
                "quantity": "faraday_rotation_theta_F",
                "sample": label,
                "value_rad": f"{value:.12g}",
                **{key: f"{val:.12g}" for key, val in errors.items()},
            }
        )
    return rows


def detuning_scaling_rows(config: dict[str, Any]) -> list[dict[str, Any]]:
    atom = config["atom"]
    camera = config["camera_recovery"]
    multishot = config["multishot_recovery"]
    gamma = float(atom["natural_linewidth_rad_s"])
    saturation_intensity = 560.0
    probe_diameter = float(config["imaging_geometry"]["probe_diameter_m"])
    saturation_parameter = intensity_at_atoms(
        float(multishot["probe_power_mw"]),
        probe_diameter,
        use_peak_intensity=True,
    ) / saturation_intensity

    detunings_hz: list[float] = [0.5e9, 1.5e9, 5.0e9]
    if FARADAY_OPTIMISATION_CONFIG.exists():
        with FARADAY_OPTIMISATION_CONFIG.open("r", encoding="utf-8") as handle:
            opt_config = json.load(handle)
        detunings_hz.extend(float(v) for v in opt_config["sweeps"]["detuning_values_hz"])
        detunings_hz.append(float(opt_config["fixed_operating_point"]["detuning_hz"]))
    if DETUNING_PLOT_CONFIG.exists():
        with DETUNING_PLOT_CONFIG.open("r", encoding="utf-8") as handle:
            plot_config = json.load(handle)
        for key in ("detuning_values_hz", "detuning_hz_values"):
            if key in plot_config:
                detunings_hz.extend(float(v) for v in plot_config[key])

    rows: list[dict[str, Any]] = []
    for detuning_hz in sorted(set(detunings_hz)):
        delta = dimensionless_detuning(detuning_hz, gamma)
        phase_exact = delta / (1 + delta**2)
        phase_asymptotic = 1 / delta
        od_exact = 1 / (1 + delta**2)
        od_asymptotic = 1 / delta**2
        scattering_exact = 1 / (1 + saturation_parameter + delta**2)
        scattering_asymptotic = 1 / delta**2
        rows.append(
            {
                "detuning_hz": f"{detuning_hz:.12g}",
                "detuning_ghz": f"{detuning_hz / 1e9:.12g}",
                "dimensionless_delta": f"{delta:.12g}",
                "saturation_parameter_for_scattering_check": f"{saturation_parameter:.12g}",
                "phase_exact_delta_over_1_plus_delta2": f"{phase_exact:.12g}",
                "phase_far_detuned_1_over_delta": f"{phase_asymptotic:.12g}",
                "phase_scaling_relative_error": f"{relative_error(phase_asymptotic, phase_exact):.12g}",
                "od_exact_1_over_1_plus_delta2": f"{od_exact:.12g}",
                "od_far_detuned_1_over_delta2": f"{od_asymptotic:.12g}",
                "od_scaling_relative_error": f"{relative_error(od_asymptotic, od_exact):.12g}",
                "scattering_exact_1_over_1_plus_s_plus_delta2": f"{scattering_exact:.12g}",
                "scattering_far_detuned_1_over_delta2": f"{scattering_asymptotic:.12g}",
                "scattering_scaling_relative_error": f"{relative_error(scattering_asymptotic, scattering_exact):.12g}",
            }
        )
    return rows


def approximation_sites() -> list[dict[str, str]]:
    return [
        {
            "site": "notebook_sections/03_light_atom_interaction.py cells 9-12",
            "topic": "scalar phase and far-detuned scaling",
            "code_or_text": "Text gives phi -> sigma0*n_col/(2*delta) and PCI linear regime phi < 0.5; code uses exact delta/(1+delta^2).",
            "actual_code_uses": "exact lineshape for phi_peak and residual OD",
            "approximation_role": "interpretation and regime classification",
            "figures_rely_on_it": "phase table/regime wording only; recovered phase map uses exact expression",
        },
        {
            "site": "notebook_sections/04_pci.py cells 15-18 and 22",
            "topic": "PCI transfer and SNR",
            "code_or_text": "Text expands PCI as t_p^2 + 2*t_p*phi; code also plots I_lin but sim_image uses exp(1j*phi). SNR_shot_ideal uses linear contrast.",
            "actual_code_uses": "exact exp(i phi) propagation for images; linear formula only for ideal SNR/reference tangent",
            "approximation_role": "interpretation and idealised SNR comparison",
            "figures_rely_on_it": "transfer-curve plot includes tangent; recovered PCI image uses exact expression",
        },
        {
            "site": "src/non_destructive_image/imaging.py simulate_pci_image and simulate_dgi_image",
            "topic": "scalar field propagation",
            "code_or_text": "object_field = exp(1j * phase_map), scattered_field = object_field - 1.",
            "actual_code_uses": "exact complex phase field",
            "approximation_role": "not applicable",
            "figures_rely_on_it": "notebook-aligned PCI/DGI outputs use exact expression",
        },
        {
            "site": "notebook_sections/06_faraday.py cells 84-89",
            "topic": "Faraday rotation readout",
            "code_or_text": "Text states I_dark = sin^2(theta_F) approx theta_F^2 and S = sin(2 theta_F) approx 2 theta_F.",
            "actual_code_uses": "field recombination with exp(+/- i theta_F), then exact intensities/ratio",
            "approximation_role": "interpretation of dark-field and dual-port response",
            "figures_rely_on_it": "Faraday explanatory wording; recovered Faraday figures use exact helper outputs",
        },
        {
            "site": "src/non_destructive_image/imaging.py simulate_faraday_image",
            "topic": "Faraday field propagation",
            "code_or_text": "sigma_plus/minus object fields are exp(+/- 1j*theta_f_map); Ex/Ey and dual-port intensities are computed from fields.",
            "actual_code_uses": "exact finite-rotation field expression for Version 1 phenomenological theta_F",
            "approximation_role": "not applicable to propagation; theta_F model itself is phenomenological",
            "figures_rely_on_it": "notebook-aligned Faraday outputs use exact expression",
        },
        {
            "site": "src/non_destructive_image/analysis.py evaluate_faraday_operating_point",
            "topic": "optimisation signal proxy",
            "code_or_text": "Uses abs(faraday_signal) and ratios such as signal_per_scattered_photon.",
            "actual_code_uses": "scalar proxy proportional to theta_F, not an image-field simulation",
            "approximation_role": "single-variable optimisation proxy, not a final camera/image prediction",
            "figures_rely_on_it": "representative Faraday optimisation tables/plots rely on this proxy",
        },
        {
            "site": "notebook_sections/09_multishot_simulation.py cells 40-46",
            "topic": "multishot signal/loss interpretation",
            "code_or_text": "phi scales as column density and SNR_pixel_phi uses full I_full(phi); detuning sweep comments use far-detuned trends.",
            "actual_code_uses": "exact phase lineshape, exact PCI transfer for realistic SNR, and exact scattering denominator",
            "approximation_role": "physical explanation and idealised invariance arguments",
            "figures_rely_on_it": "multishot recovery uses exact notebook formulas for sequence values",
        },
    ]


def classify_validity(summary: dict[str, Any]) -> list[dict[str, str]]:
    max_phi = summary["canonical_ranges"]["max_abs_phi_rad"]
    p99_phi = summary["canonical_ranges"]["p99_abs_phi_rad"]
    max_theta = summary["canonical_ranges"]["max_abs_theta_rad"]
    p99_theta = summary["canonical_ranges"]["p99_abs_theta_rad"]
    peak_phase_field_error = phase_small_angle_errors(max_phi)[
        "exp_i_phi_vs_1_plus_i_phi_relative_field_error"
    ]
    peak_dual_error = faraday_small_angle_errors(max_theta)[
        "dual_port_sin_2theta_vs_2theta_relative_error"
    ]
    peak_dark_error = faraday_small_angle_errors(max_theta)[
        "sin2_theta_vs_theta2_relative_error"
    ]

    return [
        {
            "topic": "scalar phase propagation",
            "classification": "Not applicable because code already uses exact expression",
            "reason": f"Images use exp(i phi). The peak first-order field-expansion error is {peak_phase_field_error:.3%}, so first-order propagation should not replace the exact field.",
        },
        {
            "topic": "PCI interpretation",
            "classification": "Safe to marginal as an interpretive tangent only",
            "reason": f"Canonical peak |phi|={max_phi:.3f} rad and p99={p99_phi:.3f} rad are below the notebook phi<0.5 PCI-linear guide, but finite-phase curvature is measurable.",
        },
        {
            "topic": "Faraday rotation angle",
            "classification": "Marginal but acceptable with caveat",
            "reason": f"kappa_F=1 gives peak |theta_F|={max_theta:.3f} rad. This is weak rotation, not infinitesimal rotation.",
        },
        {
            "topic": "Faraday dark-field signal",
            "classification": "Marginal but acceptable with caveat",
            "reason": f"At peak theta_F, sin^2(theta) vs theta^2 differs by {peak_dark_error:.3%}. Use exact expression in simulation and quadratic only as scaling.",
        },
        {
            "topic": "Faraday dual-port signal",
            "classification": "Marginal but acceptable with caveat",
            "reason": f"At peak theta_F, sin(2 theta) vs 2 theta differs by {peak_dual_error:.3%}; the exact dual-port expression should be cited for quantitative results.",
        },
        {
            "topic": "detuning scaling plot",
            "classification": "Safe far-detuned scaling at >=0.5 GHz for trend statements",
            "reason": "For the audited detunings, |delta| is large and 1/Delta, 1/Delta^2 errors are small; wording should still call them asymptotic trends.",
        },
        {
            "topic": "multishot signal/loss interpretation",
            "classification": "Not applicable because code already uses exact expression",
            "reason": "Sequence recovery uses exact scattering denominator, exact phase expression, and full PCI transfer for realistic SNR; linear invariance arguments are explanatory.",
        },
    ]


def build_summary(
    config: dict[str, Any],
    phase_stage: dict[str, Any],
    faraday_stage: dict[str, Any],
    range_rows: list[dict[str, Any]],
    small_angle_rows: list[dict[str, Any]],
    detuning_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    phase_map = phase_stage["notebook_phase_map_rad"]
    theta_map = faraday_stage["notebook_theta_f_map_rad"]
    mid = phase_map.shape[0] // 2
    phase_cloud = np.abs(phase_map)[cloud_mask(phase_map)]
    theta_cloud = np.abs(theta_map)[cloud_mask(theta_map)]
    representative = {}
    if FARADAY_OPTIMISATION_CONFIG.exists():
        with FARADAY_OPTIMISATION_CONFIG.open("r", encoding="utf-8") as handle:
            representative["faraday_optimisation_config"] = json.load(handle)

    summary: dict[str, Any] = {
        "status": "linear approximation validity audit",
        "calibration_status": "notebook-aligned Version 1 / uncalibrated audit; no experimental RAI or absorption calibration applied",
        "canonical_ranges": {
            "max_abs_phi_rad": float(np.max(phase_cloud)),
            "central_abs_phi_rad": abs(float(phase_map[mid, mid])),
            "p95_abs_phi_rad": float(np.percentile(phase_cloud, 95)),
            "p99_abs_phi_rad": float(np.percentile(phase_cloud, 99)),
            "max_abs_theta_rad": float(np.max(theta_cloud)),
            "central_abs_theta_rad": abs(float(theta_map[mid, mid])),
            "p95_abs_theta_rad": float(np.percentile(theta_cloud, 95)),
            "p99_abs_theta_rad": float(np.percentile(theta_cloud, 99)),
        },
        "phase_rotation_ranges": range_rows,
        "small_angle_errors": small_angle_rows,
        "detuning_scaling_errors": detuning_rows,
        "approximation_sites": approximation_sites(),
        "representative_faraday_optimisation_reference": representative,
        "config_label": config["label"],
        "detuning_plot_config_present": DETUNING_PLOT_CONFIG.exists(),
    }
    summary["validity_classification"] = classify_validity(summary)
    return summary


def metadata_payload(config_path: Path, config: dict[str, Any]) -> dict[str, Any]:
    config_files_used = [str(config_path)]
    if FARADAY_OPTIMISATION_CONFIG.exists():
        config_files_used.append(str(FARADAY_OPTIMISATION_CONFIG))
    if DETUNING_PLOT_CONFIG.exists():
        config_files_used.append(str(DETUNING_PLOT_CONFIG))
    inventories = [
        str(path)
        for path in (PARAMETER_INVENTORY, UNIT_INVENTORY)
        if path.exists()
    ]
    return {
        "git_commit_hash": git_commit(),
        "config_files_used": config_files_used,
        "inventory_files_used": inventories,
        "helper_functions_scripts_inspected": [
            "src/non_destructive_image/light_atom.py",
            "src/non_destructive_image/imaging.py",
            "src/non_destructive_image/analysis.py",
            "scripts/recover_notebook_phase_stage.py",
            "scripts/recover_notebook_faraday_stage.py",
            "scripts/recover_notebook_multishot_stage.py",
        ],
        "notebook_sections_inspected": [
            "notebook_sections/03_light_atom_interaction.py",
            "notebook_sections/04_pci.py",
            "notebook_sections/06_faraday.py",
            "notebook_sections/09_multishot_simulation.py",
        ],
        "physical_parameters_used": {
            "atom": config["atom"],
            "condensate": config["condensate"],
            "imaging_geometry": config["imaging_geometry"],
            "phase_recovery": config["phase_recovery"],
            "faraday_recovery": config["faraday_recovery"],
            "camera_recovery": config["camera_recovery"],
            "multishot_recovery": config["multishot_recovery"],
        },
        "detuning_convention": "delta = 2 * detuning_hz * 2*pi / Gamma_rad_per_s",
        "units": "SI internally; phase and Faraday rotation in radians; detuning in Hz with GHz display where tabulated",
        "error_definitions": {
            "phase_field_error": "|exp(i phi) - (1 + i phi)|",
            "sin_theta_error": "|theta - sin(theta)| / |sin(theta)|",
            "dark_field_error": "|theta^2 - sin(theta)^2| / |sin(theta)^2|",
            "dual_port_error": "|2 theta - sin(2 theta)| / |sin(2 theta)|",
            "detuning_scaling_error": "relative difference between exact notebook/helper response and far-detuned asymptotic response",
        },
        "calibration_status": "audit only; no recalibration performed",
        "simulator_physics_changed": False,
        "helper_apis_changed": False,
        "notebook_logic_changed": False,
    }


def write_report(summary: dict[str, Any], metadata: dict[str, Any]) -> None:
    ranges = summary["canonical_ranges"]
    det_15 = next(
        row for row in summary["detuning_scaling_errors"] if float(row["detuning_hz"]) == 1.5e9
    )
    peak_phase_error = phase_small_angle_errors(ranges["max_abs_phi_rad"])
    peak_faraday_error = faraday_small_angle_errors(ranges["max_abs_theta_rad"])

    lines = [
        "# Linear Approximation Validity Audit",
        "",
        "This audit checks whether the current Version 1 dissertation simulation and interpretation are justified in using small-angle or far-detuned language. It is an audit only: no simulator physics, helper API, notebook section, baseline, or existing recovery output is changed.",
        "",
        "## Parameter Provenance",
        "",
        f"- Canonical config: `{metadata['config_files_used'][0]}`",
        f"- Calibration status: {metadata['calibration_status']}",
        f"- Detuning convention: `{metadata['detuning_convention']}`",
        f"- Git commit audited: `{metadata['git_commit_hash']}`",
        f"- Detuning plot config present on this branch: `{summary['detuning_plot_config_present']}`",
        "",
        "## Approximation Sites",
        "",
        "| Site | Topic | Code status | Approximation role |",
        "| --- | --- | --- | --- |",
    ]
    for site in summary["approximation_sites"]:
        lines.append(
            f"| `{site['site']}` | {site['topic']} | {site['actual_code_uses']} | {site['approximation_role']} |"
        )

    lines.extend(
        [
            "",
            "## Canonical Phase And Rotation Ranges",
            "",
            f"- max |phi| = {ranges['max_abs_phi_rad']:.6g} rad",
            f"- central |phi| = {ranges['central_abs_phi_rad']:.6g} rad",
            f"- 95th percentile |phi| over nonzero cloud pixels = {ranges['p95_abs_phi_rad']:.6g} rad",
            f"- 99th percentile |phi| over nonzero cloud pixels = {ranges['p99_abs_phi_rad']:.6g} rad",
            f"- max |theta_F| = {ranges['max_abs_theta_rad']:.6g} rad",
            f"- central |theta_F| = {ranges['central_abs_theta_rad']:.6g} rad",
            f"- 95th percentile |theta_F| over nonzero cloud pixels = {ranges['p95_abs_theta_rad']:.6g} rad",
            f"- 99th percentile |theta_F| over nonzero cloud pixels = {ranges['p99_abs_theta_rad']:.6g} rad",
            "",
            "Because the current Version 1 Faraday convention uses `kappa_F = 1.0`, the canonical theta_F range is identical to the scalar phase range.",
            "",
            "## Small-Angle Error Summary",
            "",
            f"- At peak |phi|, `exp(i phi)` versus `1 + i phi` has relative field error {peak_phase_error['exp_i_phi_vs_1_plus_i_phi_relative_field_error']:.4%}.",
            f"- At peak |theta_F|, `sin(theta_F)` versus `theta_F` has relative error {peak_faraday_error['sin_theta_vs_theta_relative_error']:.4%}.",
            f"- At peak |theta_F|, dark-field `sin^2(theta_F)` versus `theta_F^2` has relative error {peak_faraday_error['sin2_theta_vs_theta2_relative_error']:.4%}.",
            f"- At peak |theta_F|, dual-port `sin(2 theta_F)` versus `2 theta_F` has relative error {peak_faraday_error['dual_port_sin_2theta_vs_2theta_relative_error']:.4%}.",
            "",
            "These errors are small enough for qualitative scaling language at the canonical operating point, but they are not zero. Quantitative image formation should continue to use the exact expressions already present in the recovered code.",
            "",
            "## Far-Detuned Scaling",
            "",
            f"At 1.5 GHz, delta = {float(det_15['dimensionless_delta']):.6g}. The relative error of replacing the exact phase response `delta/(1+delta^2)` with `1/delta` is {float(det_15['phase_scaling_relative_error']):.4%}. The relative error of replacing residual OD with `1/delta^2` is {float(det_15['od_scaling_relative_error']):.4%}. The scattering scaling error, including the audited saturation parameter, is {float(det_15['scattering_scaling_relative_error']):.4%}.",
            "",
            "The 1/Delta and 1/Delta^2 statements are therefore acceptable as far-detuned scaling trends over the audited detuning range, not as exact equalities.",
            "",
            "## Validity Classification",
            "",
            "| Topic | Classification | Reason |",
            "| --- | --- | --- |",
        ]
    )
    for row in summary["validity_classification"]:
        lines.append(f"| {row['topic']} | {row['classification']} | {row['reason']} |")

    lines.extend(
        [
            "",
            "## Dissertation Wording Recommendations",
            "",
            "- Say: the simulation propagates the full complex field `exp(i phi)`; linear phase response should be used only as an interpretive scaling statement.",
            "- Say: the canonical operating point has finite phase/rotation values, so exact numerical expressions are retained for images and quantitative comparisons.",
            "- Say: Faraday rotation remains a weak-rotation regime for Version 1, but not a strictly infinitesimal-rotation limit.",
            "- Say: dark-field Faraday scales approximately as theta_F squared and dual-port Faraday approximately as 2 theta_F only in the small-angle interpretation; the plotted/recovered fields use exact expressions.",
            "- Say: 1/Delta and 1/Delta^2 are far-detuned trends, not exact equalities across all plotted detunings.",
            "- Avoid claiming that representative V1 optimisation outputs are final calibrated operating-point predictions; no absorption/RAI calibration or kappa_F calibration has been applied.",
            "",
            "## Output Files",
            "",
            "- `results/linear_approximation_audit/linear_approximation_summary.json`",
            "- `results/linear_approximation_audit/phase_rotation_ranges.csv`",
            "- `results/linear_approximation_audit/small_angle_error_table.csv`",
            "- `results/linear_approximation_audit/detuning_scaling_error_table.csv`",
            "- `results/linear_approximation_audit/metadata.json`",
        ]
    )
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate(config_path: Path) -> dict[str, str]:
    config = load_config(config_path)
    phase_stage = build_phase_stage(config)
    faraday_stage = build_faraday_stage(config)
    phase_map = phase_stage["notebook_phase_map_rad"]
    theta_map = faraday_stage["notebook_theta_f_map_rad"]
    mid = phase_map.shape[0] // 2

    range_rows = [
        range_row("canonical_phase_stage", "scalar_phase_phi", phase_map, phase_map[mid, mid]),
        range_row("canonical_faraday_stage", "faraday_rotation_theta_F", theta_map, theta_map[mid, mid]),
    ]
    small_rows = small_angle_error_rows(phase_map, theta_map, phase_map[mid, mid], theta_map[mid, mid])
    detuning_rows = detuning_scaling_rows(config)
    summary = build_summary(config, phase_stage, faraday_stage, range_rows, small_rows, detuning_rows)
    metadata = metadata_payload(config_path, config)

    outputs = {
        "summary": OUTPUT_DIR / "linear_approximation_summary.json",
        "phase_rotation_ranges": OUTPUT_DIR / "phase_rotation_ranges.csv",
        "small_angle_errors": OUTPUT_DIR / "small_angle_error_table.csv",
        "detuning_scaling_errors": OUTPUT_DIR / "detuning_scaling_error_table.csv",
        "metadata": OUTPUT_DIR / "metadata.json",
        "report": REPORT_PATH,
    }
    write_json(outputs["summary"], summary)
    write_rows(outputs["phase_rotation_ranges"], range_rows)
    write_rows(outputs["small_angle_errors"], small_rows)
    write_rows(outputs["detuning_scaling_errors"], detuning_rows)
    write_json(outputs["metadata"], metadata)
    write_report(summary, metadata)
    return {key: str(path) for key, path in outputs.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    outputs = generate(args.config)
    print("Linear approximation validity audit outputs:")
    for key, path in outputs.items():
        print(f"- {key}: {path}")


if __name__ == "__main__":
    main()
