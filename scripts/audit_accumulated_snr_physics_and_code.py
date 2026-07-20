"""Physics-first and code-traceability audit for accumulated SNR and N_max."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import brentq
from scipy.special import zeta

from non_destructive_image import (
    build_thomas_fermi_state,
    reabsorption_fraction,
    scalar_phase_shift,
    scattered_photons_per_atom,
)
from scripts.recover_notebook_multishot_stage import (
    _basic_constants,
    build_blur_axis,
    pci_snr_pixel_for_phase,
    self_consistent_total_atoms,
    tf_state_for_atoms,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = REPO_ROOT / "results" / "accumulated_snr_full_physics_audit"


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _git_value(*args: str) -> str:
    commands = [["git", *args], [r"C:\Program Files\Git\cmd\git.exe", *args]]
    for command in commands:
        try:
            return subprocess.check_output(command, cwd=REPO_ROOT, text=True, stderr=subprocess.DEVNULL).strip()
        except (OSError, subprocess.CalledProcessError):
            continue
    return "unknown"


def _log_slope(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.polyfit(np.log(np.asarray(x)), np.log(np.asarray(y)), 1)[0])


def _relative_error(left: float, right: float) -> float:
    return abs(left - right) / max(abs(right), np.finfo(float).tiny)


def _parameters(notebook: dict[str, Any], plot: dict[str, Any]) -> dict[str, Any]:
    constants = notebook["constants"]
    atom = notebook["atom"]
    condensate = notebook["condensate"]
    geometry = notebook["imaging_geometry"]
    cfg = plot["accumulated_snr_invariance"]
    hbar = float(constants["hbar"])
    h_planck = 2 * np.pi * hbar
    mass = float(atom["mass_number"]) * float(constants["atomic_mass_unit"])
    scattering_length = float(condensate["scattering_length_bohr"]) * float(constants["bohr_radius_m"])
    state = build_thomas_fermi_state(
        float(condensate["atom_number"]),
        scattering_length,
        condensate["trap_frequencies_hz"],
        mass,
        hbar,
        float(constants["boltzmann_constant"]),
    )
    isat = (
        np.pi * h_planck * float(constants["speed_of_light"]) * float(atom["natural_linewidth_rad_s"])
        / (3 * float(atom["transition_wavelength_m"]) ** 3)
    )
    power_mw = float(cfg["probe_power_mw"])
    tau_s = float(cfg["exposure_time_us"]) * 1e-6
    probe_diameter = float(geometry["probe_diameter_m"])
    intensity = 2 * power_mw * 1e-3 / (np.pi * (probe_diameter / 2) ** 2)
    magnification = float(geometry["focal_length_2_m"]) / float(geometry["focal_length_1_m"])
    object_pixel = float(geometry["camera_pixel_m"]) / magnification
    photon_energy = h_planck * float(constants["speed_of_light"]) / float(atom["transition_wavelength_m"])
    nphot = intensity * object_pixel**2 * tau_s * float(cfg["quantum_efficiency"]) / photon_energy
    return {
        "hbar": hbar,
        "h_planck": h_planck,
        "mass": mass,
        "scattering_length": scattering_length,
        "tf_state": state,
        "saturation_intensity": float(isat),
        "saturation_parameter": float(intensity / isat),
        "probe_intensity_w_m2": float(intensity),
        "detected_photons_per_pixel": float(nphot),
        "power_mw": power_mw,
        "tau_s": tau_s,
        "probe_diameter_m": probe_diameter,
        "axis": int(cfg["imaging_axis"]),
        "loss_fraction": float(cfg["destruction_budget_fraction"]),
        "recoil_energy_multiplier": float(cfg["recoil_energy_multiplier"]),
        "recoil_energy_convention": str(cfg["recoil_energy_convention"]),
        "eta_coll": float(notebook["multishot_recovery"]["eta_coll"]),
        "read_noise_e": float(cfg["read_noise_electrons"]),
        "t_p": float(cfg["phase_plate_amplitude_transmittance"]),
        "gamma": float(atom["natural_linewidth_rad_s"]),
        "sigma0": float(atom["resonant_cross_section_m2"]),
        "wavelength": float(atom["transition_wavelength_m"]),
        "column_densities": np.asarray(state.column_density, dtype=float),
        "column_density_peak": float(state.column_density[int(cfg["imaging_axis"])]),
    }


def _point_physics(detuning_hz: float, notebook: dict[str, Any], params: dict[str, Any]) -> dict[str, float]:
    gamma = params["gamma"]
    delta = 2 * detuning_hz * 2 * np.pi / gamma
    ng_helper = scattered_photons_per_atom(
        detuning_hz,
        params["power_mw"],
        params["tau_s"],
        params["saturation_intensity"],
        gamma,
        params["probe_diameter_m"],
        use_peak_intensity=True,
    )
    ng_independent = (
        gamma / 2 * params["saturation_parameter"]
        / (1 + params["saturation_parameter"] + delta**2) * params["tau_s"]
    )
    phase = scalar_phase_shift(
        detuning_hz,
        params["column_density_peak"],
        params["sigma0"],
        gamma,
    )
    reabs = reabsorption_fraction(
        detuning_hz,
        params["column_densities"],
        params["sigma0"],
        gamma,
    )
    f = params["loss_fraction"]
    eta = params["eta_coll"]
    n_clean_cont = -np.log1p(-f) / (eta * ng_helper)
    p_linear = eta * ng_helper
    n_linear_multiplicative = np.log1p(-f) / np.log1p(-p_linear) if p_linear < 1 else np.nan

    constants = _basic_constants(notebook)
    thermal = self_consistent_total_atoms(notebook, constants)
    zeta3, zeta4 = float(zeta(3)), float(zeta(4))
    t0 = float(notebook["condensate"]["temperature_k"])
    tc = thermal["critical_temperature_k"]
    fc0 = 1 - (t0 / tc) ** 3
    target_t = tc * (1 - (1 - f) * fc0) ** (1 / 3)
    energy_coefficient = 3 * (zeta4 / zeta3) * float(notebook["constants"]["boltzmann_constant"]) / tc**3
    deposited_energy = (
        params["recoil_energy_multiplier"]
        * ng_helper
        * (1 + reabs)
        * constants["e_rec"]
    )
    n_heating_cont = energy_coefficient * (target_t**4 - t0**4) / deposited_energy

    nphot = params["detected_photons_per_pixel"]
    read = params["read_noise_e"]
    tp = params["t_p"]
    pci_signal = 2 * tp * abs(phase) * nphot
    pci_shot = pci_signal / np.sqrt(tp**2 * nphot)
    pci_read = pci_signal / np.sqrt(tp**2 * nphot + read**2)
    dgi_intensity = 4 * np.sin(phase / 2) ** 2
    dgi_counts = dgi_intensity * nphot
    dgi_shot = np.sqrt(dgi_counts)
    dgi_read = dgi_counts / np.sqrt(dgi_counts + read**2)
    return {
        "detuning_hz": float(detuning_hz),
        "detuning_ghz": float(detuning_hz / 1e9),
        "dimensionless_detuning": float(delta),
        "phase_rad": float(phase),
        "scattered_photons_per_atom": float(ng_helper),
        "scattering_independent_formula": float(ng_independent),
        "scattering_relative_difference": _relative_error(ng_helper, ng_independent),
        "reabsorption_fraction": float(reabs),
        "per_frame_exponential_loss_exponent": float(eta * ng_helper),
        "per_frame_linear_loss_probability_proxy": float(p_linear),
        "nmax_clean_continuous": float(n_clean_cont),
        "nmax_clean_integer_pulses_not_exceeding": int(np.floor(n_clean_cont)),
        "clean_sequence_first_crossing_index": int(np.ceil(n_clean_cont)),
        "clean_sequence_states_if_crossing_included": int(np.ceil(n_clean_cont)) + 1,
        "nmax_linear_multiplicative_probability": float(n_linear_multiplicative),
        "nmax_heating_reabs_continuous": float(n_heating_cont),
        "nmax_heating_integer_pulses_not_exceeding": int(np.floor(n_heating_cont)),
        "heating_sequence_first_crossing_index": int(np.ceil(n_heating_cont)),
        "pci_snr_shot_limit": float(pci_shot),
        "pci_snr_shot_plus_read": float(pci_read),
        "dgi_snr_shot_limit": float(dgi_shot),
        "dgi_snr_shot_plus_read": float(dgi_read),
        "pci_total_clean_initial_sqrt_n": float(pci_read * np.sqrt(n_clean_cont)),
        "dgi_total_clean_initial_sqrt_n": float(dgi_read * np.sqrt(n_clean_cont)),
        "pci_shot_total_clean_initial_sqrt_n": float(pci_shot * np.sqrt(n_clean_cont)),
        "dgi_shot_total_clean_initial_sqrt_n": float(dgi_shot * np.sqrt(n_clean_cont)),
    }


def _explicit_pci_accumulation(
    detuning_hz: float,
    notebook: dict[str, Any],
    params: dict[str, Any],
    point: dict[str, float],
    model: str,
    blur_axis: dict[int, float],
) -> dict[str, float]:
    constants = _basic_constants(notebook)
    thermal = self_consistent_total_atoms(notebook, constants)
    axis = params["axis"]
    n0_initial = float(notebook["condensate"]["atom_number"])
    if model == "clean_loss":
        allowed = int(point["nmax_clean_integer_pulses_not_exceeding"])
    else:
        allowed = int(point["nmax_heating_integer_pulses_not_exceeding"])
    allowed = max(allowed, 0)

    temperature = float(notebook["condensate"]["temperature_k"])
    tc = thermal["critical_temperature_k"]
    total_atoms = thermal["total_atoms"]
    zeta3, zeta4 = float(zeta(3)), float(zeta(4))
    energy_coefficient = 3 * (zeta4 / zeta3) * float(notebook["constants"]["boltzmann_constant"]) / tc**3
    reabs = point["reabsorption_fraction"]
    deposited = (
        params["recoil_energy_multiplier"]
        * point["scattered_photons_per_atom"]
        * (1 + reabs)
        * constants["e_rec"]
    )
    snr_values: list[float] = []
    for shot in range(allowed):
        if model == "clean_loss":
            n0_now = n0_initial * np.exp(-params["eta_coll"] * point["scattered_photons_per_atom"] * shot)
        else:
            n0_now = total_atoms * (1 - (temperature / tc) ** 3)
        state = tf_state_for_atoms(float(n0_now), constants)
        phase = scalar_phase_shift(
            detuning_hz,
            float(np.asarray(state["column_density"])[axis]),
            params["sigma0"],
            params["gamma"],
        )
        snr_values.append(
            pci_snr_pixel_for_phase(
                notebook,
                phase,
                params["power_mw"],
                axis,
                params["tau_s"],
                blur_axis,
            )
        )
        if model == "heating_reabsorption":
            temperature = (temperature**4 + deposited / energy_coefficient) ** 0.25
    accumulated = float(np.sqrt(np.sum(np.square(snr_values)))) if snr_values else 0.0
    initial = float(snr_values[0]) if snr_values else np.nan
    identical_estimate = float(initial * np.sqrt(allowed)) if allowed else 0.0
    return {
        "allowed_integer_frames": allowed,
        "initial_exact_pci_snr": initial,
        "final_exact_pci_snr": float(snr_values[-1]) if snr_values else np.nan,
        "explicit_rms_accumulated_pci_snr": accumulated,
        "initial_snr_times_sqrt_integer_frames": identical_estimate,
        "identical_frame_overestimate_fraction": (
            (identical_estimate - accumulated) / accumulated if accumulated else np.nan
        ),
    }


def _definition_rows(read_noise_e: float) -> list[dict[str, Any]]:
    return [
        {"quantity": "SNR_shot", "definition": "signal expectation divided by single-frame noise standard deviation", "units": "dimensionless", "scope": "per pixel per image in this figure", "time_dependence": "changes if density, phase, photon scale, or noise changes"},
        {"quantity": "SNR_total", "definition": "sqrt(sum_i SNR_i^2) for independent optimally combined frames", "units": "dimensionless", "scope": "accumulated sequence", "time_dependence": "cumulative"},
        {"quantity": "N_max", "definition": "model-dependent pulse budget before a stated destruction threshold", "units": "pulses or frames only after integer convention", "scope": "sequence", "time_dependence": "derived endpoint"},
        {"quantity": "destruction budget", "definition": "allowed fractional reduction of condensate atom number; configured as 0.30", "units": "fraction", "scope": "cloud", "time_dependence": "threshold"},
        {"quantity": "loss fraction", "definition": "1 - N0(frame)/N0(initial)", "units": "fraction", "scope": "condensate cloud", "time_dependence": "frame-dependent"},
        {"quantity": "scattered photons per atom per frame", "definition": "Gamma/2 * s/(1+s+delta^2) * tau", "units": "photons atom^-1 frame^-1", "scope": "per atom per pulse", "time_dependence": "constant in current fixed-probe model"},
        {"quantity": "heating contribution", "definition": "recoil_energy_multiplier*N_gamma*(1+reabsorption)*E_rec", "units": "J atom^-1 frame^-1", "scope": "per atom per pulse", "time_dependence": "constant because current reabsorption uses initial density"},
        {"quantity": "reabsorption contribution", "definition": "mean_i[1-exp(-OD_i)] multiplying recoil deposition", "units": "fraction", "scope": "angle-averaged cloud scalar", "time_dependence": "held constant in current sequence"},
        {"quantity": "atom-number depletion", "definition": "clean-loss N0=N0_initial exp(-eta N_gamma frame)", "units": "atoms", "scope": "condensate in clean-loss model", "time_dependence": "frame-dependent"},
        {"quantity": "condensate depletion", "definition": "heating model N0=N_total[1-(T/Tc)^3]", "units": "atoms", "scope": "condensate", "time_dependence": "frame-dependent"},
        {"quantity": "temperature increase", "definition": "T_next=(T^4+dE/A_E)^(1/4)", "units": "K", "scope": "cloud", "time_dependence": "frame-dependent in heating model"},
        {"quantity": "PCI signal amplitude", "definition": "figure uses 2*t_p*|phi| detected-photon contrast", "units": "electrons per pixel", "scope": "analytical peak pixel", "time_dependence": "falls as condensate depletes"},
        {"quantity": "DGI signal amplitude", "definition": "figure uses 4*sin^2(phi/2) detected photons", "units": "electrons per pixel", "scope": "ideal opaque-stop analytical peak pixel", "time_dependence": "falls quadratically with phase"},
        {"quantity": "noise variance", "definition": "photon-shot variance plus read-noise variance", "units": "electrons^2 pixel^-1 frame^-1", "scope": "per pixel per image", "time_dependence": "signal/background dependent"},
        {"quantity": "read-noise contribution", "definition": f"sigma_read^2 with sigma_read={read_noise_e:.6g} e- rms", "units": "electrons^2 pixel^-1 frame^-1", "scope": "one binned camera pixel read", "time_dependence": "constant and independent in model"},
        {"quantity": "photon-shot-noise contribution", "definition": "Poisson variance equal to expected detected count", "units": "electrons^2 pixel^-1 frame^-1", "scope": "per pixel per image", "time_dependence": "depends on total detected intensity"},
    ]


def _traceability_rows() -> list[dict[str, Any]]:
    return [
        {"physics_quantity": "dimensionless detuning", "intended_definition": "delta=2*Delta_Hz*2pi/Gamma", "code_file": "src/non_destructive_image/light_atom.py", "function_line": "dimensionless_detuning, lines 11-18", "actual_implementation": "exact intended formula", "match_status": "exact match", "issue": "none"},
        {"physics_quantity": "scalar phase", "intended_definition": "sigma0*n_col*delta/[2(1+delta^2)]", "code_file": "src/non_destructive_image/light_atom.py", "function_line": "scalar_phase_shift, lines 21-31", "actual_implementation": "exact lineshape", "match_status": "exact match", "issue": "none"},
        {"physics_quantity": "scattered photons", "intended_definition": "Gamma/2*s/(1+s+delta^2)*tau", "code_file": "src/non_destructive_image/light_atom.py", "function_line": "scattered_photons_per_atom, lines 60-83", "actual_implementation": "exact notebook formula with peak-intensity option", "match_status": "exact match", "issue": "none"},
        {"physics_quantity": "figure N_max", "intended_definition": "30% clean-loss continuous pulse budget", "code_file": "scripts/generate_accumulated_snr_invariance_plot.py", "function_line": "line 137", "actual_implementation": "-log(1-f)/(eta*N_gamma), fractional", "match_status": "exact match", "issue": "not an integer frame count; optimistic model"},
        {"physics_quantity": "heating N_max", "intended_definition": "condensate reaches 30% depletion by recoil heating plus reabsorption", "code_file": "notebook_sections/03_light_atom_interaction.py", "function_line": "Nmax_heating, lines 150-161", "actual_implementation": "closed-form T^4 threshold with initial-density reabsorption", "match_status": "approximation documented", "issue": "not used by current figure"},
        {"physics_quantity": "multishot stopping", "intended_definition": "stop at condensate loss threshold", "code_file": "src/non_destructive_image/multishot.py", "function_line": "simulate_multishot_sequence, lines 81-118", "actual_implementation": "records state, then breaks when loss>=limit", "match_status": "implementation ambiguity", "issue": "crossing state is included; array length is not N_max pulses"},
        {"physics_quantity": "RMS accumulated SNR", "intended_definition": "sqrt(cumulative sum SNR_i^2)", "code_file": "src/non_destructive_image/multishot.py", "function_line": "accumulate_snr, lines 12-15", "actual_implementation": "sqrt(nancumsum(snr^2))", "match_status": "exact match", "issue": "independence assumption remains physical caveat"},
        {"physics_quantity": "figure accumulated SNR", "intended_definition": "identical independent frames", "code_file": "scripts/generate_accumulated_snr_invariance_plot.py", "function_line": "lines 152-156", "actual_implementation": "initial analytical SNR times sqrt(fractional N_max)", "match_status": "approximation documented", "issue": "does not include depletion-driven SNR change"},
        {"physics_quantity": "PCI SNR", "intended_definition": "linear small-phase peak-pixel scaling", "code_file": "scripts/generate_accumulated_snr_invariance_plot.py", "function_line": "lines 139-143", "actual_implementation": "2*t_p*phi*Nph over carrier shot plus read variance", "match_status": "approximation documented", "issue": "not full transfer, blur, binning, or ROI"},
        {"physics_quantity": "notebook PCI SNR", "intended_definition": "full transfer peak-pixel estimate", "code_file": "scripts/recover_notebook_multishot_stage.py", "function_line": "pci_snr_pixel_for_phase, lines 180-189", "actual_implementation": "exact transfer with fixed axis blur and read noise", "match_status": "exact notebook recovery", "issue": "still peak-pixel, not integrated ROI"},
        {"physics_quantity": "DGI SNR", "intended_definition": "opaque-stop dark-field scaling", "code_file": "scripts/generate_accumulated_snr_invariance_plot.py", "function_line": "lines 145-150", "actual_implementation": "4*sin^2(phi/2), signal-photon shot variance, read variance", "match_status": "approximation documented", "issue": "does not include OD=4 leakage, pupil, binning, or matched ROI"},
        {"physics_quantity": "camera noise", "intended_definition": "Poisson detected counts plus Gaussian read noise", "code_file": "src/non_destructive_image/camera.py", "function_line": "add_camera_noise, lines 29-48", "actual_implementation": "explicit RNG; read noise per output pixel", "match_status": "exact match", "issue": "figure uses variance only, not Monte Carlo"},
        {"physics_quantity": "reabsorption", "intended_definition": "mean principal-axis absorption probability", "code_file": "src/non_destructive_image/light_atom.py", "function_line": "reabsorption_fraction, lines 102-117", "actual_implementation": "mean(1-exp(-OD))", "match_status": "exact match", "issue": "current sequence holds initial-density value fixed"},
        {"physics_quantity": "Faraday objective destructiveness", "intended_definition": "N_gamma*(1+reabsorption)", "code_file": "src/non_destructive_image/analysis.py", "function_line": "evaluate_faraday_operating_point, lines 11-81", "actual_implementation": "scalar proxy", "match_status": "approximation documented", "issue": "not used by current PCI/DGI figure"},
    ]


def _issue_rows(benchmarks: list[dict[str, Any]]) -> list[dict[str, Any]]:
    max_clean_over = max(row["clean_loss_identical_frame_overestimate_fraction"] for row in benchmarks)
    max_heat_over = max(row["heating_identical_frame_overestimate_fraction"] for row in benchmarks)
    return [
        {"severity": "Quantitative model mismatch", "issue": "Figure uses optimistic clean-loss N_max rather than heating plus reabsorption", "affected_result": "absolute N_max and accumulated-SNR magnitude", "numerical_conclusion_change": "yes; allowed frame count is lower in heating model", "required_action": "use explicit heating sequence for quantitative dissertation result; retain clean-loss only as labelled bound"},
        {"severity": "Implementation ambiguity", "issue": "N_max is fractional while frames are integer; sequence includes frame 0 and the first threshold-crossing state", "affected_result": "frame-count wording and possible off-by-one interpretation", "numerical_conclusion_change": "small for large N_max, material when N_max is only a few frames", "required_action": "define pulses, pre-pulse states, integer floor, and threshold-crossing index separately"},
        {"severity": "Quantitative model mismatch", "issue": "SNR_shot*sqrt(N_max) assumes identical frames although condensate phase and SNR decline", "affected_result": "accumulated-SNR magnitude", "numerical_conclusion_change": f"yes; benchmark overestimate reaches {max(max_clean_over, max_heat_over):.1%}", "required_action": "use sqrt(sum_i SNR_i^2) for sequence results"},
        {"severity": "Approximation requiring caveat", "issue": "Delta^2 N_max scaling is asymptotic, not exact", "affected_result": "invariance wording", "numerical_conclusion_change": "trend survives over audited range; exact flatness does not", "required_action": "say approximately/asymptotically proportional"},
        {"severity": "Approximation requiring caveat", "issue": "Figure SNR is analytical peak-pixel scaling, not full Fourier/camera/ROI simulation", "affected_result": "PCI/DGI quantitative comparison", "numerical_conclusion_change": "absolute values and prefactors may change", "required_action": "retain current scope caveat; use matched spatial ROI for quantitative comparison"},
        {"severity": "Approximation requiring caveat", "issue": "DGI opaque-stop model omits configured OD=4 leakage background", "affected_result": "high-detuning DGI shot-noise behaviour", "numerical_conclusion_change": "can matter when signal approaches leakage background", "required_action": "include leakage and total detected background in future full-pipeline calculation"},
        {"severity": "No issue", "issue": "PCI/DGI ideal prefactor near two", "affected_result": "vertical separation of ideal curves", "numerical_conclusion_change": "no", "required_action": "describe as observable-dependent homodyne prefactor, not universal efficiency"},
        {"severity": "Approximation requiring caveat", "issue": "Independent-frame RMS accumulation excludes correlated technical noise and drift", "affected_result": "long-sequence accumulated SNR", "numerical_conclusion_change": "unknown without experiment", "required_action": "add covariance/noise averaging after calibration"},
    ]


def run_audit(config_path: Path) -> dict[str, Path]:
    plot = _load_json(config_path)
    notebook_path = REPO_ROOT / plot["notebook_defaults_config"]
    notebook = _load_json(notebook_path)
    params = _parameters(notebook, plot)
    constants = _basic_constants(notebook)
    blur_axis = build_blur_axis(notebook, constants)

    benchmark_detunings = [0.75e9, 1.5e9, 3.0e9]
    points = [_point_physics(detuning, notebook, params) for detuning in benchmark_detunings]
    benchmarks: list[dict[str, Any]] = []
    nmax_rows: list[dict[str, Any]] = []
    snr_rows: list[dict[str, Any]] = []
    for point in points:
        clean = _explicit_pci_accumulation(point["detuning_hz"], notebook, params, point, "clean_loss", blur_axis)
        heating = _explicit_pci_accumulation(point["detuning_hz"], notebook, params, point, "heating_reabsorption", blur_axis)
        benchmark = dict(point)
        benchmark.update({f"clean_loss_{key}": value for key, value in clean.items()})
        benchmark.update({f"heating_{key}": value for key, value in heating.items()})
        benchmarks.append(benchmark)
        for model, value, integer_value in [
            ("clean_loss_exponential_continuous", point["nmax_clean_continuous"], point["nmax_clean_integer_pulses_not_exceeding"]),
            ("linear_probability_multiplicative", point["nmax_linear_multiplicative_probability"], int(np.floor(point["nmax_linear_multiplicative_probability"]))),
            ("heating_plus_reabsorption_continuous", point["nmax_heating_reabs_continuous"], point["nmax_heating_integer_pulses_not_exceeding"]),
        ]:
            nmax_rows.append({
                "detuning_ghz": point["detuning_ghz"], "model": model,
                "continuous_threshold_pulses": value, "integer_pulses_not_exceeding_budget": integer_value,
                "threshold": "30% condensate depletion",
            })
        for mode, noise, shot_snr, total in [
            ("PCI", "shot_noise_only", point["pci_snr_shot_limit"], point["pci_shot_total_clean_initial_sqrt_n"]),
            ("PCI", "shot_plus_read_noise", point["pci_snr_shot_plus_read"], point["pci_total_clean_initial_sqrt_n"]),
            ("DGI", "shot_noise_only", point["dgi_snr_shot_limit"], point["dgi_shot_total_clean_initial_sqrt_n"]),
            ("DGI", "shot_plus_read_noise", point["dgi_snr_shot_plus_read"], point["dgi_total_clean_initial_sqrt_n"]),
        ]:
            snr_rows.append({
                "detuning_ghz": point["detuning_ghz"], "mode": mode, "noise_model": noise,
                "per_frame_snr_initial": shot_snr, "nmax_model": "clean_loss_continuous",
                "nmax": point["nmax_clean_continuous"], "accumulated_initial_sqrt_n": total,
            })

    scan = np.geomspace(0.75, 5.0, 300)
    scan_points = [_point_physics(value * 1e9, notebook, params) for value in scan]
    ng = np.asarray([row["scattered_photons_per_atom"] for row in scan_points])
    clean_n = np.asarray([row["nmax_clean_continuous"] for row in scan_points])
    heating_n = np.asarray([row["nmax_heating_reabs_continuous"] for row in scan_points])
    ranges = {
        "0.75_to_3.0_GHz": (0.75, 3.0),
        "full_0.75_to_5.0_GHz": (0.75, 5.0),
        "low_0.75_to_1.5_GHz": (0.75, 1.5),
        "high_1.5_to_5.0_GHz": (1.5, 5.0),
    }
    exponent_rows = []
    for label, (low, high) in ranges.items():
        mask = (scan >= low) & (scan <= high)
        exponent_rows.append({
            "range": label,
            "scattered_photons_exponent": _log_slope(scan[mask], ng[mask]),
            "clean_loss_nmax_exponent": _log_slope(scan[mask], clean_n[mask]),
            "heating_reabs_nmax_exponent": _log_slope(scan[mask], heating_n[mask]),
        })

    issues = _issue_rows(benchmarks)
    summary = {
        "verdict": "B. Current figure valid only as an idealised clean-loss analytical model",
        "current_nmax_definition": "fractional continuous pulse index -log(1-0.30)/(eta_coll*N_gamma) from the optimistic clean-loss model",
        "recommended_dissertation_nmax": "integer stopping budget from the explicit heating-plus-reabsorption condensate-depletion sequence, with clean-loss retained as an optimistic comparison bound",
        "sqrt_nmax_assessment": "valid only for independent identical frames; the migrated multishot model requires sqrt(sum_i SNR_i^2)",
        "delta_squared_assessment": "asymptotic and numerically close over the audited range; not exact, and the full model changes the coefficient and slightly changes the exponent",
        "mode_independence": "physically justified for the same probe before the downstream PCI/DGI optics; detection mode does not alter upstream scattering",
        "pci_dgi_ideal_prefactor": "expected factor near two from current peak-pixel observable definitions and PCI homodyne carrier reference; not a coding normalisation bug",
        "full_multishot_invariance": "not exact in magnitude because phase and per-frame SNR decline; approximate asymptotic detuning scaling can survive self-similar evolution, but must be evaluated by explicit RMS accumulation",
        "explicit_multishot_benchmark": {
            "clean_loss_pci_accumulated_snr_0.75_to_3_GHz_ratio": (
                benchmarks[-1]["clean_loss_explicit_rms_accumulated_pci_snr"]
                / benchmarks[0]["clean_loss_explicit_rms_accumulated_pci_snr"]
            ),
            "heating_reabs_pci_accumulated_snr_0.75_to_3_GHz_ratio": (
                benchmarks[-1]["heating_explicit_rms_accumulated_pci_snr"]
                / benchmarks[0]["heating_explicit_rms_accumulated_pci_snr"]
            ),
            "identical_frame_overestimate_fraction_range": [
                min(
                    min(row["clean_loss_identical_frame_overestimate_fraction"] for row in benchmarks),
                    min(row["heating_identical_frame_overestimate_fraction"] for row in benchmarks),
                ),
                max(
                    max(row["clean_loss_identical_frame_overestimate_fraction"] for row in benchmarks),
                    max(row["heating_identical_frame_overestimate_fraction"] for row in benchmarks),
                ),
            ],
        },
        "scaling_exponents": exponent_rows,
        "parameters": {
            "config": str(config_path.relative_to(REPO_ROOT)),
            "notebook_defaults": str(notebook_path.relative_to(REPO_ROOT)),
            "power_mw": params["power_mw"],
            "exposure_time_us": params["tau_s"] * 1e6,
            "loss_fraction": params["loss_fraction"],
            "eta_coll": params["eta_coll"],
            "read_noise_e": params["read_noise_e"],
            "quantum_efficiency": plot["accumulated_snr_invariance"]["quantum_efficiency"],
            "saturation_parameter": params["saturation_parameter"],
        },
        "critical_or_quantitative_issues": [row for row in issues if row["severity"] in {"Critical physics error", "Quantitative model mismatch"}],
        "benchmark_count": len(benchmarks),
    }
    metadata = {
        "audit_type": "physics-first and code-to-physics traceability audit",
        "branch": _git_value("branch", "--show-current"),
        "git_commit": _git_value("rev-parse", "HEAD"),
        "script": str(Path(__file__).resolve().relative_to(REPO_ROOT)),
        "config": str(config_path.relative_to(REPO_ROOT)),
        "parameter_sources": [
            str(notebook_path.relative_to(REPO_ROOT)),
            "configs/dissertation_plots_v1.json",
            "results/notebook_aligned_recovery/parameter_inventory.csv",
            "results/notebook_aligned_recovery/unit_inventory.csv",
        ],
        "calibration_status": "Version 1 uncalibrated; no experimental absorption/RAI calibration",
        "figure_modified": False,
        "existing_figure_data_modified": False,
    }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    outputs = {
        "physics_definition_table": OUTPUT_DIR / "physics_definition_table.csv",
        "nmax_model_comparison": OUTPUT_DIR / "nmax_model_comparison.csv",
        "snr_model_comparison": OUTPUT_DIR / "snr_model_comparison.csv",
        "code_physics_traceability": OUTPUT_DIR / "code_physics_traceability.csv",
        "benchmark_cases": OUTPUT_DIR / "benchmark_cases.csv",
        "issue_register": OUTPUT_DIR / "issue_register.csv",
        "summary": OUTPUT_DIR / "accumulated_snr_full_audit_summary.json",
        "metadata": OUTPUT_DIR / "metadata.json",
    }
    _write_csv(outputs["physics_definition_table"], _definition_rows(float(params["read_noise_e"])))
    _write_csv(outputs["nmax_model_comparison"], nmax_rows)
    _write_csv(outputs["snr_model_comparison"], snr_rows)
    _write_csv(outputs["code_physics_traceability"], _traceability_rows())
    _write_csv(outputs["benchmark_cases"], benchmarks)
    _write_csv(outputs["issue_register"], issues)
    _write_json(outputs["summary"], summary)
    _write_json(outputs["metadata"], metadata)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("configs/dissertation_plots_v1.json"))
    args = parser.parse_args()
    config_path = args.config if args.config.is_absolute() else REPO_ROOT / args.config
    outputs = run_audit(config_path)
    for label, path in outputs.items():
        print(f"- {label}: {path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
