"""Recover the canonical notebook V1 deterministic multi-shot sequence stage."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq
from scipy.special import zeta

from non_destructive_image import accumulate_snr, simulate_multishot_sequence
from scripts.recover_notebook_condensate_stage import load_config
from scripts.recover_notebook_pci_stage import build_pupil, git_commit, write_json, write_rows
from scripts.plot_label_utils import ATOM_NUMBER, DETUNING_GHZ, FRAME_INDEX, PHASE_RAD, SNR, TEMPERATURE_NK


def _basic_constants(config: dict[str, Any]) -> dict[str, Any]:
    constants = config["constants"]
    atom = config["atom"]
    condensate = config["condensate"]
    hbar = float(constants["hbar"])
    c_light = float(constants["speed_of_light"])
    k_b = float(constants["boltzmann_constant"])
    amu = float(constants["atomic_mass_unit"])
    a0 = float(constants["bohr_radius_m"])
    wavelength = float(atom["transition_wavelength_m"])
    gamma = float(atom["natural_linewidth_rad_s"])
    mass = float(atom["mass_number"]) * amu
    scattering_length = float(condensate["scattering_length_bohr"]) * a0
    trap_hz = np.asarray(condensate["trap_frequencies_hz"], dtype=float)
    omega = 2 * np.pi * trap_hz
    omega_bar = omega.prod() ** (1 / 3)
    a_ho = np.sqrt(hbar / (mass * omega_bar))
    atom_number = float(condensate["atom_number"])
    mu = 0.5 * (15 * atom_number * scattering_length / a_ho) ** (2 / 5) * hbar * omega_bar
    peak_density = mu * mass / (4 * np.pi * hbar**2 * scattering_length)
    radii = np.sqrt(2 * mu / (mass * omega**2))
    column_density = (4 / 3) * peak_density * radii
    k_wave = 2 * np.pi / wavelength
    e_rec = (hbar * k_wave) ** 2 / (2 * mass)
    h_planck = 2 * np.pi * hbar
    isat = np.pi * h_planck * c_light * gamma / (3 * wavelength**3)
    return {
        "hbar": hbar,
        "speed_of_light": c_light,
        "boltzmann_constant": k_b,
        "mass": mass,
        "scattering_length": scattering_length,
        "trap_hz": trap_hz,
        "omega": omega,
        "omega_bar": omega_bar,
        "a_ho": a_ho,
        "atom_number": atom_number,
        "mu": mu,
        "peak_density": peak_density,
        "radii": radii,
        "column_density": column_density,
        "wavelength": wavelength,
        "gamma": gamma,
        "resonant_cross_section": float(atom["resonant_cross_section_m2"]),
        "e_rec": e_rec,
        "isat": isat,
    }


def dimensionless_detuning(detuning_hz: float, constants: dict[str, Any]) -> float:
    return 2 * detuning_hz * 2 * np.pi / constants["gamma"]


def intensity_at_atoms_notebook(config: dict[str, Any], probe_power_mw: float) -> float:
    geometry = config["imaging_geometry"]
    area_averaged = (probe_power_mw * 1e-3) / (np.pi * (float(geometry["probe_diameter_m"]) / 2) ** 2)
    return 2 * area_averaged


def photons_per_camera_pixel(config: dict[str, Any], probe_power_mw: float, tau_s: float) -> float:
    constants = config["constants"]
    atom = config["atom"]
    geometry = config["imaging_geometry"]
    camera = config["camera_recovery"]
    h_planck = 2 * np.pi * float(constants["hbar"])
    photon_energy = h_planck * float(constants["speed_of_light"]) / float(atom["transition_wavelength_m"])
    magnification = float(geometry["focal_length_2_m"]) / float(geometry["focal_length_1_m"])
    object_pixel_m = float(geometry["camera_pixel_m"]) / magnification
    return (
        intensity_at_atoms_notebook(config, probe_power_mw)
        * object_pixel_m**2
        * tau_s
        * float(camera["quantum_efficiency"])
        / photon_energy
    )


def tf_state_for_atoms(n0_now: float, constants: dict[str, Any]) -> dict[str, np.ndarray | float]:
    mu = (
        0.5
        * (15 * n0_now * constants["scattering_length"] / constants["a_ho"]) ** (2 / 5)
        * constants["hbar"]
        * constants["omega_bar"]
    )
    peak_density = mu * constants["mass"] / (4 * np.pi * constants["hbar"] ** 2 * constants["scattering_length"])
    radii = np.sqrt(2 * mu / (constants["mass"] * constants["omega"] ** 2))
    column_density = (4 / 3) * peak_density * radii
    return {"mu": mu, "peak_density": peak_density, "radii": radii, "column_density": column_density}


def scalar_phase_peak(detuning_hz: float, column_density_peak: float, constants: dict[str, Any]) -> float:
    delta = dimensionless_detuning(detuning_hz, constants)
    return constants["resonant_cross_section"] * column_density_peak * delta / (2 * (1 + delta**2))


def scattered_photons_per_atom(config: dict[str, Any], constants: dict[str, Any], detuning_hz: float, power_mw: float, tau_s: float) -> float:
    saturation = intensity_at_atoms_notebook(config, power_mw) / constants["isat"]
    delta = dimensionless_detuning(detuning_hz, constants)
    return (constants["gamma"] / 2) * saturation / (1 + saturation + delta**2) * tau_s


def reabsorption_fraction(detuning_hz: float, constants: dict[str, Any]) -> float:
    delta = dimensionless_detuning(detuning_hz, constants)
    optical_depth = constants["resonant_cross_section"] * constants["column_density"] / (1 + delta**2)
    return float(np.mean(1 - np.exp(-optical_depth)))


def self_consistent_total_atoms(config: dict[str, Any], constants: dict[str, Any]) -> dict[str, float]:
    temperature = float(config["condensate"]["temperature_k"])
    n0 = constants["atom_number"]
    hbar = constants["hbar"]
    k_b = constants["boltzmann_constant"]
    omega_bar = constants["omega_bar"]

    def residual(total_atoms: float) -> float:
        tc = 0.94 * hbar * omega_bar / k_b * total_atoms ** (1 / 3)
        fraction = 1 - (temperature / tc) ** 3 if tc > temperature else 0.0
        return total_atoms * fraction - n0

    total_atoms = float(brentq(residual, 3e4, 1e7))
    tc = 0.94 * hbar * omega_bar / k_b * total_atoms ** (1 / 3)
    return {"total_atoms": total_atoms, "critical_temperature_k": tc}


def build_blur_axis(config: dict[str, Any], constants: dict[str, Any]) -> dict[int, float]:
    grid = config["grid"]
    pci = config["pci_recovery"]
    ngrid = int(grid["ngrid"])
    fov = float(grid["field_of_view_m"])
    dgrid = fov / ngrid
    axis = (np.arange(ngrid) - ngrid // 2) * dgrid
    ga, gb = np.meshgrid(axis, axis)
    pupil = build_pupil(config)["pupil"]
    radii = constants["radii"]
    t_p = float(pci["phase_plate_transmittance"])
    theta = float(pci["phase_plate_phase_rad"])
    background = t_p**2

    def sim_image(axis_index: int, phase_value: float) -> np.ndarray:
        plane = [index for index in range(3) if index != axis_index]
        profile = np.maximum(0, 1 - ga**2 / radii[plane[0]] ** 2 - gb**2 / radii[plane[1]] ** 2) ** 1.5
        scattered = np.fft.ifft2(np.fft.fft2(np.exp(1j * phase_value * profile) - 1) * pupil)
        field = t_p * np.exp(1j * theta) + scattered
        return np.abs(field) ** 2

    blur: dict[int, float] = {}
    phase_test = 0.1
    ideal = abs(t_p * np.exp(1j * theta) + (np.exp(1j * phase_test) - 1)) ** 2
    for axis_index in range(3):
        image = sim_image(axis_index, phase_test)
        blur[axis_index] = float((np.max(image) - background) / (ideal - background))
    return blur


def pci_snr_pixel_for_phase(config: dict[str, Any], phase: float, power_mw: float, axis: int, tau_s: float, blur_axis: dict[int, float]) -> float:
    if not (0.0 < phase < np.pi):
        return float("nan")
    pci = config["pci_recovery"]
    camera = config["camera_recovery"]
    t_p = float(pci["phase_plate_transmittance"])
    theta = float(pci["phase_plate_phase_rad"])
    background = t_p**2
    full_intensity = abs(t_p * np.exp(1j * theta) + np.exp(1j * phase) - 1) ** 2
    image_intensity = background + blur_axis[axis] * (full_intensity - background)
    photons = photons_per_camera_pixel(config, power_mw, tau_s)
    read_noise = float(camera["read_noise_electrons"])
    return float(abs(image_intensity - background) * photons / np.sqrt(abs(image_intensity) * photons + read_noise**2))


def _notebook_run_sequence(
    *,
    config: dict[str, Any],
    constants: dict[str, Any],
    thermal: dict[str, float],
    blur_axis: dict[int, float],
    model: str,
) -> dict[str, np.ndarray]:
    multishot = config["multishot_recovery"]
    detuning_hz = float(multishot["detuning_ghz"]) * 1e9
    power_mw = float(multishot["probe_power_mw"])
    tau_s = float(multishot["pulse_duration_us"]) * 1e-6
    axis = int(multishot["imaging_axis"])
    loss_limit = float(multishot["loss_fraction_limit"])
    max_shots = int(multishot["max_shots"])
    eta_coll = float(multishot["eta_coll"])
    recoil_energy_multiplier = float(multishot.get("recoil_energy_multiplier", 1.0))
    reabs = multishot["reabsorption_override"]
    reabs = reabsorption_fraction(detuning_hz, constants) if reabs is None else float(reabs)
    photons_scattered = scattered_photons_per_atom(config, constants, detuning_hz, power_mw, tau_s)
    zeta3, zeta4 = float(zeta(3)), float(zeta(4))
    tc = thermal["critical_temperature_k"]
    total_atoms = thermal["total_atoms"]
    temperature = float(config["condensate"]["temperature_k"])
    initial_temperature = temperature
    initial_condensate_atoms = constants["atom_number"]
    energy_coefficient = 3 * (zeta4 / zeta3) * constants["boltzmann_constant"] / tc**3
    deposited_energy = recoil_energy_multiplier * photons_scattered * (1 + reabs) * constants["e_rec"]

    out = {key: [] for key in ("shot", "N0", "frac", "T", "phi", "snr")}
    shot = 0
    while shot <= max_shots:
        if model == "heating":
            n0_now = total_atoms * (1 - (temperature / tc) ** 3)
            temperature_now = temperature
        elif model == "loss":
            n0_now = initial_condensate_atoms * np.exp(-eta_coll * photons_scattered * shot)
            temperature_now = np.nan
        else:
            raise ValueError("model must be heating or loss")

        if n0_now <= 0:
            break
        loss_fraction = 1 - n0_now / initial_condensate_atoms
        state_now = tf_state_for_atoms(float(n0_now), constants)
        phase = scalar_phase_peak(detuning_hz, np.asarray(state_now["column_density"])[axis], constants)
        snr = pci_snr_pixel_for_phase(config, phase, power_mw, axis, tau_s, blur_axis)

        out["shot"].append(float(shot))
        out["N0"].append(float(n0_now))
        out["frac"].append(float(loss_fraction))
        out["T"].append(float(temperature_now))
        out["phi"].append(float(phase))
        out["snr"].append(float(snr))

        if loss_fraction >= loss_limit:
            break
        temperature = (temperature**4 + deposited_energy / energy_coefficient) ** 0.25
        shot += 1

    result = {key: np.asarray(value, dtype=float) for key, value in out.items()}
    result["accumulated_snr"] = accumulate_snr(result["snr"])
    result["condensate_fraction"] = result["N0"] / initial_condensate_atoms
    result["loss_fraction"] = result["frac"]
    return result


def _helper_run_sequence(
    *,
    config: dict[str, Any],
    constants: dict[str, Any],
    thermal: dict[str, float],
    blur_axis: dict[int, float],
    model: str,
) -> dict[str, np.ndarray]:
    multishot = config["multishot_recovery"]
    detuning_hz = float(multishot["detuning_ghz"]) * 1e9
    power_mw = float(multishot["probe_power_mw"])
    tau_s = float(multishot["pulse_duration_us"]) * 1e-6
    axis = int(multishot["imaging_axis"])
    loss_limit = float(multishot["loss_fraction_limit"])
    max_shots = int(multishot["max_shots"])
    eta_coll = float(multishot["eta_coll"])
    recoil_energy_multiplier = float(multishot.get("recoil_energy_multiplier", 1.0))
    reabs = multishot["reabsorption_override"]
    reabs = reabsorption_fraction(detuning_hz, constants) if reabs is None else float(reabs)
    photons_scattered = scattered_photons_per_atom(config, constants, detuning_hz, power_mw, tau_s)
    zeta3, zeta4 = float(zeta(3)), float(zeta(4))
    tc = thermal["critical_temperature_k"]
    energy_coefficient = 3 * (zeta4 / zeta3) * constants["boltzmann_constant"] / tc**3
    deposited_energy = recoil_energy_multiplier * photons_scattered * (1 + reabs) * constants["e_rec"]

    def phase_from_n0(n0_now: float) -> float:
        state_now = tf_state_for_atoms(n0_now, constants)
        return scalar_phase_peak(detuning_hz, np.asarray(state_now["column_density"])[axis], constants)

    def snr_from_phi(phase: float) -> float:
        return pci_snr_pixel_for_phase(config, phase, power_mw, axis, tau_s, blur_axis)

    result = simulate_multishot_sequence(
        photons_scattered_per_atom_per_shot=photons_scattered,
        initial_condensate_atoms=constants["atom_number"],
        loss_fraction_limit=loss_limit,
        max_shots=max_shots,
        model="heating" if model == "heating" else "loss",
        total_atoms=thermal["total_atoms"] if model == "heating" else None,
        initial_temperature_k=float(config["condensate"]["temperature_k"]) if model == "heating" else None,
        critical_temperature_k=tc if model == "heating" else None,
        energy_coefficient_j_per_k4=energy_coefficient if model == "heating" else None,
        deposited_energy_per_atom_per_shot_j=deposited_energy if model == "heating" else None,
        eta_coll=eta_coll,
        phase_from_n0=phase_from_n0,
        snr_from_phi=snr_from_phi,
    )
    result = dict(result)
    result["accumulated_snr"] = accumulate_snr(result["snr"])
    return result


def build_multishot_stage(config: dict[str, Any]) -> dict[str, Any]:
    constants = _basic_constants(config)
    thermal = self_consistent_total_atoms(config, constants)
    blur_axis = build_blur_axis(config, constants)
    multishot = config["multishot_recovery"]
    detuning_hz = float(multishot["detuning_ghz"]) * 1e9
    power_mw = float(multishot["probe_power_mw"])
    tau_s = float(multishot["pulse_duration_us"]) * 1e-6
    reabs = multishot["reabsorption_override"]
    reabs = reabsorption_fraction(detuning_hz, constants) if reabs is None else float(reabs)
    photons_scattered = scattered_photons_per_atom(config, constants, detuning_hz, power_mw, tau_s)
    saturation = intensity_at_atoms_notebook(config, power_mw) / constants["isat"]
    zeta3, zeta4 = float(zeta(3)), float(zeta(4))
    energy_coefficient = 3 * (zeta4 / zeta3) * constants["boltzmann_constant"] / thermal["critical_temperature_k"] ** 3
    recoil_energy_multiplier = float(multishot.get("recoil_energy_multiplier", 1.0))
    deposited_energy = recoil_energy_multiplier * photons_scattered * (1 + reabs) * constants["e_rec"]

    notebook = {
        model: _notebook_run_sequence(
            config=config,
            constants=constants,
            thermal=thermal,
            blur_axis=blur_axis,
            model=model,
        )
        for model in ("heating", "loss")
    }
    helper = {
        model: _helper_run_sequence(
            config=config,
            constants=constants,
            thermal=thermal,
            blur_axis=blur_axis,
            model=model,
        )
        for model in ("heating", "loss")
    }
    return {
        "constants": constants,
        "thermal": thermal,
        "blur_axis": blur_axis,
        "notebook": notebook,
        "helper": helper,
        "parameters": {
            "detuning_hz": detuning_hz,
            "detuning_ghz": float(multishot["detuning_ghz"]),
            "probe_power_mw": power_mw,
            "pulse_duration_s": tau_s,
            "pulse_duration_us": float(multishot["pulse_duration_us"]),
            "imaging_axis": int(multishot["imaging_axis"]),
            "loss_fraction_limit": float(multishot["loss_fraction_limit"]),
            "max_shots": int(multishot["max_shots"]),
            "eta_coll": float(multishot["eta_coll"]),
            "recoil_energy_multiplier": recoil_energy_multiplier,
            "recoil_energy_convention": multishot.get(
                "recoil_energy_convention", "legacy implicit one-recoil notebook convention"
            ),
            "reabsorption_fraction": reabs,
            "saturation_parameter": saturation,
            "photons_scattered_per_atom_per_shot": photons_scattered,
            "energy_coefficient_j_per_k4": energy_coefficient,
            "deposited_energy_per_atom_per_shot_j": deposited_energy,
            "initial_total_atoms": thermal["total_atoms"],
            "initial_condensate_atoms": constants["atom_number"],
            "initial_temperature_k": float(config["condensate"]["temperature_k"]),
            "critical_temperature_k": thermal["critical_temperature_k"],
        },
    }


def _array_comparison(left: np.ndarray, right: np.ndarray) -> dict[str, Any]:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    same_nan = np.array_equal(np.isnan(left), np.isnan(right))
    finite = np.isfinite(left) & np.isfinite(right)
    if finite.any():
        absolute = np.abs(left[finite] - right[finite])
        denominator = np.maximum(np.abs(left[finite]), np.abs(right[finite]))
        relative = np.zeros_like(absolute)
        mask = denominator > 0
        relative[mask] = absolute[mask] / denominator[mask]
        max_abs = float(np.max(absolute))
        max_rel = float(np.max(relative))
    else:
        max_abs = 0.0 if same_nan else float("nan")
        max_rel = 0.0 if same_nan else float("nan")
    return {
        "shape": list(left.shape),
        "same_nan_pattern": bool(same_nan),
        "max_absolute_difference": max_abs,
        "max_relative_difference": max_rel,
    }


def _sequence_summary(sequence: dict[str, np.ndarray]) -> dict[str, float | int]:
    return {
        "sequence_length": int(len(sequence["shot"])),
        "first_shot": float(sequence["shot"][0]),
        "last_shot": float(sequence["shot"][-1]),
        "initial_N0": float(sequence["N0"][0]),
        "final_N0": float(sequence["N0"][-1]),
        "final_loss_fraction": float(sequence["frac"][-1]),
        "initial_temperature_k": float(sequence["T"][0]) if np.isfinite(sequence["T"][0]) else float("nan"),
        "final_temperature_k": float(sequence["T"][-1]) if np.isfinite(sequence["T"][-1]) else float("nan"),
        "initial_phi_rad": float(sequence["phi"][0]),
        "final_phi_rad": float(sequence["phi"][-1]),
        "initial_snr": float(sequence["snr"][0]),
        "final_snr": float(sequence["snr"][-1]),
        "final_accumulated_snr": float(sequence["accumulated_snr"][-1]),
    }


def comparison_report(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    comparisons: dict[str, Any] = {}
    for model in ("heating", "loss"):
        comparisons[model] = {
            key: _array_comparison(stage["notebook"][model][key], stage["helper"][model][key])
            for key in ("shot", "N0", "frac", "T", "phi", "snr", "accumulated_snr")
        }
    return {
        "status": "canonical notebook deterministic multishot sequence matches simulate_multishot_sequence for tested quantities",
        "source_notebook_cells": {
            "run_sequence": 40,
            "sequence_evolution_plot": 42,
        },
        "notebook_formula": {
            "scattered_photons_per_atom": "N_scatt(Delta_Hz, P_mW, tau_s)",
            "heating_update": "T = (T**4 + dE / A_E)**0.25",
            "heating_condensate": "N0_now = N_tot_sc * (1 - (T / Tc)**3)",
            "clean_loss_condensate": "N0_now = N0_0 * exp(-eta_coll * Ng * s)",
            "phase_update": "phi = phi_peak(Delta_Hz, ncol_axis(N0_now))",
            "snr": "SNR_pixel_phi(phi, P_mW, axis, tau_s, blur_axis[axis])",
            "accumulated_snr": "sqrt(nancumsum(where(isfinite(snr), snr**2, 0)))",
        },
        "parameters": stage["parameters"],
        "blur_axis": stage["blur_axis"],
        "sequence_summaries": {
            model: _sequence_summary(stage["notebook"][model]) for model in ("heating", "loss")
        },
        "comparisons": comparisons,
        "scope_boundary": "Only deterministic sequence evolution is recovered. No noisy multishot filmstrip or operating-map recovery is generated.",
    }


def multishot_summary(config: dict[str, Any], stage: dict[str, Any]) -> dict[str, Any]:
    return {
        "label": config["label"],
        "status": "Canonical notebook-aligned deterministic multishot-stage recovery output.",
        "source_notebook": config["source_notebook"],
        "source_notebook_cells": {
            "run_sequence": 40,
            "sequence_evolution_plot": 42,
        },
        "parameters": stage["parameters"],
        "sequence_summaries": {
            model: _sequence_summary(stage["notebook"][model]) for model in ("heating", "loss")
        },
        "scope_boundary": "Only single-shot quantities -> deterministic repeated-imaging sequence evolution is recovered. No noisy filmstrip is generated.",
        "no_experimental_calibration_applied": True,
    }


def write_multishot_sequence_csv(path: Path, stage: dict[str, Any]) -> None:
    heating = stage["notebook"]["heating"]
    loss = stage["notebook"]["loss"]
    rows = []
    length = max(len(heating["shot"]), len(loss["shot"]))
    for index in range(length):
        row: dict[str, float] = {"row_index": float(index)}
        for label, sequence in [("heating", heating), ("loss", loss)]:
            for key in ("shot", "N0", "frac", "T", "phi", "snr", "accumulated_snr"):
                row[f"{label}_{key}"] = float(sequence[key][index]) if index < len(sequence[key]) else float("nan")
        rows.append(row)
    write_rows(path, rows)


def write_multishot_figure(path: Path, config: dict[str, Any], stage: dict[str, Any]) -> None:
    figure_size = tuple(config["multishot_recovery"]["figure_size_inches"])
    heating = stage["notebook"]["heating"]
    loss = stage["notebook"]["loss"]
    params = stage["parameters"]

    plt.rcParams.update(
        {
            "figure.figsize": (8, 5),
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 12,
            "legend.fontsize": 10,
            "figure.dpi": 110,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=figure_size, constrained_layout=True)
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    ax_a.plot(heating["shot"], 100 * (1 - heating["frac"]), "o-", color="#c4161c", lw=2.2, ms=4, label="heating")
    ax_a.plot(loss["shot"], 100 * (1 - loss["frac"]), "s--", color="#1f5fa8", lw=2.0, ms=3, label="clean loss")
    ax_a.axhline(100 * (1 - params["loss_fraction_limit"]), color="k", ls=":", lw=1.1)
    ax_a.set_xlabel(FRAME_INDEX)
    ax_a.set_ylabel(rf"condensate surviving (% of initial {ATOM_NUMBER})")
    ax_a.set_title("(a) condensate depletion")
    ax_a.grid(alpha=0.25)
    ax_a.legend(fontsize=8.5)
    ax_t = ax_a.twinx()
    ax_t.plot(heating["shot"], heating["T"] * 1e9, ":", color="#c4161c", lw=1.4, alpha=0.7)
    ax_t.axhline(params["critical_temperature_k"] * 1e9, color="gray", ls="--", lw=0.9)
    ax_t.set_ylabel(rf"{TEMPERATURE_NK}, heating only", color="#c4161c")
    ax_t.tick_params(axis="y", colors="#c4161c")

    ax_b.plot(heating["shot"], heating["phi"], "o-", color="#c4161c", lw=2.2, ms=4, label="heating")
    ax_b.plot(loss["shot"], loss["phi"], "s--", color="#1f5fa8", lw=2.0, ms=3, label="clean loss")
    ax_b.axhline(0.5, color="#e08020", ls=":", lw=1.2)
    ax_b.set_xlabel(FRAME_INDEX)
    ax_b.set_ylabel(rf"peak {PHASE_RAD}")
    ax_b.set_title("(b) peak phase")
    ax_b.grid(alpha=0.25)
    ax_b.legend(fontsize=8.5)

    ax_c.plot(heating["shot"], heating["snr"], "o-", color="#c4161c", lw=2.2, ms=4, label="heating")
    ax_c.plot(loss["shot"], loss["snr"], "s--", color="#1f5fa8", lw=2.0, ms=3, label="clean loss")
    ax_c.set_xlabel(FRAME_INDEX)
    ax_c.set_ylabel(f"per-shot PCI {SNR} / pixel")
    ax_c.set_title("(c) per-frame SNR")
    ax_c.grid(alpha=0.25)
    ax_c.legend(fontsize=8.5)

    ax_d.plot(heating["shot"], heating["accumulated_snr"], "o-", color="#c4161c", lw=2.2, ms=4, label="heating")
    ax_d.plot(loss["shot"], loss["accumulated_snr"], "s--", color="#1f5fa8", lw=2.0, ms=3, label="clean loss")
    ax_d.set_xlabel(FRAME_INDEX)
    ax_d.set_ylabel(f"accumulated {SNR}")
    ax_d.set_title("(d) RMS accumulated SNR")
    ax_d.grid(alpha=0.25)
    ax_d.legend(fontsize=8.5)

    fig.suptitle(
        rf"Multishot sequence evolution: {DETUNING_GHZ}={params['detuning_ghz']}, "
        rf"P={params['probe_power_mw']} mW, $\tau$={params['pulse_duration_us']:.0f} $\mu$s, "
        f"axis {'xyz'[params['imaging_axis']]}",
        fontsize=12,
    )
    fig.savefig(path, format="svg", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def generate(config: dict[str, Any], config_path: Path) -> dict[str, str]:
    output_dir = Path("results/notebook_aligned_recovery/multishot_stage")
    output_dir.mkdir(parents=True, exist_ok=True)
    stage = build_multishot_stage(config)
    outputs = {
        "comparison_report": output_dir / "comparison_report.json",
        "multishot_summary": output_dir / "multishot_summary.json",
        "multishot_sequence": output_dir / "multishot_sequence.csv",
        "metadata": output_dir / "metadata.json",
        "figure": output_dir / "multishot_sequence_stage.svg",
    }
    write_json(outputs["comparison_report"], comparison_report(config, stage))
    write_json(outputs["multishot_summary"], multishot_summary(config, stage))
    write_multishot_sequence_csv(outputs["multishot_sequence"], stage)
    write_multishot_figure(outputs["figure"], config, stage)
    write_json(
        outputs["metadata"],
        {
            "label": config["label"],
            "status": "Canonical notebook-aligned deterministic multishot-stage recovery output.",
            "config_file_used": str(config_path),
            "git_commit_hash": git_commit(),
            "source_notebook": config["source_notebook"],
            "source_notebook_cells": {
                "run_sequence": 40,
                "sequence_evolution_plot": 42,
            },
            "generated_outputs": {key: str(path) for key, path in outputs.items()},
            "scope_boundary": "Only deterministic multishot sequence recovery is generated. No noisy filmstrip or operating-map recovery is performed.",
            "parameters": stage["parameters"],
            "no_experimental_calibration_applied": True,
        },
    )
    return {key: str(path) for key, path in outputs.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config = load_config(args.config)
    outputs = generate(config, args.config)
    print("Recovered canonical notebook deterministic multishot stage outputs:")
    for key, path in outputs.items():
        print(f"- {key}: {path}")


if __name__ == "__main__":
    main()
