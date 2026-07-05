"""Atomic and Thomas-Fermi condensate helpers extracted from the notebook.

The reference notebook remains the source of parameter values. These helpers keep
all constants and experimental inputs explicit and only collect repeated algebra
into named functions.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class ThomasFermiState:
    """Derived Thomas-Fermi quantities for a condensate state."""

    trap_angular_frequencies: NDArray[np.floating]
    geometric_mean_frequency: float
    harmonic_oscillator_length: float
    chemical_potential: float
    chemical_potential_temperature: float
    peak_density: float
    radii: NDArray[np.floating]
    column_density: NDArray[np.floating]
    atom_number_check: float


def build_thomas_fermi_state(
    atom_number: float,
    scattering_length: float,
    trap_frequencies_hz: ArrayLike,
    atomic_mass: float,
    hbar: float,
    boltzmann_constant: float,
) -> ThomasFermiState:
    """Return the notebook's derived Thomas-Fermi condensate quantities.

    This is a direct extraction of the notebook algebra for ``omega``,
    ``omega_bar``, ``a_ho``, ``mu``, ``n_peak``, radii ``R``, column density
    ``n_col``, and the atom-number consistency check. No approximations or
    parameter values are changed.
    """

    trap_angular_frequencies = 2 * np.pi * np.asarray(trap_frequencies_hz)
    geometric_mean_frequency = float(trap_angular_frequencies.prod() ** (1 / 3))
    harmonic_oscillator_length = float(np.sqrt(hbar / (atomic_mass * geometric_mean_frequency)))
    chemical_potential = float(
        0.5
        * (15 * atom_number * scattering_length / harmonic_oscillator_length) ** (2 / 5)
        * hbar
        * geometric_mean_frequency
    )
    chemical_potential_temperature = float(chemical_potential / boltzmann_constant)
    peak_density = float(chemical_potential * atomic_mass / (4 * np.pi * hbar**2 * scattering_length))
    radii = np.sqrt(2 * chemical_potential / (atomic_mass * trap_angular_frequencies**2))
    column_density = (4 / 3) * peak_density * radii
    atom_number_check = float((8 * np.pi / 15) * peak_density * radii.prod())

    return ThomasFermiState(
        trap_angular_frequencies=trap_angular_frequencies,
        geometric_mean_frequency=geometric_mean_frequency,
        harmonic_oscillator_length=harmonic_oscillator_length,
        chemical_potential=chemical_potential,
        chemical_potential_temperature=chemical_potential_temperature,
        peak_density=peak_density,
        radii=radii,
        column_density=column_density,
        atom_number_check=atom_number_check,
    )


def recoil_quantities(
    hbar: float,
    wavevector: float,
    atomic_mass: float,
    boltzmann_constant: float,
) -> tuple[float, float, float]:
    """Return recoil energy, recoil temperature, and recoil velocity.

    The formulas match the notebook definitions ``E_rec = (hbar*k)**2/(2*m)``,
    ``T_rec = E_rec/kB``, and ``v_rec = hbar*k/m``.
    """

    recoil_energy = (hbar * wavevector) ** 2 / (2 * atomic_mass)
    recoil_temperature = recoil_energy / boltzmann_constant
    recoil_velocity = hbar * wavevector / atomic_mass
    return float(recoil_energy), float(recoil_temperature), float(recoil_velocity)
