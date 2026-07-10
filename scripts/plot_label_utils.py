"""Shared Matplotlib-safe labels for notebook-aligned recovery figures."""

from __future__ import annotations

LENGTH_UM = r"$\mu\mathrm{m}$"
DENSITY_CM3 = r"$\mathrm{cm}^{-3}$"
COLUMN_DENSITY_CM2 = r"$\mathrm{cm}^{-2}$"
DENSITY_M3 = r"$\mathrm{m}^{-3}$"
COLUMN_DENSITY_M2 = r"$\mathrm{m}^{-2}$"
PHASE_RAD = r"$\phi$ (rad)"
FARADAY_ANGLE_RAD = r"$\theta_F$ (rad)"
NORMALISED_INTENSITY = r"$I/I_0$"
DARK_FIELD_INTENSITY = r"$I_{\mathrm{dark}}/I_0$"
DUAL_PORT_SIGNAL = r"$S=(I_v-I_u)/(I_v+I_u)$"
ATOM_NUMBER = r"$N_0$"
TEMPERATURE_NK = r"$T$ (nK)"
SNR = "SNR"
LOSS_FRACTION = "Loss fraction"
FRAME_INDEX = "frame index"
POWER_MW = "Power (mW)"
EXPOSURE_US = r"Exposure ($\mu$s)"
READ_NOISE = r"Read noise ($e^-$)"
TRAP_FREQUENCY = r"$\bar{\omega}/2\pi$ (Hz)"
CHEMICAL_POTENTIAL = r"$\mu/k_B$ (nK)"
HARMONIC_OSCILLATOR_LENGTH = r"$a_{\mathrm{ho}}$ ($\mu\mathrm{m}$)"
TF_RADII = r"$R_i$ ($\mu\mathrm{m}$)"
PEAK_DENSITY = r"$n_0$ ($\mathrm{m}^{-3}$)"
PEAK_COLUMN_DENSITY = r"$\tilde{n}_i$ ($\mathrm{m}^{-2}$)"
DETUNING_HZ = r"$\Delta$ (Hz)"
DETUNING_GHZ = r"$\Delta$ (GHz)"


def coordinate_label(axis: str | None = None) -> str:
    """Return a coordinate-axis label in micrometres."""

    if axis:
        return rf"${axis}$ ({LENGTH_UM})"
    return rf"position ({LENGTH_UM})"


def cut_label(axis: str) -> str:
    """Return a legend label for a central lineout."""

    return rf"cut along ${axis}$"


def radius_legend_label(axis_label: str, radius_um: float) -> str:
    """Return a density-cut legend label with a radius in micrometres."""

    return rf"{axis_label}, $R={radius_um:.2f}\,\mu\mathrm{{m}}$"


def peak_column_density_symbol(integrated_axis: str) -> str:
    """Return the peak column-density scalar symbol for one integration axis."""

    return rf"$\tilde{{n}}_{integrated_axis}$"


def column_density_distribution_label(
    plane_axis_a: str,
    plane_axis_b: str,
    display_unit: str = "cm^-2",
) -> str:
    """Return a full 2D column-density distribution label.

    Use this for maps. Reserve ``peak_column_density_symbol(...)`` for scalar
    peak values such as those in the thesis parameter table.
    """

    unit = COLUMN_DENSITY_CM2 if display_unit == "cm^-2" else COLUMN_DENSITY_M2
    return rf"$n_{{\mathrm{{col}}}}({plane_axis_a},{plane_axis_b})$ ({unit})"
