"""Independent synthetic density morphologies for reconstruction stress tests.

These maps are not fitted object families and are never used as priors.  They
exercise resolved features that a reconstruction may encounter without
claiming to be equilibrium solutions of a dipolar Gross--Pitaevskii equation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class SyntheticMorphology:
    """Named non-negative truth map used only in synthetic assessment."""

    name: str
    column_density_m2: NDArray[np.floating]
    description: str
    feature_class: str

    def __post_init__(self) -> None:
        density = np.asarray(self.column_density_m2, dtype=float)
        if density.ndim != 2 or np.any(~np.isfinite(density)):
            raise ValueError("synthetic density must be a finite two-dimensional map")
        if np.any(density < 0) or not np.any(density > 0):
            raise ValueError("synthetic density must be non-negative and non-empty")
        object.__setattr__(self, "column_density_m2", density)


@dataclass(frozen=True)
class SyntheticMorphologySplit:
    """Disjoint analytic instances used for selection and frozen assessment."""

    calibration: tuple[SyntheticMorphology, ...]
    held_out: tuple[SyntheticMorphology, ...]

    def __post_init__(self) -> None:
        calibration_names = {case.name for case in self.calibration}
        held_out_names = {case.name for case in self.held_out}
        if not self.calibration or not self.held_out:
            raise ValueError("both morphology splits must be non-empty")
        overlap = calibration_names.intersection(held_out_names)
        if overlap:
            raise ValueError(
                "calibration and held-out morphology names overlap: "
                + ", ".join(sorted(overlap))
            )


def _validated_grids(
    y_grid_m: ArrayLike,
    z_grid_m: ArrayLike,
) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
    y = np.asarray(y_grid_m, dtype=float)
    z = np.asarray(z_grid_m, dtype=float)
    if y.ndim != 2 or y.shape != z.shape:
        raise ValueError("synthetic morphology grids must be same-shape 2D arrays")
    if np.any(~np.isfinite(y)) or np.any(~np.isfinite(z)):
        raise ValueError("synthetic morphology coordinates must be finite")
    return y * 1e6, z * 1e6


def _normalise_peak(
    density: NDArray[np.floating],
    peak_column_density_m2: float,
) -> NDArray[np.floating]:
    if not np.isfinite(peak_column_density_m2) or peak_column_density_m2 <= 0:
        raise ValueError("peak column density must be finite and positive")
    maximum = float(np.max(density))
    if maximum <= 0:
        raise ValueError("synthetic morphology has zero support")
    return peak_column_density_m2 * density / maximum


def _compact_envelope(
    y_um: NDArray[np.floating],
    z_um: NDArray[np.floating],
    *,
    centre_y_um: float,
    centre_z_um: float,
    radius_y_um: float,
    radius_z_um: float,
) -> NDArray[np.floating]:
    if radius_y_um <= 0 or radius_z_um <= 0:
        raise ValueError("synthetic morphology radii must be positive")
    q = (
        ((y_um - centre_y_um) / radius_y_um) ** 2
        + ((z_um - centre_z_um) / radius_z_um) ** 2
    )
    return np.clip(1.0 - q, 0.0, None) ** 1.5


def build_reference_morphology_suite(
    y_grid_m: ArrayLike,
    z_grid_m: ArrayLike,
    *,
    peak_column_density_m2: float,
    centre_y_um: float = 0.0,
    centre_z_um: float = 0.0,
    radius_y_um: float = 20.0,
    radius_z_um: float = 6.0,
) -> tuple[SyntheticMorphology, ...]:
    """Return analytic, model-mismatched maps spanning five feature classes.

    The suite covers a smooth reference, asymmetric deformation, weak resolved
    modulation, a compact fragmented chain and a local density notch.  The
    shapes are controlled stress tests, not predictions for a specific erbium
    state.
    """

    y_um, z_um = _validated_grids(y_grid_m, z_grid_m)
    envelope = _compact_envelope(
        y_um,
        z_um,
        centre_y_um=centre_y_um,
        centre_z_um=centre_z_um,
        radius_y_um=radius_y_um,
        radius_z_um=radius_z_um,
    )

    smooth = _normalise_peak(envelope, peak_column_density_m2)

    y_scaled = (y_um - centre_y_um) / radius_y_um
    z_scaled = (z_um - centre_z_um) / radius_z_um
    asymmetric = envelope * np.clip(1.0 + 0.45 * y_scaled - 0.20 * y_scaled * z_scaled, 0.0, None)
    asymmetric = _normalise_peak(asymmetric, peak_column_density_m2)

    modulation = 1.0 + 0.22 * np.cos(
        2.0 * np.pi * (y_um - centre_y_um) / (0.55 * radius_y_um) + 0.35
    )
    modulated = _normalise_peak(envelope * modulation, peak_column_density_m2)

    fragmented = np.zeros_like(envelope)
    for displacement, weight, width_scale in (
        (-0.48, 0.72, 1.00),
        (0.00, 1.00, 0.88),
        (0.46, 0.82, 1.06),
    ):
        fragmented += weight * np.exp(
            -0.5
            * (
                ((y_um - centre_y_um - displacement * radius_y_um) / (0.14 * radius_y_um * width_scale)) ** 2
                + ((z_um - centre_z_um) / (0.28 * radius_z_um)) ** 2
            )
        )
    fragmented *= envelope > 0
    fragmented = _normalise_peak(fragmented, peak_column_density_m2)

    notch = 1.0 - 0.65 * np.exp(
        -0.5
        * (
            ((y_um - centre_y_um - 0.12 * radius_y_um) / (0.10 * radius_y_um)) ** 2
            + ((z_um - centre_z_um) / (0.40 * radius_z_um)) ** 2
        )
    )
    notched = _normalise_peak(envelope * notch, peak_column_density_m2)

    return (
        SyntheticMorphology(
            "smooth_reference",
            smooth,
            "Compact smooth single cloud used as a low-complexity control.",
            "smooth",
        ),
        SyntheticMorphology(
            "asymmetric_single_cloud",
            asymmetric,
            "Smooth compact cloud with a skewed density distribution.",
            "asymmetry",
        ),
        SyntheticMorphology(
            "weak_longitudinal_modulation",
            modulated,
            "Compact cloud carrying a weak resolved longitudinal modulation.",
            "modulation",
        ),
        SyntheticMorphology(
            "fragmented_three_peak_chain",
            fragmented,
            "Unequal compact peaks within the reference cloud support.",
            "fragmentation",
        ),
        SyntheticMorphology(
            "local_density_notch",
            notched,
            "Smooth compact cloud with a local density depletion.",
            "local_defect",
        ),
    )


def build_calibration_held_out_morphology_suite(
    y_grid_m: ArrayLike,
    z_grid_m: ArrayLike,
    *,
    peak_column_density_m2: float,
    centre_y_um: float = 0.0,
    centre_z_um: float = 0.0,
    radius_y_um: float = 20.0,
    radius_z_um: float = 6.0,
) -> SyntheticMorphologySplit:
    """Build disjoint instances spanning the same five morphology classes.

    Candidate selection must see the *types* of structure the reconstruction
    is expected to represent.  Held-out assessment must nevertheless use maps
    that did not participate in that selection.  This function therefore
    supplies different analytic instances of smooth displacement, asymmetry,
    longitudinal modulation, compact fragmentation and a local density defect
    on the two sides of the split.  None is generated from the fitting basis or
    asserted to be an equilibrium dipolar state.
    """

    y_um, z_um = _validated_grids(y_grid_m, z_grid_m)

    def envelope(y_shift: float = 0.0, z_shift: float = 0.0) -> NDArray[np.floating]:
        return _compact_envelope(
            y_um,
            z_um,
            centre_y_um=centre_y_um + y_shift,
            centre_z_um=centre_z_um + z_shift,
            radius_y_um=radius_y_um,
            radius_z_um=radius_z_um,
        )

    def normalised(values: NDArray[np.floating]) -> NDArray[np.floating]:
        return _normalise_peak(values, peak_column_density_m2)

    base = envelope()
    y_scaled = (y_um - centre_y_um) / radius_y_um
    z_scaled = (z_um - centre_z_um) / radius_z_um

    calibration_smooth = normalised(base)
    calibration_asymmetric = normalised(
        base * np.clip(1.0 + 0.40 * y_scaled - 0.18 * y_scaled * z_scaled, 0.0, None)
    )
    calibration_modulated = normalised(
        base
        * (
            1.0
            + 0.20
            * np.cos(2.0 * np.pi * (y_um - centre_y_um) / (0.58 * radius_y_um) + 0.25)
        )
    )
    calibration_fragmented = np.zeros_like(base)
    for displacement, weight, width_scale in (
        (-0.30, 0.82, 1.00),
        (0.32, 1.00, 0.92),
    ):
        calibration_fragmented += weight * np.exp(
            -0.5
            * (
                (
                    (y_um - centre_y_um - displacement * radius_y_um)
                    / (0.17 * radius_y_um * width_scale)
                )
                ** 2
                + ((z_um - centre_z_um) / (0.30 * radius_z_um)) ** 2
            )
        )
    calibration_fragmented = normalised(calibration_fragmented * (base > 0))
    calibration_notched = normalised(
        base
        * (
            1.0
            - 0.55
            * np.exp(
                -0.5
                * (
                    ((y_um - centre_y_um) / (0.12 * radius_y_um)) ** 2
                    + ((z_um - centre_z_um) / (0.48 * radius_z_um)) ** 2
                )
            )
        )
    )

    shifted_envelope = envelope(0.08 * radius_y_um, -0.05 * radius_z_um)
    held_out_smooth = normalised(shifted_envelope)
    held_out_asymmetric = normalised(
        shifted_envelope
        * np.clip(1.0 - 0.34 * y_scaled + 0.22 * y_scaled * z_scaled, 0.0, None)
    )
    held_out_modulated = normalised(
        base
        * (
            1.0
            + 0.24
            * np.cos(2.0 * np.pi * (y_um - centre_y_um) / (0.46 * radius_y_um) + 1.05)
        )
    )
    held_out_fragmented = np.zeros_like(base)
    for displacement, weight, width_scale in (
        (-0.50, 0.70, 1.04),
        (-0.02, 1.00, 0.86),
        (0.45, 0.80, 1.10),
    ):
        held_out_fragmented += weight * np.exp(
            -0.5
            * (
                (
                    (y_um - centre_y_um - displacement * radius_y_um)
                    / (0.14 * radius_y_um * width_scale)
                )
                ** 2
                + ((z_um - centre_z_um) / (0.28 * radius_z_um)) ** 2
            )
        )
    held_out_fragmented = normalised(held_out_fragmented * (base > 0))
    held_out_notched = normalised(
        base
        * (
            1.0
            - 0.68
            * np.exp(
                -0.5
                * (
                    (
                        (y_um - centre_y_um - 0.18 * radius_y_um)
                        / (0.10 * radius_y_um)
                    )
                    ** 2
                    + ((z_um - centre_z_um + 0.08 * radius_z_um) / (0.38 * radius_z_um))
                    ** 2
                )
            )
        )
    )

    def case(
        name: str,
        density: NDArray[np.floating],
        description: str,
        feature_class: str,
    ) -> SyntheticMorphology:
        return SyntheticMorphology(name, density, description, feature_class)

    return SyntheticMorphologySplit(
        calibration=(
            case(
                "calibration_smooth_reference",
                calibration_smooth,
                "Centred smooth compact cloud used as the low-complexity control.",
                "smooth",
            ),
            case(
                "calibration_asymmetric_cloud",
                calibration_asymmetric,
                "Smooth cloud with a declared skew and weak transverse coupling.",
                "asymmetry",
            ),
            case(
                "calibration_weak_modulation",
                calibration_modulated,
                "Weak resolved longitudinal density modulation.",
                "modulation",
            ),
            case(
                "calibration_two_peak_fragment",
                calibration_fragmented,
                "Unequal compact two-peak chain inside the common support.",
                "fragmentation",
            ),
            case(
                "calibration_central_notch",
                calibration_notched,
                "Smooth cloud with a central resolution-scale density depletion.",
                "local_defect",
            ),
        ),
        held_out=(
            case(
                "held_out_shifted_smooth_cloud",
                held_out_smooth,
                "Translated smooth cloud not used for candidate selection.",
                "smooth",
            ),
            case(
                "held_out_opposite_asymmetry",
                held_out_asymmetric,
                "Translated cloud with a different asymmetry sign and strength.",
                "asymmetry",
            ),
            case(
                "held_out_modulation_phase",
                held_out_modulated,
                "Longitudinal modulation with unseen wavelength and phase.",
                "modulation",
            ),
            case(
                "held_out_three_peak_fragment",
                held_out_fragmented,
                "Unequal three-peak chain with unseen positions and widths.",
                "fragmentation",
            ),
            case(
                "held_out_off_centre_notch",
                held_out_notched,
                "Off-centre local depletion with unseen position and depth.",
                "local_defect",
            ),
        ),
    )
