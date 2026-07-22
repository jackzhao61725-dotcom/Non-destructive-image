"""Physical parameters and fit-coordinate transforms for reconstruction."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
from numpy.typing import ArrayLike, NDArray


PARAMETER_NAMES = (
    "column_density_peak_m2",
    "y0_um",
    "z0_um",
    "radius_y_um",
    "radius_z_um",
)


@dataclass(frozen=True)
class SmoothTFParameters:
    """Parameters of a projected, free-radius Thomas-Fermi-like profile."""

    column_density_peak_m2: float
    y0_um: float
    z0_um: float
    radius_y_um: float
    radius_z_um: float

    def __post_init__(self) -> None:
        values = np.asarray(tuple(self.as_dict().values()), dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError("all smooth-TF parameters must be finite")
        if self.column_density_peak_m2 <= 0:
            raise ValueError("column_density_peak_m2 must be positive")
        if self.radius_y_um <= 0 or self.radius_z_um <= 0:
            raise ValueError("smooth-TF radii must be positive")

    def as_dict(self) -> dict[str, float]:
        return {
            "column_density_peak_m2": float(self.column_density_peak_m2),
            "y0_um": float(self.y0_um),
            "z0_um": float(self.z0_um),
            "radius_y_um": float(self.radius_y_um),
            "radius_z_um": float(self.radius_z_um),
        }


@dataclass(frozen=True)
class SmoothTFBounds:
    """Truth-independent physical bounds for a smooth-TF reconstruction."""

    lower: SmoothTFParameters
    upper: SmoothTFParameters

    def __post_init__(self) -> None:
        lower = to_internal(self.lower)
        upper = to_internal(self.upper)
        if np.any(lower >= upper):
            invalid = [
                name
                for name, lo, hi in zip(PARAMETER_NAMES, lower, upper, strict=True)
                if lo >= hi
            ]
            raise ValueError(f"lower bounds must be below upper bounds: {invalid}")

    @classmethod
    def from_mapping(cls, values: Mapping[str, ArrayLike]) -> "SmoothTFBounds":
        """Build bounds from absolute physical ranges.

        The mapping must contain each parameter name and two absolute values.
        Keys that encode a fraction of a simulated truth are rejected so that
        production fit construction cannot accidentally depend on ground truth.
        """

        forbidden = [key for key in values if "truth" in key.lower()]
        if forbidden:
            raise ValueError(
                "reconstruction bounds must not depend on simulated truth: "
                + ", ".join(forbidden)
            )
        missing = [name for name in PARAMETER_NAMES if name not in values]
        if missing:
            raise KeyError(f"missing smooth-TF bounds: {missing}")
        ranges: dict[str, tuple[float, float]] = {}
        for name in PARAMETER_NAMES:
            pair = np.asarray(values[name], dtype=float)
            if pair.shape != (2,):
                raise ValueError(f"{name} bounds must contain exactly two values")
            ranges[name] = (float(pair[0]), float(pair[1]))
        return cls(
            lower=SmoothTFParameters(*(ranges[name][0] for name in PARAMETER_NAMES)),
            upper=SmoothTFParameters(*(ranges[name][1] for name in PARAMETER_NAMES)),
        )


def to_internal(parameters: SmoothTFParameters) -> NDArray[np.floating]:
    """Transform positive scale parameters to logarithmic fit coordinates."""

    return np.asarray(
        [
            np.log(parameters.column_density_peak_m2),
            parameters.y0_um,
            parameters.z0_um,
            np.log(parameters.radius_y_um),
            np.log(parameters.radius_z_um),
        ],
        dtype=float,
    )


def from_internal(vector: ArrayLike) -> SmoothTFParameters:
    """Transform a five-element optimiser vector to physical parameters."""

    values = np.asarray(vector, dtype=float)
    if values.shape != (5,):
        raise ValueError("the smooth-TF optimiser vector must contain five values")
    return SmoothTFParameters(
        column_density_peak_m2=float(np.exp(values[0])),
        y0_um=float(values[1]),
        z0_um=float(values[2]),
        radius_y_um=float(np.exp(values[3])),
        radius_z_um=float(np.exp(values[4])),
    )


def covariance_to_physical(
    parameters: SmoothTFParameters,
    covariance_internal: ArrayLike,
) -> NDArray[np.floating]:
    """Propagate a local covariance from fit coordinates to physical units."""

    covariance = np.asarray(covariance_internal, dtype=float)
    if covariance.shape != (5, 5):
        raise ValueError("the smooth-TF covariance must have shape (5, 5)")
    transform = np.diag(
        [
            parameters.column_density_peak_m2,
            1.0,
            1.0,
            parameters.radius_y_um,
            parameters.radius_z_um,
        ]
    )
    return transform @ covariance @ transform.T
