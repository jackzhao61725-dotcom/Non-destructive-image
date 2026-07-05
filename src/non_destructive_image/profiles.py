"""Thomas-Fermi profile helpers extracted from the reference notebook."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def thomas_fermi_profile_2d(
    coordinate_a: ArrayLike,
    coordinate_b: ArrayLike,
    radius_a: float,
    radius_b: float,
) -> NDArray[np.floating]:
    """Return the notebook's 2D Thomas-Fermi column profile.

    This is a direct extraction of the repeated notebook expression
    ``max(0, 1 - a**2/Ra**2 - b**2/Rb**2)**1.5``. The function does not change
    the exponent, normalisation, grid convention, or edge behaviour.
    """

    coordinate_a = np.asarray(coordinate_a)
    coordinate_b = np.asarray(coordinate_b)
    return np.maximum(
        0,
        1 - coordinate_a**2 / radius_a**2 - coordinate_b**2 / radius_b**2,
    ) ** 1.5
