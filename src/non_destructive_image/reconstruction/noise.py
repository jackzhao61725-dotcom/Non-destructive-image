"""Detector-noise helpers for reconstruction likelihoods and simulations."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def poisson_gaussian_variance(
    expected_electrons: ArrayLike,
    read_noise_electrons: float,
) -> NDArray[np.floating]:
    """Return the Gaussian-approximation variance in calibrated electrons."""

    expected = np.asarray(expected_electrons, dtype=float)
    if read_noise_electrons < 0:
        raise ValueError("read noise cannot be negative")
    return np.clip(expected, 0.0, None) + float(read_noise_electrons) ** 2


def simulate_poisson_gaussian_counts(
    expected_electrons: ArrayLike,
    *,
    read_noise_electrons: float,
    rng: np.random.Generator,
) -> NDArray[np.floating]:
    """Draw Poisson photoelectrons and independent Gaussian read noise."""

    expected = np.asarray(expected_electrons, dtype=float)
    if read_noise_electrons < 0:
        raise ValueError("read noise cannot be negative")
    return (
        rng.poisson(np.clip(expected, 0.0, None))
        + rng.normal(0.0, read_noise_electrons, expected.shape)
    ).astype(float)


def gaussian_quasi_deviance(
    observed_electrons: ArrayLike,
    expected_electrons: ArrayLike,
    *,
    read_noise_electrons: float,
) -> float:
    """Return a Poisson-Gaussian Gaussian-quasi-likelihood objective.

    This is retained as the Phase-1 approximation used to rank multi-start
    solutions. It must be benchmarked against a mixed Poisson-Gaussian
    likelihood before low-count dark-field results are promoted.
    """

    observed = np.asarray(observed_electrons, dtype=float)
    expected = np.asarray(expected_electrons, dtype=float)
    if observed.shape != expected.shape:
        raise ValueError("observed and expected count vectors must have the same shape")
    variance = poisson_gaussian_variance(expected, read_noise_electrons)
    return float(np.sum((observed - expected) ** 2 / variance + np.log(variance)))
