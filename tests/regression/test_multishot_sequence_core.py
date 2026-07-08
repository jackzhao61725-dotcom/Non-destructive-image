from __future__ import annotations

import pytest

np = pytest.importorskip("numpy")

from non_destructive_image import accumulate_snr, simulate_multishot_sequence


def test_accumulate_snr_matches_notebook_rms_convention() -> None:
    snr = np.array([3.0, 4.0, np.nan, 12.0])

    np.testing.assert_allclose(accumulate_snr(snr), np.array([3.0, 5.0, 5.0, 13.0]))


def test_simulate_multishot_sequence_clean_loss_matches_hand_calculation() -> None:
    result = simulate_multishot_sequence(
        photons_scattered_per_atom_per_shot=0.1,
        initial_condensate_atoms=100.0,
        loss_fraction_limit=0.35,
        max_shots=10,
        model="loss",
        eta_coll=2.0,
        phase_from_n0=lambda n0: n0 / 100.0,
        snr_from_phi=lambda phi: 2.0 * phi,
    )

    expected_shots = np.arange(4, dtype=float)
    expected_n0 = 100.0 * np.exp(-2.0 * 0.1 * expected_shots)
    expected_loss = 1 - expected_n0 / 100.0

    np.testing.assert_allclose(result["shot"], expected_shots)
    np.testing.assert_allclose(result["N0"], expected_n0)
    np.testing.assert_allclose(result["loss_fraction"], expected_loss)
    np.testing.assert_allclose(result["frac"], expected_loss)
    np.testing.assert_allclose(result["condensate_fraction"], expected_n0 / 100.0)
    np.testing.assert_allclose(result["phi"], expected_n0 / 100.0)
    np.testing.assert_allclose(result["snr"], 2.0 * expected_n0 / 100.0)
    assert result["loss_fraction"][-1] >= 0.35


def test_simulate_multishot_sequence_heating_update_matches_hand_calculation() -> None:
    total_atoms = 1000.0
    initial_temperature = 1.0
    critical_temperature = 10.0
    initial_condensate_atoms = total_atoms * (1 - (initial_temperature / critical_temperature) ** 3)
    deposited_over_coefficient = 15.0

    result = simulate_multishot_sequence(
        photons_scattered_per_atom_per_shot=0.2,
        initial_condensate_atoms=initial_condensate_atoms,
        loss_fraction_limit=0.02,
        max_shots=10,
        model="heating",
        total_atoms=total_atoms,
        initial_temperature_k=initial_temperature,
        critical_temperature_k=critical_temperature,
        energy_coefficient_j_per_k4=2.0,
        deposited_energy_per_atom_per_shot_j=2.0 * deposited_over_coefficient,
        phase_from_n0=lambda n0: n0 ** (3 / 5),
        snr_from_phi=lambda phi: phi / 10.0,
    )

    expected_temperatures = []
    expected_n0 = []
    temperature = initial_temperature
    for _shot in result["shot"]:
        expected_temperatures.append(temperature)
        expected_n0.append(total_atoms * (1 - (temperature / critical_temperature) ** 3))
        temperature = (temperature**4 + deposited_over_coefficient) ** 0.25

    expected_temperatures = np.asarray(expected_temperatures)
    expected_n0 = np.asarray(expected_n0)

    np.testing.assert_allclose(result["T"], expected_temperatures)
    np.testing.assert_allclose(result["N0"], expected_n0)
    np.testing.assert_allclose(result["phi"], expected_n0 ** (3 / 5))
    np.testing.assert_allclose(result["snr"], expected_n0 ** (3 / 5) / 10.0)
    assert result["loss_fraction"][-1] >= 0.02


def test_simulate_multishot_sequence_outputs_have_consistent_lengths() -> None:
    result = simulate_multishot_sequence(
        photons_scattered_per_atom_per_shot=0.05,
        initial_condensate_atoms=50.0,
        loss_fraction_limit=0.2,
        max_shots=20,
        model="loss",
        eta_coll=1.5,
    )

    lengths = {len(value) for value in result.values()}
    assert lengths == {len(result["shot"])}


def test_simulate_multishot_sequence_is_deterministic_for_same_inputs() -> None:
    kwargs = {
        "photons_scattered_per_atom_per_shot": 0.07,
        "initial_condensate_atoms": 80.0,
        "loss_fraction_limit": 0.25,
        "max_shots": 20,
        "model": "loss",
        "eta_coll": 1.2,
        "phase_from_n0": lambda n0: 0.01 * n0,
        "snr_from_phi": lambda phi: phi**2,
    }

    first = simulate_multishot_sequence(**kwargs)
    second = simulate_multishot_sequence(**kwargs)

    for key in first:
        np.testing.assert_allclose(first[key], second[key])
