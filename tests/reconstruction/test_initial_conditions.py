from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from scipy import special

from non_destructive_image.reconstruction.initial_conditions import (
    InitialConditionTruth,
    build_initial_condition_truth,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
SUITE_CONFIG = (
    REPO_ROOT / "configs" / "dpfi_initial_condition_suite_v1_orca_fusion_m10.json"
)


@pytest.fixture(scope="module")
def suite_config() -> dict[str, object]:
    return json.loads(SUITE_CONFIG.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def active_object_grid() -> tuple[np.ndarray, np.ndarray]:
    ngrid = 306
    field_of_view_m = 100e-6
    axis_m = (np.arange(ngrid, dtype=float) - ngrid // 2) * (
        field_of_view_m / ngrid
    )
    return np.meshgrid(axis_m, axis_m)


def _conditions_by_id(config: dict[str, object]) -> dict[str, dict[str, object]]:
    return {
        condition["id"]: condition
        for condition in config["initial_conditions"]  # type: ignore[index,union-attr]
    }


def _build(
    condition: dict[str, object],
    config: dict[str, object],
    grid: tuple[np.ndarray, np.ndarray],
) -> InitialConditionTruth:
    return build_initial_condition_truth(
        grid[0],
        grid[1],
        condition,
        config["atomic_constants"],  # type: ignore[arg-type]
    )


def test_declared_suite_generators_produce_immutable_nonnegative_truths(
    suite_config: dict[str, object],
    active_object_grid: tuple[np.ndarray, np.ndarray],
) -> None:
    conditions = suite_config["initial_conditions"]
    truths = [
        _build(condition, suite_config, active_object_grid)
        for condition in conditions  # type: ignore[union-attr]
    ]

    assert len(truths) == 7
    assert {truth.generator_name for truth in truths} == {
        "contact_tf",
        "bimodal_ideal_gas",
        "tf_roton_modulation",
        "gaussian_chain",
    }
    for truth in truths:
        density = truth.column_density_m2
        assert truth.morphology.name == truth.condition_id
        assert density.shape == active_object_grid[0].shape
        assert np.all(np.isfinite(density))
        assert np.all(density >= 0.0)
        assert np.any(density > 0.0)
        assert not density.flags.writeable
        assert all(not component.flags.writeable for component in truth.components.values())
        np.testing.assert_array_equal(
            sum(truth.components.values(), start=np.zeros_like(density)),
            density,
        )


def test_reference_contact_tf_matches_active_repository_condensate(
    suite_config: dict[str, object],
    active_object_grid: tuple[np.ndarray, np.ndarray],
) -> None:
    condition = _conditions_by_id(suite_config)[
        "oxford_reference_pure_bec"
    ]
    truth = _build(condition, suite_config, active_object_grid)
    metadata = truth.derived_metadata

    assert condition["trap_frequencies_hz"] == [293.0, 14.0, 233.0]
    assert metadata["tf_peak_column_density_m2"] == pytest.approx(
        5.3759624525784675e14, rel=2e-8
    )
    radii_um = metadata["tf_radii_um"]
    assert radii_um[1] == pytest.approx(24.84, rel=2e-3)  # type: ignore[index]
    assert radii_um[2] == pytest.approx(1.49, rel=3e-3)  # type: ignore[index]
    assert metadata["tf_atom_number_check"] == pytest.approx(
        condition["condensate_atoms"], rel=1e-12
    )
    assert metadata["grid_captured_fraction_of_declared_total"] == pytest.approx(
        1.0, abs=1e-3
    )
    assert truth.source_metadata["relationship_to_existing_reference"].startswith(
        "This is the existing reference condensate"
    )


def test_bimodal_halo_uses_declared_population_and_saturated_bose_shape(
    suite_config: dict[str, object],
    active_object_grid: tuple[np.ndarray, np.ndarray],
) -> None:
    condition = _conditions_by_id(suite_config)["oxford_thermal_4650ms"]
    truth = _build(condition, suite_config, active_object_grid)
    components = truth.components
    metadata = truth.derived_metadata

    assert np.count_nonzero(components["condensate_core"]) == 0
    np.testing.assert_array_equal(components["thermal_halo"], truth.column_density_m2)
    assert metadata["declared_condensate_atoms"] == 0.0
    assert metadata["grid_captured_fraction_of_declared_total"] == pytest.approx(
        0.820115, rel=3e-5
    )

    constants = suite_config["atomic_constants"]
    mass = constants["mass_number"] * constants["atomic_mass_unit_kg"]  # type: ignore[index,operator]
    temperature_k = condition["temperature_nk"] * 1e-9  # type: ignore[operator]
    beta = 1.0 / (constants["boltzmann_constant_j_k"] * temperature_k)  # type: ignore[index,operator]
    _, omega_y, omega_z = 2.0 * np.pi * np.asarray(
        condition["trap_frequencies_hz"], dtype=float
    )
    expected_centre = (
        condition["thermal_atoms"]  # type: ignore[operator]
        * beta
        * mass
        * omega_y
        * omega_z
        / (2.0 * np.pi * special.zeta(3.0, 1.0))
        * (np.pi**2 / 6.0)
    )
    centre = active_object_grid[0].shape[0] // 2
    assert truth.column_density_m2[centre, centre] == pytest.approx(
        expected_centre, rel=1e-12
    )


def test_subresolution_modulation_stress_preserves_core_population_and_declares_period(
    suite_config: dict[str, object],
    active_object_grid: tuple[np.ndarray, np.ndarray],
) -> None:
    condition = _conditions_by_id(suite_config)[
        "er_subresolution_modulation_stress"
    ]
    truth = _build(condition, suite_config, active_object_grid)
    metadata = truth.derived_metadata

    assert metadata["modulation_axis"] == "y"
    assert metadata["modulation_period_um"] == pytest.approx(1.42607, rel=1e-5)
    assert metadata["grid_captured_modulated_condensate_core_atoms"] == pytest.approx(
        condition["condensate_atoms"], rel=2e-15
    )
    assert metadata["modulation_grid_population_correction"] == pytest.approx(
        1.0, abs=5e-3
    )
    assert metadata["grid_peak_column_density_m2"] > 3.0e15
    assert np.min(truth.components["modulated_condensate_core"]) >= 0.0


def test_resolved_three_peak_stress_is_population_normalised_with_explicit_positions(
    suite_config: dict[str, object],
    active_object_grid: tuple[np.ndarray, np.ndarray],
) -> None:
    condition_id = "er_resolved_three_peak_stress"
    expected_peak_m2 = 4.18655e15
    condition = _conditions_by_id(suite_config)[condition_id]
    truth = _build(condition, suite_config, active_object_grid)
    metadata = truth.derived_metadata

    assert metadata["grid_captured_condensate_chain_atoms"] == pytest.approx(
        condition["condensate_atoms"], rel=2e-12
    )
    assert metadata["grid_peak_column_density_m2"] == pytest.approx(
        expected_peak_m2, rel=2e-5
    )
    assert metadata["weighted_centre_y_um"] == pytest.approx(0.0, abs=1e-14)
    assert metadata["peak_displacements_um"] == tuple(
        condition["peak_displacements_um"]  # type: ignore[arg-type]
    )


def test_dispatch_rejects_unknown_generator(
    suite_config: dict[str, object],
    active_object_grid: tuple[np.ndarray, np.ndarray],
) -> None:
    condition = dict(
        _conditions_by_id(suite_config)["oxford_reference_pure_bec"]
    )
    condition["generator"] = "unimplemented_equilibrium_solver"
    with pytest.raises(ValueError, match="unsupported initial-condition generator"):
        _build(condition, suite_config, active_object_grid)


def test_generator_validation_rejects_silently_discarded_population(
    suite_config: dict[str, object],
    active_object_grid: tuple[np.ndarray, np.ndarray],
) -> None:
    condition = dict(
        _conditions_by_id(suite_config)["oxford_reference_pure_bec"]
    )
    condition["thermal_atoms"] = 1.0
    with pytest.raises(ValueError, match="does not represent a non-zero thermal"):
        _build(condition, suite_config, active_object_grid)
