from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest

from non_destructive_image.reconstruction.studies.initial_condition_suite import (
    aggregate_condition_rows,
    derive_trial_seed,
    validate_initial_condition_suite,
)
from non_destructive_image.reconstruction.studies.io import load_json
from scripts.run_dpfi_initial_condition_suite import (
    DEFAULT_CONFIG,
    validate_only,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
ACTIVE_CONFIG = (
    REPO_ROOT / "configs" / "dpfi_initial_condition_suite_v1_orca_fusion_m10.json"
)


def test_active_config_validates_without_fitting_and_has_the_frozen_dimensions() -> None:
    summary = validate_only(ACTIVE_CONFIG)

    assert summary["condition_count"] == 7
    assert tuple(summary["condition_ids"]) == (
        "oxford_reference_pure_bec",
        "oxford_bimodal_300ms",
        "oxford_thermal_4650ms",
        "er_ordinary_elongated_bec",
        "er_subresolution_modulation_stress",
        "er_same_trap_smooth_control",
        "er_resolved_three_peak_stress",
    )
    assert len(set(summary["condition_ids"])) == 7
    assert summary["fluence_mw_us"] == (90.0, 150.0, 300.0)
    assert summary["realizations_per_condition"] == 1
    assert summary["expected_trial_count"] == 21
    assert summary["parameter_count"] == 85
    assert summary["supported_cell_count"] == 7425


def test_trial_seed_mapping_is_stable_under_condition_reordering() -> None:
    config = load_json(ACTIVE_CONFIG)
    condition_ids = [condition["id"] for condition in config["initial_conditions"]]

    def seed_mapping(identifiers: list[str]) -> dict[tuple[float, str, int], int]:
        mapping: dict[tuple[float, str, int], int] = {}
        ensemble = config["ensemble"]
        for fluence in ensemble["fluence_mw_us"]:
            fluence_value = float(fluence)
            base_seed = int(
                ensemble["base_seed_by_fluence"][f"{fluence_value:g}"]
            )
            for condition_id in identifiers:
                for realization_index in range(
                    int(ensemble["realizations_per_condition"])
                ):
                    mapping[(fluence_value, condition_id, realization_index)] = (
                        derive_trial_seed(
                            base_seed,
                            condition_id,
                            realization_index,
                        )
                    )
        return mapping

    original = seed_mapping(condition_ids)
    reordered = seed_mapping(list(reversed(condition_ids)))

    assert original == reordered
    assert len(original) == 21
    assert len(set(original.values())) == 21


@pytest.mark.parametrize(
    ("key", "drifted_value"),
    (
        ("integrated_response_relative_error_tolerance", 0.11),
        ("centroid_position_error_um_tolerance", 0.651),
        ("major_rms_width_relative_error_tolerance", 0.11),
    ),
)
def test_suite_config_rejects_drift_from_each_frozen_usability_limit(
    key: str,
    drifted_value: float,
) -> None:
    config = deepcopy(load_json(ACTIVE_CONFIG))
    config["observable_usability"][key] = drifted_value

    with pytest.raises(ValueError, match=key):
        validate_initial_condition_suite(config, REPO_ROOT)


def _fake_trial_row(
    *,
    c_a: float | str,
    c_r: float | str,
    c_w: float | str,
    integral_passed: bool,
    centroid_passed: bool,
    width_passed: bool,
) -> dict[str, object]:
    return {
        "condition_id": "reference",
        "family": "test family",
        "feature_class": "smooth_condensate",
        "generator": "contact_tf",
        "observable_semantics": "support_restricted_column_density",
        "fit_success": True,
        "c_A": c_a,
        "c_r": c_r,
        "c_w": c_w,
        "integrated_response_passed": integral_passed,
        "centroid_position_passed": centroid_passed,
        "major_rms_width_passed": width_passed,
    }


def test_condition_aggregate_keeps_observable_decisions_separate_without_overall() -> None:
    rows = [
        _fake_trial_row(
            c_a=0.2,
            c_r=0.4,
            c_w=0.6,
            integral_passed=True,
            centroid_passed=True,
            width_passed=True,
        ),
        _fake_trial_row(
            c_a=0.8,
            c_r=1.2,
            c_w="",
            integral_passed=True,
            centroid_passed=False,
            width_passed=False,
        ),
    ]

    aggregate = aggregate_condition_rows(rows)

    assert len(aggregate) == 1
    result = aggregate[0]
    assert result["trial_count"] == 2
    assert result["maximum_c_A"] == pytest.approx(0.8)
    assert result["integrated_response_passed_all_trials"] is True
    assert result["maximum_c_r"] == pytest.approx(1.2)
    assert result["centroid_position_passed_all_trials"] is False
    assert result["supported_c_w_trials"] == 1
    assert result["major_rms_width_passed_all_trials"] is False
    assert not any("overall" in key.lower() for key in result)


def test_cli_default_is_the_active_initial_condition_config() -> None:
    assert DEFAULT_CONFIG == ACTIVE_CONFIG
    assert DEFAULT_CONFIG.is_file()
