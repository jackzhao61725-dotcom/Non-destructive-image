from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

import pytest

from scripts.audit_thesis_numerical_consistency import build_audit


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "configs" / "thesis_numerical_contract_v1.json"


@lru_cache(maxsize=1)
def _audit() -> dict:
    return build_audit(json.loads(CONTRACT.read_text(encoding="utf-8")))


def test_legacy_table_exposure_mismatch_is_detected() -> None:
    audit = _audit()
    legacy = next(row for row in audit["legacy_table_rows"] if row["Delta_over_2pi_GHz"] == 1.5)
    corrected = next(row for row in audit["corrected_legacy_rows"] if row["Delta_over_2pi_GHz"] == 1.5)

    assert legacy["SNR_exposure_us"] == 100.0
    assert legacy["N_max_exposure_us"] == 40.0
    assert legacy["PCI_SNR_total"] == pytest.approx(171.7, abs=0.1)
    assert corrected["SNR_exposure_us"] == corrected["N_max_exposure_us"] == 40.0
    assert corrected["PCI_SNR_total"] == pytest.approx(108.6, abs=0.1)


def test_canonical_fig_and_table_share_one_parameter_set() -> None:
    audit = _audit()
    rows = audit["accumulated_snr_rows"]
    assert all(row["power_mW"] == 3.5 for row in rows)
    assert all(row["exposure_us"] == 40.0 for row in rows)
    assert all(row["QE"] == 0.4 for row in rows)
    assert all(row["read_noise_e_rms"] == 7.0 for row in rows)
    assert max(row["PCI_SNR_total_shot_only"] for row in rows) - min(
        row["PCI_SNR_total_shot_only"] for row in rows
    ) < 0.05


def test_mislabelled_15us_budget_is_corrected() -> None:
    audit = _audit()
    budgets = {row["parameter_set"]: row for row in audit["budget_rows"]}
    stated_15 = budgets["legacy_label_claim_P2_tau15"]
    actual_40 = budgets["actual_legacy_calls_P2_tau40"]

    assert actual_40["N_max_clean_continuous"] == pytest.approx(52.0, abs=0.5)
    assert actual_40["N_stop_heating_continuous"] == pytest.approx(25.0, abs=0.5)
    assert actual_40["N_stop_heating_reabs_continuous"] == pytest.approx(24.0, abs=0.5)
    assert stated_15["N_max_clean_continuous"] == pytest.approx(
        actual_40["N_max_clean_continuous"] * 40 / 15,
        rel=1e-12,
    )
    assert stated_15["N_stop_heating_reabs_continuous"] == pytest.approx(
        actual_40["N_stop_heating_reabs_continuous"] * 40 / 15,
        rel=1e-12,
    )


def test_route_a_counting_and_reference_values_are_explicit() -> None:
    audit = _audit()
    route = next(row for row in audit["budget_rows"] if row["parameter_set"] == "canonical_route_a_P3p5_tau40")
    along = {row["Delta_over_2pi_GHz"]: row for row in audit["along_cigar_rows"]}

    assert route["N_max_clean_continuous"] == pytest.approx(29.581988680261396, rel=1e-12)
    assert route["N_stop_heating_reabs_continuous"] == pytest.approx(13.757594626407487, rel=1e-12)
    assert route["strict_integer_clean_frames"] == 29
    assert route["strict_integer_heating_reabs_frames"] == 13
    assert audit["constants"]["dimensionless_detuning_at_1p5GHz"] == pytest.approx(101.69491525423727)
    assert along[1.5]["peak_phase_rad"] == pytest.approx(4.247, abs=0.001)
    assert along[13.0]["peak_phase_rad"] == pytest.approx(0.490, abs=0.001)
    assert along[13.0]["N_max_clean_continuous"] == pytest.approx(3888, abs=1)


def test_full_model_reference_uses_same_route_a_but_distinct_observable() -> None:
    audit = _audit()
    rows = audit["full_model_rows"]
    values = {(row["mode"], row["noise_model"]): float(row["snr_total_full"]) for row in rows}

    assert all(int(float(row["n_frames_full"])) == 13 for row in rows)
    assert all(row["observable"] == "fixed-ROI diagonal-covariance matched-filter SNR" for row in rows)
    assert values[("PCI", "shot_plus_read_noise")] == pytest.approx(117.78942028852353, rel=1e-12)
    assert values[("DGI", "shot_plus_read_noise")] == pytest.approx(21.92223446100225, rel=1e-12)


def test_operating_point_snr_is_recomputed_with_explicit_40us() -> None:
    audit = _audit()
    rows = {(row["Delta_over_2pi_GHz"], row["axis"]): row for row in audit["operating_point_snr_rows"]}
    across = rows[(1.5, "x")]
    along = rows[(13.0, "y")]

    assert across["exposure_us"] == along["exposure_us"] == 40.0
    assert across["ideal_QE1_peak_pixel_SNR"] == pytest.approx(15.093850636466073, rel=1e-12)
    assert across["realistic_QE0p4_read7_peak_pixel_SNR"] == pytest.approx(5.781195504124932, rel=1e-12)
    assert across["realistic_QE0p4_read7_resolution_element_SNR"] == pytest.approx(11.691905893082325, rel=1e-12)
    assert along["ideal_QE1_peak_pixel_SNR"] == pytest.approx(36.45263879372922, rel=1e-12)
    assert along["realistic_QE0p4_read7_peak_pixel_SNR"] == pytest.approx(6.0149162768925635, rel=1e-12)
    assert along["realistic_QE0p4_read7_resolution_element_SNR"] == pytest.approx(8.525022820339919, rel=1e-12)
