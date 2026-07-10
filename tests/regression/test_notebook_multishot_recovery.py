from __future__ import annotations

import json
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

from scripts.recover_notebook_multishot_stage import build_multishot_stage, comparison_report


CONFIG_PATH = Path("configs/notebook_v1_defaults.json")


def _load_config() -> dict:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def test_notebook_multishot_recovery_matches_helper_sequences() -> None:
    config = _load_config()
    stage = build_multishot_stage(config)
    report = comparison_report(config, stage)

    for model in ["heating", "loss"]:
        for key in ["shot", "N0", "frac", "T", "phi", "snr", "accumulated_snr"]:
            comparison = report["comparisons"][model][key]
            assert comparison["same_nan_pattern"] is True
            assert comparison["max_absolute_difference"] == 0.0
            assert comparison["max_relative_difference"] == 0.0


def test_notebook_multishot_recovery_stable_parameters() -> None:
    config = _load_config()
    stage = build_multishot_stage(config)
    params = stage["parameters"]

    assert params["detuning_ghz"] == 1.5
    assert params["probe_power_mw"] == 3.5
    assert params["pulse_duration_us"] == 40.0
    assert params["imaging_axis"] == 0
    assert params["loss_fraction_limit"] == 0.3
    assert params["max_shots"] == 400
    assert params["eta_coll"] == 1.3
    assert params["photons_scattered_per_atom_per_shot"] == pytest.approx(
        0.009274742967987243,
        rel=1e-12,
    )
    assert params["reabsorption_fraction"] == pytest.approx(0.029708652968257532, rel=1e-12)
    assert params["critical_temperature_k"] == pytest.approx(2.168262079928517e-07, rel=1e-12)


def test_notebook_multishot_recovery_stable_sequence_outputs() -> None:
    config = _load_config()
    stage = build_multishot_stage(config)
    heating = stage["notebook"]["heating"]
    loss = stage["notebook"]["loss"]

    assert len(heating["shot"]) == 15
    assert len(loss["shot"]) == 31
    assert heating["shot"][-1] == 14.0
    assert loss["shot"][-1] == 30.0

    assert heating["N0"][0] == pytest.approx(25000.0, rel=1e-12)
    assert heating["N0"][-1] == pytest.approx(17369.583984817014, rel=1e-12)
    assert heating["frac"][-1] == pytest.approx(0.30521664060731946, rel=1e-12)
    assert heating["T"][-1] == pytest.approx(2.0543100787139788e-07, rel=1e-12)
    assert heating["phi"][0] == pytest.approx(0.20294165287929006, rel=1e-12)
    assert heating["phi"][-1] == pytest.approx(0.16311009404814675, rel=1e-12)
    assert heating["snr"][0] == pytest.approx(7.755454842736431, rel=1e-12)
    assert heating["accumulated_snr"][-1] == pytest.approx(27.281738255565195, rel=1e-12)

    assert loss["N0"][0] == pytest.approx(25000.0, rel=1e-12)
    assert loss["N0"][-1] == pytest.approx(17412.021337141854, rel=1e-12)
    assert loss["frac"][-1] == pytest.approx(0.3035191465143259, rel=1e-12)
    assert np.isnan(loss["T"]).all()
    assert loss["phi"][0] == pytest.approx(0.20294165287929014, rel=1e-12)
    assert loss["phi"][-1] == pytest.approx(0.16334908360737055, rel=1e-12)
    assert loss["snr"][0] == pytest.approx(7.755454842736431, rel=1e-12)
    assert loss["accumulated_snr"][-1] == pytest.approx(39.020132203959946, rel=1e-12)


def test_notebook_multishot_recovery_accumulated_snr_convention() -> None:
    config = _load_config()
    stage = build_multishot_stage(config)
    for model in ["heating", "loss"]:
        sequence = stage["notebook"][model]
        expected = np.sqrt(np.nancumsum(np.where(np.isfinite(sequence["snr"]), sequence["snr"] ** 2, 0.0)))
        np.testing.assert_allclose(sequence["accumulated_snr"], expected)
