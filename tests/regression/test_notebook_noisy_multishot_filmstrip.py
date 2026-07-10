from __future__ import annotations

import json
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

from scripts.recover_notebook_noisy_multishot_filmstrip import (
    build_noisy_multishot_filmstrip,
    comparison_report,
)


CONFIG_PATH = Path("configs/notebook_v1_defaults.json")


def _load_config() -> dict:
    return json.loads(CONFIG_PATH.read_text(encoding="utf-8"))


def test_noisy_multishot_filmstrip_matches_seeded_helper_replay() -> None:
    config = _load_config()
    stage = build_noisy_multishot_filmstrip(config)
    report = comparison_report(config, stage)

    assert report["primary_recovered_cell"] == 93
    assert report["selected_frame_indices"] == [0, 5, 10, 14]
    assert report["related_cell44_frame_indices"] == [0, 7, 14]
    assert report["rng_policy"]["rng_seed"] == 7
    assert report["rng_policy"]["exact_full_notebook_rng_reproduction"] is False

    for frame_report in report["frame_comparisons"]:
        for key in ["binned_image", "noisy_counts", "noisy_frame"]:
            assert frame_report[key]["max_absolute_difference"] == 0.0
            assert frame_report[key]["max_relative_difference"] == 0.0


def test_noisy_multishot_filmstrip_stable_selected_sequence_values() -> None:
    config = _load_config()
    stage = build_noisy_multishot_filmstrip(config)
    frames = stage["frames"]

    assert [frame["frame_index"] for frame in frames] == [0, 5, 10, 14]
    assert stage["photons_per_pixel"] == pytest.approx(1072.6265629146646, rel=1e-12)
    assert stage["read_noise_electrons"] == 7.0

    expected_n0 = [
        24999.999999999978,
        22251.24175853835,
        19529.04745904292,
        17369.583984817014,
    ]
    expected_loss = [
        8.881784197001252e-16,
        0.10995032965846607,
        0.2188381016382832,
        0.30521664060731946,
    ]
    expected_phi = [
        0.20294165287929006,
        0.18924296307881183,
        0.1749910449421481,
        0.16311009404814675,
    ]
    expected_accumulated = [
        7.755454842736431,
        18.383674274183647,
        24.047846590361974,
        27.281738255565195,
    ]
    for index, frame in enumerate(frames):
        assert frame["N0"] == pytest.approx(expected_n0[index], rel=1e-12)
        assert frame["loss_fraction"] == pytest.approx(expected_loss[index], rel=1e-12)
        assert frame["phi_rad"] == pytest.approx(expected_phi[index], rel=1e-12)
        assert frame["accumulated_snr"] == pytest.approx(expected_accumulated[index], rel=1e-12)


def test_noisy_multishot_filmstrip_frame_statistics_are_stable() -> None:
    config = _load_config()
    stage = build_noisy_multishot_filmstrip(config)
    frames = stage["frames"]

    expected_means = [
        0.9048341400329197,
        0.9045503168872816,
        0.9044529944801176,
        0.9043769748075294,
    ]
    expected_residual_std = [
        0.030427061235380748,
        0.030517108798712705,
        0.02956380887046506,
        0.029429505453540096,
    ]

    for index, frame in enumerate(frames):
        noisy = frame["noisy_frame"]
        residual = noisy - frame["binned_image"]
        assert noisy.shape == (68, 68)
        assert np.isfinite(noisy).all()
        assert np.mean(noisy) == pytest.approx(expected_means[index], rel=1e-12)
        assert np.std(residual) == pytest.approx(expected_residual_std[index], rel=1e-12)


def test_noisy_multishot_filmstrip_is_reproducible_for_same_config() -> None:
    config = _load_config()
    first = build_noisy_multishot_filmstrip(config)
    second = build_noisy_multishot_filmstrip(config)
    for left, right in zip(first["frames"], second["frames"]):
        np.testing.assert_allclose(left["noisy_frame"], right["noisy_frame"])
        np.testing.assert_allclose(left["noisy_counts"], right["noisy_counts"])
