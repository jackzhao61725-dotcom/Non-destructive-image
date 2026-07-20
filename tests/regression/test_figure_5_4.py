from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from scripts.generate_figure_5_4_snr_panel import DEFAULT_CONFIG, generate


@pytest.fixture(scope="module")
def generated(tmp_path_factory: pytest.TempPathFactory) -> tuple[dict[str, Path], dict, list[dict[str, str]]]:
    output = tmp_path_factory.mktemp("figure_5_4")
    outputs = generate(DEFAULT_CONFIG, output)
    metadata = json.loads(outputs["metadata"].read_text(encoding="utf-8"))
    with outputs["csv"].open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return outputs, metadata, rows


def test_figure_5_4_sequence_counts_match_screen(generated: tuple[dict[str, Path], dict, list[dict[str, str]]]) -> None:
    _, metadata, rows = generated
    assert metadata["detuning_ghz"] == pytest.approx(1.5)
    assert metadata["displayed_frames"] == 25
    assert metadata["sequence_counts"] == {
        "30": {
            "depletion_limited_frames": 32,
            "usable_frames": 2,
            "photoelectrons_per_incident_i0_pixel": pytest.approx(344.77282379399935),
        },
        "50": {
            "depletion_limited_frames": 19,
            "usable_frames": 19,
            "photoelectrons_per_incident_i0_pixel": pytest.approx(574.621372989999),
        },
        "90": {
            "depletion_limited_frames": 10,
            "usable_frames": 10,
            "photoelectrons_per_incident_i0_pixel": pytest.approx(1034.318471381998),
        },
        "150": {
            "depletion_limited_frames": 6,
            "usable_frames": 6,
            "photoelectrons_per_incident_i0_pixel": pytest.approx(1723.8641189699968),
        },
    }
    assert len(rows) == 4 * 25


def test_figure_5_4_status_bands_follow_the_two_limits(generated: tuple[dict[str, Path], dict, list[dict[str, str]]]) -> None:
    _, _, rows = generated
    statuses = {
        fluence: [
            row["status"]
            for row in rows
            if float(row["fluence_mw_us"]) == float(fluence)
        ]
        for fluence in (30, 50, 90, 150)
    }
    assert statuses[30] == ["usable"] * 2 + ["quality_limited"] * 23
    assert statuses[50][:19] == ["usable"] * 19
    assert statuses[90][:10] == ["usable"] * 10
    assert statuses[150][:6] == ["usable"] * 6
    for fluence in (50, 90, 150):
        assert "depletion_limited" in statuses[fluence]
        assert "both_limits" in statuses[fluence]


def test_figure_5_4_outputs_are_complete(generated: tuple[dict[str, Path], dict, list[dict[str, str]]]) -> None:
    outputs, _, _ = generated
    assert all(path.exists() and path.stat().st_size > 0 for path in outputs.values())
    svg = outputs["svg"].read_text(encoding="utf-8")
    assert "Image number" in svg
    assert "high-N band" in svg
    assert "loss above 30% only" in svg
    assert "both limits exceeded" in svg
