from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from scripts.generate_figure_5_4_snr_panel import DEFAULT_CONFIG, generate


ROOT = Path(__file__).resolve().parents[2]
CANONICAL_OUTPUT = ROOT / "results" / "dissertation_plots_v2_orca_fusion" / "figure_5_4"


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
    assert metadata["displayed_frames"] == 15
    assert metadata["sequence_counts"] == {
        "90": {
            "depletion_limited_frames": 10,
            "usable_frames": 3,
            "photoelectrons_per_incident_i0_pixel": pytest.approx(220.58087277528466),
        },
        "150": {
            "depletion_limited_frames": 6,
            "usable_frames": 6,
            "photoelectrons_per_incident_i0_pixel": pytest.approx(367.63478795880775),
        },
        "300": {
            "depletion_limited_frames": 3,
            "usable_frames": 3,
            "photoelectrons_per_incident_i0_pixel": pytest.approx(735.2695759176155),
        },
    }
    assert metadata["atomic_response"]["kappa_F"] == pytest.approx(-45 / 91)
    assert metadata["optics"]["effective_numerical_aperture"] == pytest.approx(0.13)
    assert metadata["camera"]["read_noise_electrons_per_pixel_per_port"] == pytest.approx(0.7)
    assert metadata["result_status"] == {
        "status": "frozen_version_1_screening_output",
        "current_oxford_prediction": False,
        "canonical_regeneration_allowed": False,
        "regeneration_blocked_pending": "docs/multiframe_heating_model_optimisation.md",
    }
    assert len(rows) == 3 * 15


def test_figure_5_4_status_bands_follow_the_two_limits(generated: tuple[dict[str, Path], dict, list[dict[str, str]]]) -> None:
    _, _, rows = generated
    statuses = {
        fluence: [
            row["status"]
            for row in rows
            if float(row["fluence_mw_us"]) == float(fluence)
        ]
        for fluence in (90, 150, 300)
    }
    assert statuses[90] == ["usable"] * 3 + ["quality_limited"] * 7 + ["both_limits"] * 5
    assert statuses[150] == ["usable"] * 6 + ["depletion_limited"] * 2 + ["both_limits"] * 7
    assert statuses[300] == ["usable"] * 3 + ["depletion_limited"] * 4 + ["both_limits"] * 8


def test_figure_5_4_outputs_are_complete(generated: tuple[dict[str, Path], dict, list[dict[str, str]]]) -> None:
    outputs, _, _ = generated
    assert all(path.exists() and path.stat().st_size > 0 for path in outputs.values())
    svg = outputs["svg"].read_text(encoding="utf-8")
    assert "Image number" in svg
    assert "SNR below 8 only" in svg
    assert "loss above 30% only" in svg
    assert "both limits exceeded" in svg


def test_figure_5_4_canonical_regeneration_is_blocked() -> None:
    with pytest.raises(RuntimeError, match="canonical Figure 5.4 is frozen"):
        generate(DEFAULT_CONFIG)


def test_figure_5_4_canonical_metadata_records_retrospective_provenance() -> None:
    metadata = json.loads((CANONICAL_OUTPUT / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["result_status"] == {
        "status": "frozen_version_1_screening_output",
        "current_oxford_prediction": False,
        "canonical_regeneration_allowed": False,
        "regeneration_blocked_pending": "docs/multiframe_heating_model_optimisation.md",
    }

    correction = metadata["provenance_correction"]
    assert correction["generation_base_commit"] == "df00b42ae509921e663cccc0adc0a7a7a240c2a8"
    assert correction["generation_base_commit_status"] == "pre-rewrite local checkpoint identifier"
    assert correction["generation_worktree_dirty"] is True
    assert correction["original_dirty_diff_sealed"] is False
    assert correction["metadata_annotation_only"] is True
    assert len(correction["retrospective_source_snapshot_commit"]) == 40
    assert correction["retrospective_source_snapshot_status"] == (
        "post-generation preserved source/config snapshot; not a clean-generation claim"
    )
