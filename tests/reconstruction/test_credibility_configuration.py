from __future__ import annotations

from pathlib import Path

from non_destructive_image.reconstruction.studies import (
    build_morphology_study_context,
    load_json,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
CREDIBILITY_CONFIG = (
    REPO_ROOT / "configs" / "reconstruction_credibility_v2_orca_fusion_m10.json"
)


def test_orca_credibility_config_uses_frozen_v4_candidates() -> None:
    credibility = load_json(CREDIBILITY_CONFIG)
    benchmark_path = REPO_ROOT / credibility["source_benchmark_config"]
    benchmark = load_json(benchmark_path)
    context = build_morphology_study_context(benchmark)
    candidate_labels = {candidate.label for candidate in context.candidates}

    assert credibility["frozen_candidates"] == {
        "dual_port": "resolution_matched_17x5__curvature_30_um2",
        "dark_field": "resolution_matched_17x5__curvature_30_um2",
    }
    assert set(credibility["frozen_candidates"].values()) <= candidate_labels
    assert credibility["representative"]["morphology_name"] in {
        case.name for case in context.morphology_split.held_out
    }
    assert credibility["representative"]["fluence_mw_us"] in benchmark["ensemble"][
        "held_out_fluence_mw_us"
    ]

    for check_path in credibility["curvature_range_evidence"].values():
        check = load_json(REPO_ROOT / check_path)
        assert check["source_benchmark_config"] == credibility[
            "source_benchmark_config"
        ]
        assert check["reference_weight_um2"] == 30.0
