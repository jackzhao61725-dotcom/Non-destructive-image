from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest

from non_destructive_image.reconstruction.studies import (
    build_morphology_study_context,
    load_json,
)
from non_destructive_image.reconstruction.studies.morphology import (
    _representative_specification,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
STUDY_CONFIG = REPO_ROOT / "configs" / "reconstruction_morphology_benchmark_v3.json"


def test_representative_realization_must_exist_in_held_out_ensemble() -> None:
    config = deepcopy(load_json(STUDY_CONFIG))
    config["representative"]["realization_index"] = int(
        config["ensemble"]["held_out_realizations_per_morphology"]
    )
    context = build_morphology_study_context(config)

    with pytest.raises(ValueError, match="held-out ensemble range"):
        _representative_specification(context)
