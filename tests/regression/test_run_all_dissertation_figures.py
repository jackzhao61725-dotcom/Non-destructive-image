from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path


SCRIPT_PATH = Path("scripts/run_all_dissertation_figures.py")


def _load_module():
    spec = importlib.util.spec_from_file_location("run_all_dissertation_figures", SCRIPT_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_run_all_plan_contains_only_approved_main_generators() -> None:
    module = _load_module()
    scripts = {step.script for step in module.RUN_STEPS}

    assert "scripts/validate_notebook_sections.py" in scripts
    assert "scripts/recover_notebook_condensate_stage.py" in scripts
    assert "scripts/recover_notebook_phase_stage.py" in scripts
    assert "scripts/recover_notebook_pci_stage.py" in scripts
    assert "scripts/recover_notebook_dgi_stage.py" in scripts
    assert "scripts/recover_notebook_faraday_stage.py" in scripts
    assert "scripts/recover_notebook_camera_stage.py" in scripts
    assert "scripts/recover_notebook_noisy_camera_stage.py" in scripts
    assert "scripts/recover_notebook_multishot_stage.py" in scripts
    assert "scripts/recover_notebook_noisy_multishot_filmstrip.py" in scripts
    assert "scripts/generate_condensate_three_view.py" in scripts
    assert "scripts/generate_dissertation_results.py" in scripts
    assert "scripts/generate_detuning_tradeoff_plot.py" in scripts
    assert "scripts/audit_linear_approximation_validity.py" in scripts


def test_run_all_has_no_pending_approved_outputs() -> None:
    module = _load_module()

    assert module.PENDING_ITEMS == ()


def test_run_all_dry_run_succeeds_without_manifest(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.json"
    result = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "--dry-run", "--manifest", str(manifest)],
        check=True,
        text=True,
        capture_output=True,
    )

    assert "recover_condensate_stage" in result.stdout
    assert "generate_detuning_tradeoff_plot" in result.stdout
    assert "audit_linear_approximation_validity" in result.stdout
    assert not manifest.exists()
