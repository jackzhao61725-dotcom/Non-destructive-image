# Milestone Report

## Objective

Establish a numerical-regression and validation layer for the original ultracold-atom imaging notebook before any further refactoring or migration occurs. The original notebook remains the authoritative scientific implementation.

## Work Completed

- Added a stored-output regression baseline at `regression/baseline/notebook_outputs.json`.
- Added `scripts/extract_notebook_output_baseline.py` to regenerate the stored-output baseline from the reference notebook JSON.
- Added regression tests under `tests/regression/`:
  - `test_notebook_output_baseline.py` checks stored notebook outputs against the checked-in baseline.
  - `test_helper_regression_status.py` compares selected helper outputs against original notebook formulas when `numpy` is available.
- Added `docs/migration_status.md` with:
  - helper validation statuses,
  - dependency graph,
  - notebook-section migration readiness,
  - extraction percentage estimate,
  - next migration gate.
- Reviewed the current validation state and documented limitations explicitly.

## Files Modified

- `docs/migration_status.md`
- `regression/baseline/notebook_outputs.json`
- `scripts/extract_notebook_output_baseline.py`
- `tests/regression/test_helper_regression_status.py`
- `tests/regression/test_notebook_output_baseline.py`
- `docs/milestone_2_report.md`

## New Files

- `docs/migration_status.md`
- `docs/milestone_2_report.md`
- `regression/baseline/notebook_outputs.json`
- `scripts/extract_notebook_output_baseline.py`
- `tests/regression/test_helper_regression_status.py`
- `tests/regression/test_notebook_output_baseline.py`

## Validation Performed

PASS — `python3 scripts/extract_notebook_output_baseline.py`

- Regenerated `regression/baseline/notebook_outputs.json` from outputs already stored in the notebook.

PASS — `python3 scripts/validate_notebook_sections.py`

- Confirmed exported notebook section files remain in sync with notebook JSON and are syntactically valid.

PASS — `python3 -m py_compile src/non_destructive_image/*.py tests/test_helpers.py tests/regression/*.py scripts/*.py`

- Confirmed helper modules, tests, and scripts are Python-syntax valid.

PASS — `python3 -m pytest tests/regression/test_notebook_output_baseline.py -q`

- Confirmed current stored notebook outputs match the checked-in stored-output baseline.

SKIPPED — `python3 -m pytest -q`

- Overall pytest command completed with one passing stored-output baseline test and two skipped numpy-dependent tests.
- Skips occurred because `numpy` is not installed in this execution environment.

SKIPPED — full fresh notebook execution

- The full original notebook was not freshly executed in this environment because required scientific dependencies such as `numpy` are unavailable.
- The current baseline is therefore a stored-output baseline extracted from the already executed notebook, not a fresh runtime baseline.

## Behaviour Changes

No observable scientific behaviour was changed.

- Notebook code was not modified.
- Equations were not changed.
- Algorithms were not optimised or simplified.
- Parameter values were not changed.
- New files are validation, baseline, reporting, and test infrastructure only.

## Remaining Issues

- Full numerical array baselines are not yet available for PCI, DGI, Faraday, shot-noise images, or multi-shot arrays.
- Full helper-vs-notebook numerical tests require `numpy` and, for complete notebook execution, likely `scipy`, `matplotlib`, and notebook execution tooling.
- The current stored-output baseline captures scalar text and rich-output hashes, but not raw arrays.
- Helper modules are not yet wired into notebook sections.
- PCI, DGI, Faraday, shot-noise, and multi-shot logic remain notebook-local and are not verified as migrated modules.

## Recommended Next Milestone

Create a full executed numerical baseline in an environment with the scientific Python stack installed, saving representative raw arrays as `.npz` files for phase maps, PCI/DGI/Faraday images, camera images, shot-noise outputs, and multi-shot evolution.

## Questions Requiring Review

None.

## Git Summary

- Current branch: `work`
- Latest commit hash at report creation: `fff5b82`
- Latest commit message at report creation: `Add regression baseline and migration status`
- PR title: `Add regression baseline and migration status`
- PR status: Created for review

## Executive Summary

The repository now has a Milestone 2 verification layer, but it should be interpreted conservatively. The original notebook remains the only authoritative runnable scientific workflow. A stored-output baseline has been extracted from the already executed notebook, and regression tests now detect drift in those stored outputs. Helper modules for atomic model, Thomas-Fermi profiles, Fourier propagation, camera handling, and light-atom formulas exist and have formula-equivalence tests, but the numpy-dependent tests could not run in this environment. The main remaining work is to execute the original notebook in a complete scientific Python environment and save raw numerical array baselines before migrating any notebook section to helper modules.

## Review Request

Please review the Milestone 2 validation layer before any further refactoring. Confirm that the stored-output baseline is an acceptable interim baseline, and pay close attention to the documented limitation that a fresh numerical notebook execution was not possible in this environment. The highest risk for the next milestone is migrating notebook code before raw array baselines exist for PCI, DGI, Faraday, shot-noise, and multi-shot outputs. The next milestone should therefore focus only on creating full `.npz` numerical baselines in a complete scientific Python environment, not on architectural cleanup.
