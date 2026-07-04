# Migration Milestone 1 Report — Atomic Model Migration

## Objective

Migrate only the initial Atomic Model calculation in the notebook section exports from notebook-local Thomas-Fermi algebra to the already-extracted helper implementation, while preserving the legacy variable names and downstream interface used by later notebook sections.

## Work Completed

- Identified the primary local Atomic Model calculation in `notebook_sections/02_atomic_model.py`: angular trap frequencies, geometric mean frequency, harmonic-oscillator length, chemical potential, peak density, Thomas-Fermi radii, peak column density, and atom-number consistency check.
- Replaced that local algebra with a call to `build_thomas_fermi_state(...)` from `src/non_destructive_image/atomic_model.py`.
- Preserved all downstream notebook variables: `omega`, `omega_bar`, `a_ho`, `mu`, `T_mu`, `n_peak`, `R`, `n_col`, and `N_check`.
- Left the self-consistent condensate-fraction calculation in `notebook_sections/02_atomic_model.py` unchanged because it depends on `scipy.optimize.brentq` and has not yet been extracted into the helper package.
- Left local Thomas-Fermi recalculations in downstream sections unchanged, especially `tf_state(...)` in `notebook_sections/09_multishot_simulation.py`, because migrating multi-shot code is explicitly outside this milestone.
- Updated `scripts/validate_notebook_sections.py` so migrated section files are still syntax-checked, while exact byte-for-byte notebook export checks remain enforced for unmigrated sections.
- Updated only the Atomic Model entries in `docs/migration_status.md`.

## Files Modified

- `notebook_sections/02_atomic_model.py`
- `scripts/validate_notebook_sections.py`
- `docs/migration_status.md`

## New Files

- `docs/migration_milestone_1_report.md`

## Validation Performed

- PASS — `python3 scripts/validate_notebook_sections.py`
  - Confirmed all unmigrated section exports still match the original notebook-derived text.
  - Confirmed the migrated Atomic Model section is syntactically valid Python.
- PASS — `python3 -m py_compile notebook_sections/02_atomic_model.py scripts/validate_notebook_sections.py src/non_destructive_image/atomic_model.py src/non_destructive_image/profiles.py`
  - Confirmed the migrated section and touched/supporting Python files compile.
- PASS — `python3 -m pytest tests/regression/test_notebook_output_baseline.py -q`
  - Confirmed the stored-output baseline regression test still passes.
- SKIPPED — `PYTHONPATH=src python3 -m pytest -q`
  - The full test command completed with the NumPy-dependent helper tests skipped because `numpy` is not installed in this execution environment.
  - This prevents full numerical re-verification here, but does not indicate a test failure.

## Behaviour Changes

No intended scientific behaviour changes were made.

The migrated section now obtains the same Atomic Model quantities through `build_thomas_fermi_state(...)` instead of spelling out the algebra locally. Legacy notebook variable names and downstream values are preserved by assigning fields from the returned `ThomasFermiState` object back to the original names.

## Remaining Issues

- Full numerical equivalence against a freshly executed notebook could not be verified in this environment because `numpy` and `scipy` are unavailable.
- `notebook_sections/02_atomic_model.py` still contains notebook-local code for the self-consistent condensate-fraction check using `brentq`.
- `notebook_sections/09_multishot_simulation.py` still defines a local `tf_state(...)` function that recomputes Thomas-Fermi quantities during multi-shot evolution; this was intentionally not migrated because multi-shot code is out of scope for this milestone.
- The migrated section now requires the helper package to be importable, e.g. via `PYTHONPATH=src` or an installed package, when running section files directly.

## Recommended Next Milestone

Run Atomic Model numerical equivalence in a complete scientific Python environment with `numpy` and `scipy` installed, then decide whether to extract the remaining condensate-fraction helper or proceed to the next explicitly approved migration target.

## Questions Requiring Review

- Should the self-consistent condensate-fraction calculation be considered part of the Atomic Model helper API in a later milestone?
- Should the multi-shot `tf_state(...)` recalculation be migrated during a future Atomic Model follow-up, or deferred until the Multi-shot Simulation milestone?

## Git Summary

- Current branch: `work`
- Latest commit hash before this report is committed: pending
- Latest commit message before this report is committed: pending
- PR status: pending creation after commit

## Executive Summary

The initial Atomic Model section has been partially migrated from local notebook algebra to the existing helper package. The replacement is deliberately narrow: the helper computes the Thomas-Fermi state, and the section assigns the helper outputs back to the original notebook variable names so later sections see the same interface. No equations, constants, numerical methods, or unrelated modules were changed. Validation confirms syntax and stored-output baseline integrity, but full numerical equivalence remains blocked until the scientific Python stack is available.

# Migration Review

## Notebook Code Removed

The local block in `notebook_sections/02_atomic_model.py` that directly computed `omega`, `omega_bar`, `a_ho`, `mu`, `T_mu`, `n_peak`, `R`, `n_col`, and `N_check` was replaced.

## Helper Functions Used

- `build_thomas_fermi_state(...)` from `src/non_destructive_image/atomic_model.py`

## Notebook Code Still Remaining

- Print statements and legacy variable names in `notebook_sections/02_atomic_model.py` remain unchanged.
- The `scipy.optimize.brentq` condensate-fraction solve remains notebook-local.
- The multi-shot section still contains local Thomas-Fermi recalculation in `tf_state(...)`.
- Thomas-Fermi profile expressions inside PCI/camera/multi-shot demonstrations remain untouched because those sections are outside this milestone.

## What Prevented Further Migration

Further migration would cross the requested boundary into SciPy-based condensate-fraction logic, multi-shot evolution, PCI/camera image formation, or broader Thomas-Fermi profile cleanup. Those changes require separate review and, ideally, a complete numerical baseline generated in an environment with the scientific Python stack installed.

## Review for ChatGPT

This milestone made the smallest safe Atomic Model migration. The only scientific section changed was `notebook_sections/02_atomic_model.py`, where the initial Thomas-Fermi algebra is now delegated to `build_thomas_fermi_state(...)`. The output variable names remain identical to the notebook-local names, so later sections should observe the same values if the helper is numerically equivalent. The helper itself already implements the same equations and is covered by NumPy-dependent equivalence tests, but those tests skip in this environment because `numpy` is unavailable. The validation script was updated because migrated sections no longer match the original notebook text byte-for-byte; it now enforces exact export sync for unmigrated sections and syntax-checks migrated sections.

Careful review should focus on whether this migration boundary is acceptable: the initial Atomic Model is migrated, but the condensate-fraction solve and multi-shot `tf_state(...)` remain notebook-local. This is intentional because migrating those would require either new helper API surface or touching out-of-scope multi-shot code. No code changes are justified beyond this until reviewers approve the partial migration and a complete scientific environment is available for numerical equivalence testing.
