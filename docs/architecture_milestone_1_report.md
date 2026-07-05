# Architecture Milestone 1 Report — API Review & Architecture Freeze

## Objective

Review the current repository architecture after Atomic Model and Light-Atom Interaction migrations, determine whether the helper API is stable enough for future Imaging migration, and document the architecture without changing physics, helper implementations, or notebook migration state.

## Work Completed

- Reviewed the current logical layers: imports/runtime setup, parameters, Atomic Model, Light-Atom Interaction, Imaging, Camera, Shot Noise, Multi-shot Simulation, and Analysis.
- Reviewed every current public helper exported through the helper package.
- Reviewed the scientific data flow from atomic state through interaction quantities, imaging, camera/noise, multi-shot evolution, and analysis.
- Reviewed future extensibility for PCI, DGI, Faraday, and optimisation over detuning, intensity, and exposure time.
- Reviewed module boundaries and documented functions that are acceptable but potentially near future boundaries, especially `reabsorption_fraction(...)` and `thomas_fermi_profile_2d(...)`.
- Rewrote `docs/architecture.md` as the architecture freeze document with current architecture, recommended stable architecture, block diagrams, helper API review, data-flow review, coupling risks, and freeze recommendation.

## Files Modified

- `docs/architecture.md`

## New Files

- `docs/architecture_milestone_1_report.md`

## Validation Performed

- PASS — Documentation-only review; no executable scientific code was changed.
- PASS — `python3 scripts/validate_notebook_sections.py`
  - Confirmed the existing section-validation workflow still passes after documentation updates.
- PASS — `PYTHONPATH=src:. python3 -m pytest tests/regression/test_notebook_output_baseline.py -q`
  - Confirmed the stored notebook-output baseline test still passes.
- SKIPPED — Full numerical equivalence remains blocked by missing NumPy/SciPy in this environment and was not required for this documentation-only milestone.

## Behaviour Changes

None.

No helper implementations, notebook sections, equations, APIs, tests, scripts, or physics calculations were changed in this milestone.

## Remaining Issues

- Full numerical scientific baselines are still unavailable in this environment.
- Notebook sections still depend on globals and execution order.
- Imaging code is not yet migrated and still mixes simulation, propagation, plotting, and analysis.
- `kappa_f` remains a documented Faraday calibration placeholder, not a derived physical calibration.
- Future optimisation may benefit from parameter objects, but adding them now would be premature.

## Recommended Next Milestone

Begin the Imaging Layer migration only after review approval of this architecture freeze. The next migration should be explicitly scoped, preferably PCI/DGI shared Fourier-imaging structure, without changing the frozen lower-level Atomic Model and Light-Atom APIs.

## Questions Requiring Review

- Should future imaging helpers be grouped in one `imaging.py` module or split into `pci.py`, `dgi.py`, and `faraday.py`?
- Should optimisation ergonomics later introduce parameter dataclasses while keeping the current formula helpers stable?

## Git Summary

- Current branch: `work`
- Latest commit hash before this report is committed: pending
- Latest commit message before this report is committed: pending
- PR status: pending creation after commit

## Executive Summary

The architecture is ready to freeze for the next migration phase. The helper API has clear physical layers: Atomic Model produces state and column densities; Light-Atom Interaction consumes column density and optical constants; Fourier and Camera helpers are independent lower-level utilities; future Imaging helpers should sit above those layers. The main remaining risk is not API design but scientific validation: full numerical baselines still need a complete scientific Python environment.

# Architecture Review

## Strengths

- The helper package keeps physics responsibilities separated: atomic state, profiles, light-atom formulas, Fourier propagation, and camera operations are distinct.
- Current public helpers use explicit inputs rather than notebook globals.
- The notebook wrappers preserve legacy names, allowing incremental migration without breaking downstream sections.
- Faraday calibration remains explicit through `kappa_f`, so the placeholder is not hidden in the architecture.
- The architecture supports future PCI/DGI/Faraday migration by adding imaging-level helpers above the current formula helpers.

## Weaknesses

- Notebook section exports still rely on global variables and execution order.
- Plotting, analysis, and physics calculations remain interleaved in several sections.
- Full numerical baselines are not available in this environment.
- The API is formula-level and may become verbose for optimisation workflows.

## Unnecessary Coupling

- `03_light_atom_interaction.py` binds helper calls to notebook globals through wrapper functions.
- `06_faraday.py` still mixes scalar Faraday rotation setup, field propagation, and plotting.
- `09_multishot_simulation.py` still contains direct phase/OD calculations and should not be considered migrated.

## Future Risks

- Renaming helper APIs during Imaging migration would destabilise already-reviewed layers.
- Moving `reabsorption_fraction(...)` prematurely could create churn without changing physics.
- Treating `kappa_f` as calibrated rather than a placeholder would be a physics change and should require a separate review.
- Adding optimisation dataclasses before numerical baselines exist could obscure equivalence checks.

# Architecture Freeze Recommendation

**APPROVE ARCHITECTURE**

The current public helper API is sufficiently clean and stable for future migration to proceed without changing public APIs. Remaining issues are validation and migration-scope concerns, not architecture blockers. The freeze should apply to the current low-level helper names and module responsibilities; future imaging helpers should be added above these layers rather than modifying them.

## Review for ChatGPT

This milestone was documentation and architecture review only. No notebook code, helper code, equations, APIs, or physics calculations were changed. The updated `docs/architecture.md` now defines the current architecture, logical layers, helper responsibilities, public API review, data flow, coupling risks, future extensibility, recommended stable architecture, and freeze recommendation.

The recommended architecture is: Atomic Model helpers produce Thomas-Fermi state and column densities; Light-Atom helpers consume column density and optical constants to produce scalar phase, OD, scattering, reabsorption, and phenomenological Faraday rotation; Fourier and Camera helpers remain lower-level utilities; future Imaging helpers should sit above these layers and orchestrate PCI, DGI, and Faraday simulations.

The API is judged stable enough for Imaging migration. The main limitations are known and documented: no full numerical baselines in this environment, notebook wrappers still bind to globals, `reabsorption_fraction(...)` is near the heating/loss boundary, and future optimisation may want parameter objects. None of these require changing current public helper functions before proceeding.

Recommendation: approve the architecture freeze, then begin a separately scoped Imaging Layer migration after review. The next milestone should not rename or redesign the Atomic Model or Light-Atom Interaction APIs.
