# Migration Milestone 2 Report — Light–Atom Interaction Layer

## Objective

Migrate only the Light–Atom Interaction layer from notebook-local formulas to helper calls while preserving the original notebook variable names and downstream interfaces. This milestone covers the conversion from atomic density/column density to scalar optical quantities, photon scattering quantities, and the phenomenological Faraday rotation scalar. It does not migrate Fourier propagation, imaging, camera, shot noise, multi-shot evolution, or optimisation.

## Work Completed

- Identified notebook-local Light–Atom calculations in `notebook_sections/03_light_atom_interaction.py`: `delta_of`, `phi_peak`, `od_resonant_equiv`, `intensity_at_atoms`, `N_scatt`, `reabs_frac`, and `theta_F_peak`.
- Identified the duplicated `theta_F_peak` definition in `notebook_sections/06_faraday.py` as a scalar Light–Atom/Faraday-rotation calibration calculation, not a Jones-propagation migration.
- Added helper functions to `src/non_destructive_image/light_atom.py` for the missing extracted formulas:
  - `faraday_rotation_angle(...)`
  - `reabsorption_fraction(...)`
- Re-exported the new helpers from `src/non_destructive_image/__init__.py`.
- Replaced notebook-local formulas in `notebook_sections/03_light_atom_interaction.py` with calls to helper functions while preserving the legacy wrapper names used downstream.
- Replaced the duplicated scalar `theta_F_peak(...)` formula in `notebook_sections/06_faraday.py` with a call to `faraday_rotation_angle(...)` while leaving all Jones propagation and Faraday imaging logic untouched.
- Updated the migration validator so `03_light_atom_interaction.py` and `06_faraday.py` are treated as migrated sections: syntax-checked, but no longer byte-for-byte compared to the original notebook export.
- Extended existing helper tests to cover Faraday rotation scaling and reabsorption fraction formulas.
- Updated only the Light–Atom/Faraday-rotation entries in `docs/migration_status.md`.

## Files Modified

- `src/non_destructive_image/light_atom.py`
- `src/non_destructive_image/__init__.py`
- `notebook_sections/03_light_atom_interaction.py`
- `notebook_sections/06_faraday.py`
- `scripts/validate_notebook_sections.py`
- `tests/test_helpers.py`
- `tests/regression/test_helper_regression_status.py`
- `docs/migration_status.md`

## New Files

- `docs/migration_milestone_2_report.md`

## Validation Performed

- PASS — `python3 scripts/validate_notebook_sections.py`
  - Confirmed unmigrated section exports still match the notebook-derived text.
  - Confirmed migrated sections `02_atomic_model.py`, `03_light_atom_interaction.py`, and `06_faraday.py` are syntactically valid.
- PASS — `python3 -m py_compile src/non_destructive_image/__init__.py src/non_destructive_image/light_atom.py notebook_sections/03_light_atom_interaction.py notebook_sections/06_faraday.py scripts/validate_notebook_sections.py tests/test_helpers.py tests/regression/test_helper_regression_status.py`
  - Confirmed touched Python files compile.
- PASS — `PYTHONPATH=src python3 -m pytest tests/regression/test_notebook_output_baseline.py -q`
  - Confirmed the stored notebook-output baseline test still passes.
- SKIPPED — `PYTHONPATH=src python3 -m pytest -q`
  - The full test command completed with NumPy-dependent tests skipped because `numpy` is not installed in this environment.
  - The new formula-equivalence checks are present, but full numerical execution remains blocked until the scientific Python stack is available.

## Behaviour Changes

No intended scientific behaviour changes were made.

The migrated wrappers preserve the original notebook function names and argument conventions. The changed section code delegates to helpers but returns the same formulas:

- `delta_of(Delta_Hz)` delegates to `dimensionless_detuning(Delta_Hz, Gamma)`.
- `phi_peak(Delta_Hz, n_col_peak)` delegates to `scalar_phase_shift(Delta_Hz, n_col_peak, sigma0, Gamma)`.
- `od_resonant_equiv(Delta_Hz, n_col_peak)` delegates to `residual_optical_depth(Delta_Hz, n_col_peak, sigma0, Gamma)`.
- `intensity_at_atoms(P_mW)` delegates to the helper using `D_probe` and `use_peak_intensity`.
- `N_scatt(...)` delegates to `scattered_photons_per_atom(...)`.
- `reabs_frac(Delta_Hz)` delegates to `reabsorption_fraction(Delta_Hz, n_col, sigma0, Gamma)`.
- `theta_F_peak(...)` delegates to `faraday_rotation_angle(...)` with the notebook's `kappa_F` placeholder.

## Remaining Issues

- Full numerical equivalence against a freshly executed notebook is still blocked in this environment because `numpy` and `scipy` are unavailable.
- Loss-budget wrappers such as `Nmax_loss`, `Nmax_cleanloss`, `Nmax_heating`, and pulse-duration cap calculations remain notebook-local because they combine Light–Atom outputs with destruction/heating model logic outside the core interaction layer.
- `notebook_sections/09_multishot_simulation.py` still directly uses `phi_peak`, `N_scatt`, `intensity_at_atoms`, and a few inline OD/phase expressions. This section was intentionally not modified because multi-shot migration is out of scope.
- PCI, DGI, Faraday field propagation, Jones propagation, camera, shot-noise, and optimisation code remain notebook-local.

## Recommended Next Milestone

Run the Light–Atom formula-equivalence tests and migrated section execution in a complete scientific Python environment with NumPy/SciPy installed, then review whether the next approved migration target should be Fourier/imaging propagation or camera/shot-noise helpers.

## Questions Requiring Review

- Should the future public API use the already-established helper names (`scalar_phase_shift`, `scattered_photons_per_atom`, `faraday_rotation_angle`) or add aliases matching the requested `calculate_*` naming style?
- Should reabsorption remain in the Light–Atom helper module, or should it later move to a destruction/heating-budget module because it feeds the multi-shot loss model?

## Git Summary

- Current branch: `work`
- Latest commit hash before this report is committed: pending
- Latest commit message before this report is committed: pending
- PR status: pending creation after commit

## Executive Summary

The Light–Atom Interaction layer has been partially migrated. The core scalar formulas now live in the helper package and the notebook section wrappers call them while keeping the original notebook API available to downstream code. The Faraday scalar rotation calculation has also been moved to the helper layer, but Jones propagation and imaging remain untouched. The migration preserves equations and constants, but full numerical validation still requires a complete scientific Python environment.

# API Review

## Strengths

- Helper functions keep all physical constants and experimental parameters explicit, avoiding hidden notebook-global dependencies inside the helper package.
- The existing names are specific and traceable to notebook formulas: `dimensionless_detuning`, `scalar_phase_shift`, `residual_optical_depth`, `intensity_at_atoms`, `scattered_photons_per_atom`, `faraday_rotation_angle`, and `reabsorption_fraction`.
- The interface can support PCI and DGI because both consume the scalar phase returned by `scalar_phase_shift` through the preserved `phi_peak` wrapper.
- The interface can support the current Faraday model because `faraday_rotation_angle` takes `kappa_f` explicitly and does not hide the calibration placeholder.
- The interface can support future optimisation because detuning, power, pulse duration, intensity convention, saturation intensity, cross section, and linewidth are explicit inputs.

## Weaknesses

- The API is formula-level rather than object-level; future optimisation may want a parameter bundle to reduce repeated argument lists.
- `faraday_rotation_angle` preserves the current phenomenological `kappa_f * scalar_phase_shift` model but does not resolve the physical calibration of `kappa_f`.
- `reabsorption_fraction` sits at the boundary between Light–Atom interaction and heating/loss modelling, so its long-term module ownership should be reviewed.
- No full notebook execution baseline is available in this environment, so helper correctness remains formula-tested but not fully end-to-end verified here.

## Future Extensibility

The current helper API is sufficient for incremental migration of PCI, DGI, Faraday, and optimisation workflows without changing physics, provided those future modules continue to receive explicit constants/parameters and preserve the notebook's scalar-phase and phenomenological-Faraday definitions. A future ergonomics milestone may introduce parameter dataclasses or `calculate_*` aliases, but that should not happen until numerical baselines are available.

## Architecture Review

The repository has a cleaner separation between Atomic Model and Light–Atom Interaction, but it is not complete.

Current separation:

```text
Atomic Model helper
  -> produces Thomas-Fermi state / column densities
  -> notebook variables preserve n_col and related legacy names
Light-Atom helper
  -> consumes column density + optical constants
  -> produces delta, scalar phase, residual OD, scattering, reabsorption, Faraday rotation
Notebook imaging sections
  -> still consume wrapper names such as phi_peak, N_scatt, theta_F_peak
```

Remaining coupling:

- Notebook sections still rely on global names such as `sigma0`, `Gamma`, `D_probe`, `Isat`, `tau`, `use_peak_intensity`, `n_col`, and `kappa_F`.
- The helper package is decoupled from notebook globals, but the notebook wrappers still bind helper arguments to those globals.
- Multi-shot logic still contains local phase/OD expressions and should not be considered migrated.

## Review for ChatGPT

This milestone migrated the core Light–Atom scalar formulas without touching imaging or propagation. The main scientific section changed was `notebook_sections/03_light_atom_interaction.py`, where local formulas are now thin wrappers around helper calls. The duplicated `theta_F_peak(...)` in `notebook_sections/06_faraday.py` was also changed because it is a scalar Light–Atom/Faraday-rotation calibration formula; the Jones field propagation and Faraday imaging code were intentionally left untouched.

The helper module now contains the existing detuning, scalar phase, OD, intensity, and scattering functions plus new direct extractions for `faraday_rotation_angle(...)` and `reabsorption_fraction(...)`. These additions do not introduce new physics: `faraday_rotation_angle` implements the notebook's documented `kappa_F * phi_peak` convention, and `reabsorption_fraction` implements the existing `mean(1 - exp(-OD))` calculation.

Numerical verification is still limited by the execution environment. The formula-equivalence tests have been extended but are skipped when NumPy is unavailable. The stored-output regression test still passes, and all touched files compile. A reviewer should verify that migrating `theta_F_peak` in `06_faraday.py` is acceptable within this milestone, because it touches the Faraday section but only for the scalar rotation calculation explicitly included in scope. The migration should be approved if the reviewer accepts that boundary. The next milestone should not begin until these tests are run in a complete scientific Python environment or the reviewer explicitly accepts the environment limitation.
