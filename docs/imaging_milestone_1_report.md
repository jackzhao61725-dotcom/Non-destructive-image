# Imaging Milestone 1 Report — Shared PCI/DGI Fourier-Imaging Core

## Objective

Migrate only the shared coherent Fourier-imaging core used by PCI and DGI while preserving the original notebook FFT convention, mask multiplication, reference-field recombination, and intensity calculation. This milestone does not migrate camera noise, shot noise, multi-shot logic, Faraday Jones propagation, scalar phase physics, Atomic Model helpers, Light-Atom helpers, or `kappa_F`.

## Work Completed

- Added `src/non_destructive_image/imaging.py` with one helper, `simulate_fourier_image(...)`.
- The helper implements the shared notebook core:
  - take an object field after atom interaction;
  - form the scattered component as `object_field - 1`;
  - propagate it using the existing FFT/pupil helper;
  - add a caller-supplied reference field;
  - return either `abs(E)**2` intensity or the complex image-plane field.
- Re-exported `simulate_fourier_image(...)` from the package API.
- Updated `notebook_sections/04_pci.py` so the existing `sim_image(...)` wrapper still owns PCI/DGI mode selection, Thomas-Fermi profile construction, and notebook variable names, but delegates the common coherent Fourier core to `simulate_fourier_image(...)`.
- Left `notebook_sections/05_dgi.py` unchanged because it already uses the `sim_image(...)` wrapper defined in section 04.
- Left Faraday field/Jones propagation unchanged because Faraday does not share the PCI/DGI Fourier-plane mask model in this notebook.
- Updated the section validator migrated-section list to include `04_pci.py`.
- Added a helper test that compares `simulate_fourier_image(...)` against the original notebook-equivalent PCI/DGI calculation for a representative object field, pupil, PCI reference, and DGI reference.
- Updated `docs/migration_status.md` and `docs/architecture.md` to document the new shared imaging-core helper and remaining notebook-local imaging responsibilities.

## Files Modified

- `src/non_destructive_image/__init__.py`
- `notebook_sections/04_pci.py`
- `scripts/validate_notebook_sections.py`
- `tests/test_helpers.py`
- `docs/migration_status.md`
- `docs/architecture.md`

## New Files

- `src/non_destructive_image/imaging.py`
- `docs/imaging_milestone_1_report.md`

## Validation Acceptance Status

Scope accepted; architecture accepted; physics boundary accepted; numerical acceptance pending full scientific-environment test execution.

The migrated implementation is structurally accepted for review because `simulate_fourier_image(...)` is limited to the shared PCI/DGI Fourier-core machinery and does not change scalar phase, mask, reference-field, or camera/noise physics. Final numerical acceptance remains pending until the NumPy-backed equivalence tests are run in Google Colab or another complete scientific Python environment.

## Validation Performed

- PASS — `python3 scripts/validate_notebook_sections.py`
  - Confirmed unmigrated section exports still match the notebook-derived text.
  - Confirmed migrated sections `02_atomic_model.py`, `03_light_atom_interaction.py`, `04_pci.py`, and `06_faraday.py` are syntactically valid.
- PASS — `python3 -m py_compile src/non_destructive_image/__init__.py src/non_destructive_image/fourier.py src/non_destructive_image/imaging.py notebook_sections/04_pci.py scripts/validate_notebook_sections.py tests/test_helpers.py`
  - Confirmed touched helper, section, validation, and test files compile.
- PASS — `PYTHONPATH=src python3 -m pytest tests/regression/test_notebook_output_baseline.py -q`
  - Confirmed the stored notebook-output baseline test still passes.
- SKIPPED — `PYTHONPATH=src python3 -m pytest -q`
  - The full test command completed with NumPy-dependent tests skipped because `numpy` is not installed in this environment.
  - The new numerical equivalence test for `simulate_fourier_image(...)` is present and should run normally in Colab or another complete scientific Python environment.
  - The milestone should only be fully closed after those NumPy-backed tests pass in that environment.

## Behaviour Changes

No intended scientific behaviour changes were made.

The notebook-facing `sim_image(...)` wrapper still returns `(intensity, prof)` and still chooses the same PCI, DGI, and fallback reference fields. Only the common propagation/recombination expression was moved into `simulate_fourier_image(...)`.

## Remaining Issues

- Full numerical equivalence against a freshly executed notebook remains blocked in this environment because NumPy/SciPy are unavailable; numerical acceptance is pending Colab/scientific-environment execution.
- PCI-specific phase-plate configuration, SNR analysis, operating-point scans, and plotting remain notebook-local.
- DGI-specific display and analysis remain notebook-local.
- Faraday imaging and Jones propagation remain notebook-local and are deliberately out of scope.
- Camera, shot noise, and multi-shot pathways remain notebook-local except for previously extracted lower-level helpers.

## Recommended Next Milestone

Run the full NumPy-backed helper and regression tests in Google Colab or another complete scientific Python environment. Only after those tests pass should this milestone be considered numerically closed and the next explicitly scoped Imaging migration be reviewed.

## Questions Requiring Review

- Should the next imaging helper split PCI and DGI into separate `simulate_pci(...)` and `simulate_dgi(...)` functions, or keep one mode-selecting wrapper until more numerical baselines exist?
- Should future imaging helpers return complex fields by default, intensities by default, or both via explicit options as this helper does?

## Git Summary

- Current branch: `work`
- Latest commit hash before this report is committed: pending
- Latest commit message before this report is committed: pending
- PR status: pending creation after commit

## Executive Summary

The shared coherent Fourier-imaging core used by PCI and DGI has been migrated into a single helper, `simulate_fourier_image(...)`. The helper is intentionally small: it only handles `object_field - 1`, FFT/pupil propagation, reference-field recombination, and optional intensity conversion. PCI/DGI-specific mode selection remains in the notebook wrapper. This keeps the migration conservative and reversible while creating a stable core for future PCI and DGI migration.

# Migration Review

## What was migrated

The expression equivalent to

```python
Esc = ifft2(fft2(exp(1j * phi_peak_val * prof) - 1) * pupil)
E = reference + Esc
I = abs(E)**2
```

was moved into `simulate_fourier_image(...)`, with the object field and reference field supplied by the caller.

## What was deliberately left unchanged

- Thomas-Fermi profile construction inside `sim_image(...)`.
- PCI reference-field choice `t_p * exp(1j * theta)`.
- DGI reference-field choice `10**(-OD/2)`.
- The fallback reference field `1`.
- Blur-factor calculation and SNR code.
- All plotting and analysis.
- Faraday Jones propagation.
- Camera/noise/multi-shot code.

## Numerical equivalence testing

A new helper test constructs a representative profile, object field, pupil, PCI reference, and DGI reference, then compares `simulate_fourier_image(...)` to the original notebook-equivalent calculation using `np.testing.assert_allclose(...)`. This test is NumPy-dependent and therefore skipped in the current environment, but it is present and ready for a complete scientific environment.

## Discrepancies found

None in available validation. Full numerical comparison could not execute here due to missing NumPy.

## Frozen API confirmation

No frozen Atomic Model or Light-Atom helper APIs were modified. This milestone added a new imaging helper above the frozen lower layers.

## Review for ChatGPT

This milestone migrated only the shared PCI/DGI coherent Fourier-imaging core. It did not migrate PCI-specific analysis, DGI-specific analysis, Faraday Jones propagation, camera noise, shot noise, multi-shot logic, or scalar phase physics. The new helper `simulate_fourier_image(...)` takes an already-formed object field, a pupil/mask, and a mode-specific reference field. This is deliberately minimal and matches the notebook structure: the notebook wrapper still constructs `object_field = exp(1j * phi_peak_val * prof)` and chooses the PCI or DGI reference.

The numerical test compares the helper to the original notebook expression for both PCI and DGI references. It cannot run here because NumPy is missing, but it will run automatically in a complete scientific environment. Reviewers should check that the helper boundary is acceptable: it returns intensity by default, with an option to return the complex image field for future migration. This choice keeps the current notebook `sim_image(...)` behaviour unchanged while allowing later imaging code to work with fields if needed.

Recommendation: approve this migration if the reviewer accepts the minimal helper boundary and the environment-limited test status. The next milestone should run the NumPy-backed test suite, then consider migrating PCI-specific or DGI-specific orchestration separately.

## Validation Support Update

A follow-up validation-support update added `docs/colab_validation.md` and `scripts/run_validation.py`. These do not change physics or helper APIs; they only document and automate how to run the existing validation checks in an environment with NumPy/SciPy available.
