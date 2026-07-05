# Migration Status

This document tracks verification status before any further migration from notebook-local code to the helper package. The original notebook remains the authoritative reference implementation.

## Dependency graph

```text
Original Notebook
  |
  |  stored outputs + notebook JSON
  v
Regression Baseline
  |
  |  validates output hashes/text and section export sync
  v
Notebook Section Exports
  |
  |  expose duplicated patterns for mechanical replacement
  v
Helper Modules
  |
  |  equivalence tests compare helpers to original notebook formulas
  v
Regression Tests
  |
  |  gate future replacements and section migration
  v
Future Refactoring
```

Current direct dependencies:

- `notebook_sections/` still contains notebook-local code and notebook-global variables.
- `src/non_destructive_image/` does not depend on notebook globals; helper inputs are explicit.
- `tests/regression/test_notebook_output_baseline.py` depends on stored notebook outputs and `regression/baseline/notebook_outputs.json`.
- `tests/regression/test_helper_regression_status.py` depends on `numpy`; it is skipped when `numpy` is unavailable.

## Helper validation status

| Area | Module / location | Status | Explanation |
|---|---|---|---|
| Atomic model | `src/non_destructive_image/atomic_model.py` | PARTIAL | Initial TF state algebra is now called by `notebook_sections/02_atomic_model.py` and covered by formula-equivalence tests when `numpy` is available; full executed-notebook numerical regression is still blocked by missing scientific dependencies. |
| Thomas-Fermi profile | `src/non_destructive_image/profiles.py` | PARTIAL | 2D profile expression extracted and tested against the original expression when `numpy` is available. |
| Fourier propagation | `src/non_destructive_image/fourier.py` | PARTIAL | FFT/pupil operation extracted and tested against the original FFT expression when `numpy` is available. |
| Camera model | `src/non_destructive_image/camera.py` | PARTIAL | Binning, noise, and normalisation recipes extracted and tested against original expressions when `numpy` is available; stochastic image-level regression still needs a full dependency environment. |
| Light-atom interaction | `src/non_destructive_image/light_atom.py` | PARTIAL | Detuning, scalar phase, residual OD, intensity, scattering, Faraday-rotation scaling, and reabsorption helpers are extracted; `notebook_sections/03_light_atom_interaction.py` and the Faraday rotation definition in `06_faraday.py` now call these helpers while preserving notebook wrapper names. |
| PCI | `notebook_sections/04_pci.py`, `src/non_destructive_image/imaging.py` | PARTIAL | Shared PCI/DGI Fourier-core propagation now uses `simulate_fourier_image`; PCI-specific phase-plate setup, SNR, plotting, and operating-point analysis remain notebook-local. |
| DGI | `notebook_sections/04_pci.py`, `notebook_sections/05_dgi.py`, `src/non_destructive_image/imaging.py` | PARTIAL | DGI uses the migrated shared Fourier-core via the existing `sim_image` wrapper; DGI-specific demonstrations and downstream plotting remain notebook-local. |
| Faraday | Notebook sections only | NOT VERIFIED | Faraday is exported and notebook outputs are baselined, but no Faraday helper/module has been migrated yet. |
| Multi-shot simulation | Notebook sections only | NOT VERIFIED | Multi-shot evolution remains notebook-local; only stored notebook outputs are baselined. |

## Notebook section migration status

| Target section | Current source | Extracted? | Validated? | Still depends on notebook code? | Ready to migrate completely? | Notes |
|---|---|---:|---:|---:|---:|---|
| 00 Imports | `notebook_sections/00_imports.py` | No | Export validated | Yes | No | Imports and plotting setup remain notebook-style. |
| 01 Parameters | `notebook_sections/01_parameters.py` | No | Export validated | Yes | No | Parameter values should remain centralised and unchanged. |
| 02 Atomic Model | `notebook_sections/02_atomic_model.py`, `atomic_model.py`, `profiles.py` | Partial | Partial | Partial | Partial | Initial TF state values now come from `build_thomas_fermi_state`; the section still performs the self-consistent condensate-fraction solve locally and downstream sections still contain local TF recalculations. |
| 03 Light-Atom Interaction | `notebook_sections/03_light_atom_interaction.py`, `light_atom.py` | Partial | Partial | Partial | Partial | Core detuning, phase, OD, scattering, reabsorption, and Faraday-rotation scalar calculations now call helper functions; plotting, loss-budget wrappers, and downstream demonstrations remain notebook-local. |
| 04 PCI | `notebook_sections/04_pci.py`, `imaging.py`, `fourier.py`, `camera.py` | Partial | Partial | Partial | Partial | The shared coherent Fourier-imaging core now calls `simulate_fourier_image`; PCI-specific setup, SNR, plots, and analysis remain notebook-local. |
| 05 DGI | `notebook_sections/05_dgi.py`, `notebook_sections/04_pci.py`, `imaging.py`, `fourier.py`, `camera.py` | Partial | Partial | Partial | Partial | DGI shares the migrated `sim_image` Fourier-core path from section 04; DGI displays and analysis remain notebook-local. |
| 06 Faraday | `notebook_sections/06_faraday.py`, `light_atom.py` | Partial | Partial | Yes | No | The scalar `theta_F_peak` calculation now calls the light-atom helper, but Jones propagation, dual-port detection, and Faraday imaging remain notebook-local. |
| 07 Camera | `notebook_sections/07_camera.py`, `camera.py` | Partial | Partial | Yes | No | Camera helper exists; section still mixes demonstrations with helper definitions. |
| 08 Shot Noise | `notebook_sections/08_shot_noise.py` | No | Baseline only | Yes | No | SNR helpers remain notebook-local. |
| 09 Multi-shot Simulation | `notebook_sections/09_multishot_simulation.py` | No | Baseline only | Yes | No | Heating/loss sequence code remains notebook-local. |
| 10 Analysis | `notebook_sections/10_analysis.py` | No | Baseline only | Yes | No | Plotting and pedagogical analysis remain notebook-local. |

## Regression baseline

A lightweight stored-output baseline exists at `regression/baseline/notebook_outputs.json`.

What it captures:

- Stream text from executed notebook cells, including important scalar outputs such as recoil quantities, Thomas-Fermi radii, chemical potential, peak density, column density, phase/scattering summaries, and printed comparison tables.
- Hashes and byte sizes for rich display payloads already stored in the notebook, including generated figure outputs.

What it does not yet capture:

- Raw intermediate arrays such as full PCI/DGI/Faraday image arrays.
- A freshly executed notebook run in the current environment.
- `.npz` array baselines for helper-to-notebook numerical comparisons.

## Current extraction percentage

Conservative status estimate:

- Safely extracted into helper modules: 5 of 11 target sections have at least partial helper coverage (`02`, `03`, `04`, `05`, `07`) = about 45% structurally covered.
- Fully migrated and verified sections: 0 of 11 = 0%.
- Full notebook replacement readiness: 0%; the notebook remains authoritative.

This percentage is intentionally conservative. The helper modules are useful and ready for incremental migration, but no section should be considered fully migrated until it has been wired to helpers and checked against a numerical baseline in a complete scientific Python environment.

## Next migration gate

Before migrating any section completely:

1. Install the dependency stack (`numpy`, `scipy`, `matplotlib`, `pytest`, and notebook execution tools as needed).
2. Execute the original notebook to produce fresh scalar and array baselines.
3. Add `.npz` baselines for representative arrays: phase maps, PCI/DGI/Faraday images, camera images, and multi-shot sequences.
4. Replace one duplicated implementation with one helper.
5. Compare the replacement against the original baseline using tolerances appropriate to floating-point and stochastic calculations.
