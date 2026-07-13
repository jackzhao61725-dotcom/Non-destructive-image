# Non-destructive Imaging Simulator for Ultracold 166Er Condensates

A notebook-aligned, calibration-aware simulation framework for continuous
non-destructive imaging of an ultracold `166Er` Bose-Einstein condensate. The
repository connects condensate physics, dispersive light-atom interaction,
coherent imaging, camera noise, repeated measurement, and
information-versus-destruction analysis in one reproducible workflow.

The project is intended to answer a practical MSc research question:

> How much useful information can PCI, DGI, or Faraday imaging extract from a
> condensate before photon scattering, heating, reabsorption, and detector
> noise make continued observation destructive or uninformative?

## Scientific Status

The original notebook
`1 calculations revised 2  multishot  6  extended.ipynb` is the historical
computational reference for Version 1 notebook-aligned behaviour. It is not a
final calibrated theory and it is not treated as experimental ground truth.

The maintained implementation is the regression-tested package in
`src/non_destructive_image/`, together with explicit configs, recovery scripts,
tests, and stored baselines. The migration preserved the notebook equations and
FFT conventions; it did not introduce a new microscopic model.

Current status:

- the Version 1 Atomic, Light-Atom, Imaging, Camera, and deterministic
  Multi-shot layers are migrated and tested;
- PCI, DGI, dark-field Faraday, and dual-port Faraday orchestration are present;
- camera binning, Poisson photon noise, and Gaussian read noise are available;
- deterministic one-variable Faraday analysis helpers are available;
- notebook-aligned recovery outputs and numerical consistency audits are
  reproducible;
- absorption/RAI preprocessing helpers exist for future calibration;
- no experimental calibration has yet been applied;
- `kappa_F = 1.0` remains a phenomenological Version 1 placeholder.

The detailed model, equations, assumptions, verified values, and limitations
are collected in
[`docs/simulation_based_physics_report.md`](docs/simulation_based_physics_report.md).

The `dissertation-v1-clean` branch is the compact dissertation working surface.
It retains the simulator, tests, baselines, canonical recovery workflows,
validated results, and thesis-facing documentation while omitting superseded
milestone records and generated bundle archives. The full development history
remains available through Git and `main`.

## Architecture

```text
Explicit configs / future calibration files
                    |
                    v
Atomic Model -> Light-Atom Interaction -> Imaging -> Camera
                                                    |
                                                    v
                                              Multi-shot
                                                    |
                                                    v
                                      Analysis / Calibration
```

| Layer | Main module | Scope |
| --- | --- | --- |
| Atomic Model | `atomic_model.py`, `profiles.py` | Thomas-Fermi state, radii, densities, recoil scales, projected profiles. |
| Light-Atom | `light_atom.py` | Detuning, scalar phase, residual OD, scattering, phenomenological Faraday rotation, reabsorption. |
| Imaging | `fourier.py`, `imaging.py` | Notebook FFT/pupil convention; PCI, DGI, dark-field and dual-port Faraday fields. |
| Camera | `camera.py` | Deterministic binning/count scaling and explicit-RNG Poisson plus read noise. |
| Multi-shot | `multishot.py` | Clean-loss or heating sequence bookkeeping and RMS SNR accumulation. |
| Analysis | `analysis.py` | One-point Faraday proxies, one-dimensional sweeps, and summaries. |
| Calibration | `calibration.py` | Absorption OD preprocessing and cloud observables; no fitting yet. |

## Canonical Version 1 Parameters

The notebook-aligned source is `configs/notebook_v1_defaults.json`. Important
defaults include:

| Quantity | Version 1 value |
| --- | ---: |
| Species / transition | `166Er`, 401 nm |
| Natural linewidth | `Gamma/2pi = 29.5 MHz` |
| Condensate atom number | `N0 = 2.5e4` |
| Scattering length | `72 a0` |
| Trap frequencies `(x,y,z)` | `(293, 14, 233) Hz` |
| Initial cloud temperature | `200 nK` |
| Probe diameter | `24 mm` |
| Existing-arm numerical aperture | `0.080` |
| Camera quantum efficiency | `0.40` |
| Camera read noise | `7 e- rms/pixel` |
| PCI phase plate | amplitude `0.95`, phase `pi/2` |
| DGI stop | `OD = 4` |
| Faraday calibration factor | `kappa_F = 1.0` placeholder |

These are historical Version 1 inputs, not fitted laboratory constants.

## Verified Reference Results

The following numbers are retained because their parameter definitions are
explicit and regression-tested.

At `|Delta|/2pi = 1.5 GHz` along the across-cigar `x` imaging axis:

- dimensionless detuning: `delta = 101.6949`;
- peak scalar phase: `0.2029417 rad`;
- peak phenomenological Faraday rotation: `0.2029417 rad` only because
  `kappa_F = 1.0`;
- continuous clean-loss 30% budget at `P=3.5 mW`, `tau=40 us`: `29.582`
  pulses;
- heating-plus-reabsorption continuous crossing: `13.758` pulses;
- strict accepted full-model frames: `13`.

Two SNR result families must not be mixed:

| Result family at 1.5 GHz | PCI shot only | PCI shot + read | DGI shot only | DGI shot + read |
| --- | ---: | ---: | ---: | ---: |
| Idealised Fig. 3.2: analytical peak pixel, clean loss, identical frames | `72.30` | `70.54` | `36.09` | `24.83` |
| Evolving full model: fixed matched ROI, 13 heating-aware frames | `120.33` | `117.79` | `62.82` | `21.92` |

The first row is a scaling demonstration before NA/PSF and without spatial
summation. The second row uses Fourier propagation, camera binning, a fixed
228-pixel ROI, frame-dependent Thomas-Fermi states, heating, reabsorption, and
RMS accumulation. Their absolute values are not directly comparable.

Legacy values `171.7` and `52/25/24 at 15 us` are deprecated because notebook
defaults mixed exposure times or labels. See
[`docs/thesis_numerical_consistency_correction_report.md`](docs/thesis_numerical_consistency_correction_report.md).

## Installation

Windows PowerShell:

```powershell
python -m venv .venv
.venv\Scripts\activate
python -m pip install --upgrade pip
pip install -r requirements.txt
$env:PYTHONPATH="src;."
$env:PYTHONUTF8="1"
```

Linux/macOS:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
export PYTHONPATH="src:."
export PYTHONUTF8=1
```

## Validation

```powershell
pytest -q
python scripts\validate_notebook_sections.py
```

The tests cover helper behaviour, notebook-derived scalar outputs, PCI/DGI and
Faraday array baselines, camera noise reproducibility, multishot bookkeeping,
figure/result generators, linear-approximation checks, and thesis numerical
consistency.

## Reproduce Results

Run all approved notebook-aligned and dissertation-facing generators:

```powershell
python scripts\run_all_dissertation_figures.py
```

This writes `results/reproducibility_manifest.json` and includes the numerical
consistency audit. Key focused commands are:

```powershell
python scripts\generate_accumulated_snr_invariance_plot.py --config configs\dissertation_plots_v1.json
python scripts\generate_full_multishot_accumulated_snr.py --config configs\dissertation_plots_v1.json
python scripts\audit_thesis_numerical_consistency.py --config configs\thesis_numerical_contract_v1.json
python scripts\generate_dissertation_results.py --config configs\dissertation_results_v1.json
```

The representative Faraday optimisation outputs use a separate placeholder
config (`configs/dissertation_results_v1.json`). They are useful for workflow
testing but are not canonical notebook-cloud predictions.

## Numerical Provenance Contract

Every thesis-facing number must identify its parameter set. A generated
conclusion should record, at minimum:

```text
quantity | value | |Delta|/2pi | P | tau | imaging axis |
normalisation | N_max model | QE/read | repository path
```

Additional mandatory distinctions:

- `per pixel`, `per resolution element`, and `matched ROI` are different
  observables;
- `continuous threshold`, `strict accepted frames`, and sequence-array length
  are different counts;
- `clean loss` is an optimistic analytical bound;
- `heating plus reabsorption` is the current quantitative repeated-imaging
  model;
- exact complex fields are used for images; weak-phase and small-angle formulas
  are interpretive scalings;
- every output must state its calibration status.

The machine-readable contract is
`configs/thesis_numerical_contract_v1.json`.

## Repository Map

```text
src/non_destructive_image/    maintained simulator helpers
tests/                        unit and regression tests
notebook_sections/            exported historical notebook sections
scripts/                      recovery, validation, audit, and result generators
configs/                      explicit notebook, plot, result, and thesis contracts
regression/                   stored numerical baselines
results/                      generated numerical and dissertation outputs
docs/                         physics, architecture, calibration, and thesis notes
deliverables/                 bundle instructions; generated zip is ignored
```

Start with `docs/README.md` for the documentation index and
`docs/reproducibility.md` for the output workflow.

## Current Scientific Boundaries

The repository does not yet provide:

- experimentally fitted atom number, density, magnification, gain, or noise;
- a calibrated Er vector-polarisability model or fitted `kappa_F`;
- a microscopic multi-level Faraday calculation;
- correlated probe noise, reference-image noise, drift, or pixel covariance;
- density-updated reabsorption throughout every sequence;
- a full likelihood or Fisher-information estimator;
- calibrated multi-parameter optimisation;
- droplet, supersolid, mixture, or external GPE state providers.

Accordingly, current outputs are Version 1 representative and uncalibrated.
They support model comparison, implementation validation, and dissertation
methodology, but not a final experimental operating-point claim.

## Future Development

Priority order:

1. **Experimental calibration.** Ingest absorption/RAI data; fit atom number,
   widths, magnification, camera response, and optical-depth scale; validate on
   held-out images.
2. **Faraday calibration.** Replace the placeholder `kappa_F` with a documented
   Er vector-polarisability calculation and/or experimental fit, with
   uncertainty.
3. **Measurement realism.** Add reference-frame noise, technical intensity and
   polarisation noise, drift, covariance, and density-updated reabsorption.
4. **Information metric.** Move from scalar signal proxies to a calibrated
   likelihood or Fisher-information objective under a fixed destruction budget.
5. **Optimisation.** Add constrained multi-parameter searches only after the
   objective and calibration are fixed.
6. **Beyond Thomas-Fermi.** Add phenomenological droplets/supersolids, external
   density maps, and later component-resolved mixtures as additive state
   providers.
7. **Publication.** Select a license, pin a release environment, archive a
   validated snapshot, and add the Zenodo DOI to `CITATION.cff`.

## Citation And License

`CITATION.cff` and `.zenodo.json` are present. A Zenodo DOI has not yet been
issued. No repository license file is currently present, so reuse rights must
not be assumed until a license is selected.
