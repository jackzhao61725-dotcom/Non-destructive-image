# Non-destructive imaging of an ultracold 166Er condensate

This repository contains the reproducible forward model used to screen
non-destructive imaging conditions for an ultracold `166Er` condensate. It
connects a reference Thomas-Fermi state to scalar and Faraday image formation,
finite-aperture propagation, camera detection and repeated-exposure evolution.

The model supports the dissertation and experimental commissioning. It is not
yet an experimentally calibrated prediction.

## Current status

The Dissertation v2 snapshot is maintained through the tested package in
`src/non_destructive_image/`, together with explicit configs, generators and
stored result metadata. The original notebook remains a historical
computational reference.

Implemented readouts:

- phase-contrast imaging (PCI);
- dark-ground imaging (DGI);
- dark-field Faraday imaging;
- dual-port Faraday imaging.

Current scientific boundaries:

- `kappa_F = 1` is an uncalibrated structural placeholder;
- the contact-interaction Thomas-Fermi state omits dipolar mean-field effects;
- the 401 nm light-atom treatment is an effective two-level approximation to an
  open transition;
- the repeated-exposure model assumes rethermalisation and a fixed
  initial-density reabsorption fraction;
- detector, aperture and magnification values remain screening inputs until
  measured on the apparatus.

## Active screening contract

The authoritative parameter record is
[`docs/simulation_reference_parameters.md`](docs/simulation_reference_parameters.md).

| Quantity | Active value | Status |
| --- | ---: | --- |
| Initial condensate population | `2.5e4` | reference input |
| Reference detuning | `|Delta|/2pi = 1.5 GHz` | scan reference |
| Reference power / exposure | `1.0 mW / 90 us` | one division of the fluence |
| Reference fluence | `90 mW us` | screening reference |
| Numerical aperture | `0.080` | provisional |
| Magnification | `4` | provisional |
| Physical camera pixel | `5.86 um` | DCC3260M hardware value |
| Object-plane pixel | `1.465 um` | provisional `M=4` sampling |
| Quantum efficiency | `0.60` | provisional |
| Read noise | `3 e- rms` per pixel and readout | provisional |
| Recoil-energy convention | `2 E_rec` per scattering cycle | model choice |
| Condensate-loss limit | `30%` | screening criterion |
| Faraday coefficient | `kappa_F = 1` | uncalibrated placeholder |

The historical notebook defaults remain in
`configs/notebook_v1_defaults.json`. They must not be mixed with the active
screening contract in `configs/dissertation_v2_dcc3260m.json`.

## Canonical performance gate

At the reference point, the heating-plus-initial-density-reabsorption model
accepts 10 strict frames. The post-sequence condensate loss is `0.28060`; the
next pulse would raise it to `0.30829` and is excluded.

The current accumulated matched-template SNR values use the same 228-pixel
object-space ROI and the evolving 10-frame sequence:

| Readout | Shot noise only | Shot + read noise |
| --- | ---: | ---: |
| PCI | 104.10 | 103.67 |
| DGI | 54.36 | 34.30 |
| Faraday dark-field | 55.78 | 34.71 |
| Faraday dual-port | 107.48 | 106.56 |

These values are estimator-specific screening results. In particular, they are
not interchangeable with the single-frame central-`3x3` SNR used in the
Chapter 5 image-quality figures. The Faraday rows are not calibrated absolute
predictions.

Machine-readable evidence is stored in
`results/performance_validation_v1/`, with the contract in
`configs/performance_validation_v1.json`.

## Repository map

```text
src/non_destructive_image/    maintained physical and imaging helpers
configs/                      active and historical parameter contracts
scripts/                      generators, audits and validation entry points
tests/                        unit and regression tests
regression/                   stored notebook-derived numerical baselines
results/                      maintained generated data and dissertation figures
docs/                         active model, provenance and laboratory documents
notebook_sections/            exported historical notebook sections
```

Start with:

- [documentation index](docs/README.md);
- [active parameter contract](docs/simulation_reference_parameters.md);
- [figure and data index](docs/figure_index.md);
- [experimental measurement plan](docs/experimental_measurement_plan.md);
- [reproducibility guide](docs/reproducibility.md).

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

## Validation and reproduction

Run the test suite:

```powershell
pytest -q
```

Check the maintained canonical outputs without regenerating them:

```powershell
python scripts\run_performance_validation.py --skip-prerequisites
```

Regenerate the approved result set and write
`results/reproducibility_manifest.json`:

```powershell
python scripts\run_all_dissertation_figures.py
```

Use `--dry-run` to inspect the generation plan first.

## Provenance rules

Every dissertation-facing numerical result must identify:

```text
quantity | value | |Delta|/2pi | P | tau | imaging axis |
normalisation | estimator/ROI | sequence model | QE/read | config | output
```

The following distinctions are mandatory:

- single-frame `SNR_3x3` and accumulated matched-template SNR are different
  observables;
- continuous clean-loss budgets and strict accepted-frame counts are different
  quantities;
- random camera frames illustrate the noise model, while reported SNR values
  are calculated from the expected counts and analytic variance;
- `kappa_F=1` Faraday values are structural comparisons, not calibrated
  erbium predictions.

## Experimental hand-off

Commissioning should first measure the quantities that currently enter as
provisional inputs: optical path and polarisation behaviour, magnification and
effective aperture, camera offset/gain/read noise, delivered power and pulse
timing, the effective Faraday response and the net disturbance of repeated
imaging. The detailed task order and acceptance criteria are recorded in
`docs/experimental_measurement_plan.md`.

Measured values used to set model parameters constitute calibration. A
held-out operating condition is required before agreement can be described as
experimental validation.

## Citation and licence

`CITATION.cff` and `.zenodo.json` are present. A Zenodo DOI has not yet been
issued. No repository licence file is currently present, so reuse rights must
not be assumed until a licence is selected.
