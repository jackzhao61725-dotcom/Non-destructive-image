# Repository architecture

The maintained workflow separates the physical model, configuration, generated
evidence and experimental inputs. The original notebook is retained as a
historical computational reference; new dissertation results are generated from
the tested helper package and explicit configs.

## Data flow

```text
reference condensate and apparatus parameters
                    |
                    v
Thomas-Fermi state and projected column density
                    |
                    v
scalar phase / phenomenological Faraday rotation / scattering
                    |
                    v
finite-aperture PCI, DGI or Faraday image formation
                    |
                    v
camera sampling, photoelectron counts and noise variance
                    |
                    v
single-frame estimator or repeated-exposure state evolution
                    |
                    v
figures, CSV data, metadata and canonical gate
```

## Maintained package

| Module | Responsibility |
| --- | --- |
| `atomic_model.py`, `profiles.py` | reference Thomas-Fermi state, radii, densities and projected profiles |
| `light_atom.py` | detuning, scalar phase, residual absorption estimate, scattering, reabsorption and `theta_F=kappa_F phi` |
| `fourier.py`, `imaging.py` | finite circular pupil and PCI/DGI/dark-field/dual-port fields |
| `camera.py` | object-plane binning, photon scale, Poisson noise and Gaussian read noise |
| `multishot.py` | clean-loss and heating-aware sequence bookkeeping |
| `analysis.py` | lightweight operating-point and one-dimensional sweep helpers |
| `calibration.py` | absorption-image preprocessing utilities for future measured inputs |

The package does not contain a microscopic erbium Faraday calculation or an
experimental state-reconstruction pipeline.

## Configuration roles

- `configs/notebook_v1_defaults.json` preserves historical notebook-aligned
  values for regression.
- `configs/dissertation_v2_dcc3260m.json` is the active provisional apparatus
  and detector contract.
- `configs/dissertation_plots_v1.json` defines the shared screening and
  accumulated-SNR scans.
- `configs/figure_*.json` define individual dissertation figures.
- `configs/thesis_numerical_contract_v1.json` records reporting and legacy
  correction rules.
- `configs/performance_validation_v1.json` freezes the current canonical gate.

A figure generator may refer to a historical config only when the output is
explicitly labelled as notebook-aligned. Active screening outputs must use the
current dissertation config.

## Evidence layers

```text
tests/ + regression/              numerical behaviour and notebook baselines
results/notebook_aligned_recovery historical stage-by-stage recovery
results/dissertation_plots_v1     current dissertation figures and data
results/performance_validation_v1 canonical provenance gate
docs/                             active interpretation and lab hand-off
```

Random camera frames are presentation outputs generated from the same expected
counts and variance model. They do not replace analytic SNR calculation.

## Extension rules

1. Keep every new parameter in a config, not as a hidden script default.
2. Give every dissertation-facing output a metadata record and regression test.
3. Regenerate dependent outputs when sampling, detector or disturbance
   contracts change.
4. Keep `kappa_F` explicit until paired experimental data determine it.
5. Treat measured calibration inputs and held-out validation data as separate
   datasets.
6. Add new condensate state providers without changing the retained
   Thomas-Fermi regression path.
