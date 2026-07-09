# Version 1 Migrated Core Summary

## Current Version 1 Status

The repository now contains a closed Version 1 migrated simulator core for the
continuous non-destructive imaging project.

The migrated core covers:

- Atomic Model;
- Light-Atom Interaction;
- Imaging;
- Camera;
- Stochastic Camera Noise;
- Deterministic Multi-shot Core.

The original notebook remains the authoritative scientific implementation. The
helper package preserves notebook-equivalent behaviour for the migrated core
and supports maintainable, regression-tested development.

No physics redesign was introduced during Version 1 migration.

## Notebook Authority

The authoritative scientific reference remains:

```text
1 calculations revised 2  multishot  6  extended.ipynb
```

The helper modules in `src/non_destructive_image/` are conservative extractions
of notebook-equivalent formulas and orchestration patterns. They are not an
independent replacement model.

The migration strategy was:

- preserve notebook equations and conventions;
- migrate one layer at a time;
- add focused regression tests;
- keep notebook sections available for validation and audit;
- avoid new physics, optimisation, or calibration changes during core migration.

## Layer-By-Layer Summary

### Atomic Model

Implemented helpers cover:

- Thomas-Fermi state construction;
- recoil quantities;
- reusable Thomas-Fermi profile expression.

These helpers keep atomic-state calculations separate from optical interaction
and imaging code.

### Light-Atom Interaction

Implemented helpers cover:

- dimensionless detuning;
- scalar phase shift;
- residual optical depth;
- detected/probe intensity at atoms;
- scattered photons per atom;
- phenomenological Faraday rotation scaling;
- reabsorption fraction.

The current Faraday model remains the notebook's phenomenological convention:

```text
theta_F = kappa_F * phi_peak
kappa_F = 1.0
```

### Imaging

Implemented helpers cover:

- shared Fourier propagation core;
- shared coherent Fourier image recombination;
- PCI orchestration;
- DGI orchestration;
- Faraday orchestration.

PCI/DGI helpers are tested against the PCI/DGI imaging baseline. The Faraday
helper is tested against the Faraday imaging baseline and preserves the
opposite circular-component phase convention and Ex/Ey recombination convention
from the notebook.

### Camera

Implemented helpers cover:

- camera-pixel binning;
- deterministic normalisation;
- deterministic camera pipeline orchestration;
- stochastic camera noise using Poisson photon noise and Gaussian read noise.

The stochastic helper requires an explicit `np.random.Generator`. It does not
create a hidden global RNG and does not hard-code a seed.

### Multi-shot

Implemented helpers cover:

- deterministic multi-shot sequence bookkeeping;
- heating-model condensate update;
- clean-loss model update;
- stop condition at the configured condensate-loss fraction;
- RMS accumulated-SNR convention.

Noisy frame rendering, detuning sweeps, Faraday dual-port frame sequences, and
optimisation maps remain outside the Version 1 migrated core.

## Validation Summary

Current validation state:

```text
pytest -q: 37 passed
notebook section validation: passed
```

Current stored baselines:

```text
regression/baseline/notebook_outputs.json
regression/baseline/imaging/pci_dgi_imaging_baseline_v1.npz
regression/baseline/imaging/faraday_imaging_baseline_v1.npz
```

The baseline files support regression checks while the notebook remains the
scientific authority.

## Remaining Notebook-Local Work

The following work remains notebook-local:

- noisy frame rendering and filmstrips;
- Faraday dual-port frame sequence;
- detuning sweep and operating maps;
- plotting and figure generation;
- optimisation logic;
- broader analysis narrative in `notebook_sections/10_analysis.py`.

These workflows are suitable future migration or report-preparation targets,
but they are not part of the closed Version 1 migrated core.

## MSc Report Relevance

The Version 1 migrated core supports the MSc project by making the simulator
more maintainable and auditable while preserving scientific reproducibility.

The migration separates:

- atomic-state calculations;
- light-atom interaction formulas;
- imaging orchestration;
- camera handling;
- stochastic noise handling;
- deterministic multi-shot bookkeeping.

This separation makes future optimisation more controlled because changes can
be made and tested layer-by-layer. It also supports later information-versus-
destruction analysis by keeping scattering, camera noise, imaging signal, and
multi-shot depletion in distinct components.

For report writing, the migrated core can be described as a reproducibility and
software-engineering contribution: it preserves the original notebook's
scientific behaviour while reducing duplication and exposing tested reusable
simulation layers.

## Recommended Next Directions

Possible next directions, without implementing them here:

- optimisation readiness planning;
- noisy frame-sequence migration;
- report writing and figure preparation;
- RAI-data-based calibration extension;
- beyond-Thomas-Fermi state models such as droplets, supersolids, and mixtures.

The recommended immediate next step is MSc report integration: document the
closed Version 1 core, select figures, and explain how this migration supports
future optimisation of information gained per scattered photon.
