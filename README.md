# Non-destructive imaging of an ultracold 166Er condensate

- **Status:** dissertation research code and synthetic screening model
- **Calibration boundary:** not an experimentally calibrated apparatus prediction
- **Update trigger:** a change to the active model contract, public entry points,
  installation procedure, retained result status or release information

This repository models non-destructive imaging of an ultracold `166Er`
condensate. It connects a declared atomic state to scalar and Faraday optical
fields, finite-aperture propagation, camera counts, repeated-exposure screening
and physics-informed inference of low-order condensate observables.

Implemented readouts include phase-contrast imaging (PCI), dark-ground imaging
(DGI), dark-field Faraday imaging (DFFI) and dual-port Faraday imaging (DPFI).

Start with the [documentation index](docs/README.md). The
[current Codex hand-off](CODEX_HANDOFF_CURRENT.md) records branch-specific work
and the immediate task; it is an operational document, not a scientific source.

## Scientific scope

The primary inverse result is

```text
q = (A, y_c, z_c, w_major),
```

where `A` is the supported integrated response, `(y_c,z_c)` is the two-
dimensional centroid and `w_major` is the major-axis rms width. Each component
has an uncertainty, data-consistency interval or unsupported status. The fitted
density field is an internal nuisance representation required by the nonlinear
forward operator; the project does not claim super-resolved recovery of a
unique density image.

Current scientific boundaries are:

- the approved isolated-transition, fully polarised axial estimate is the
  signed conversion `kappa_F=-45/91`; the effective apparatus response remains
  uncalibrated;
- the compact `N0=2.5e4` condensate is an optical core surrogate, not the
  thermodynamic initial state of the approved Oxford multiframe replacement;
- effective aperture, detector and sampling values are screening inputs until
  measured on the installed apparatus;
- the 401 nm atom-light treatment is an effective model of an open transition;
- sealed inverse evidence generated with `NA=0.080`, `kappa_F=1` and
  `sigma_r=1.4 e- rms` is historical method-development evidence;
- the stored Figures 5.2 and 5.4 and their `N_dep`/`N_use` values are frozen
  Version 1 screening outputs pending the approved heating replacement; their
  canonical directories are excluded from default regeneration and protected
  against accidental overwrite.

The active parameter authority is
[Simulation reference parameters](docs/simulation_reference_parameters.md).
The minimum approved multiframe replacement is specified in
[Multiframe heating model optimisation](docs/multiframe_heating_model_optimisation.md).
Result status and provenance are maintained in the
[figure and data index](docs/figure_index.md).

## Repository map

```text
src/non_destructive_image/                 forward model and analysis helpers
src/non_destructive_image/reconstruction/  inverse model and observables
configs/                                   active and historical contracts
scripts/                                   generators, recovery and validation
tests/                                     unit and regression tests
results/                                   retained data and figures
docs/                                      active model and evidence documents
notebook_sections/                         historical notebook exports
regression/baseline/                       maintained regression baselines
```

The [notebook export index](notebook_sections/README.md) and
[regression-baseline contract](regression/baseline/README.md) describe the two
historical compatibility surfaces. They do not supersede the maintained package
or active configs.

## Installation

The verified Windows workflow uses a repository-local Python 3.12 environment:

```powershell
$repoRoot = (Resolve-Path -LiteralPath ".").Path
$projectPython = "$repoRoot\.venv\Scripts\python.exe"
& $projectPython -m pip install -e '.[dev]'
& $projectPython -m pip check
& $projectPython -c "import non_destructive_image; print(non_destructive_image.__file__)"
```

The printed package path must lie under this checkout. Do not rely on a Windows
Store `python` alias, a global `pytest` executable or an unrelated editable
clone. Rebuild `.venv` with a verified Python 3.12 interpreter when required.

For Linux or macOS:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
python -m pip check
```

Detailed environment, dependency-order and long-run instructions are maintained
only in [Reproducibility](docs/reproducibility.md).

## Validation and regeneration

Run the complete test suite:

```powershell
& $projectPython -m pytest -q
```

Inspect the maintained generation plan before running it:

```powershell
& $projectPython scripts\run_all_dissertation_figures.py --dry-run
```

Then regenerate the approved default result set with:

```powershell
& $projectPython scripts\run_all_dissertation_figures.py
```

The command writes `results/reproducibility_manifest.json`. Frozen Figures 5.2
and 5.4, historical notebook outputs and long reconstruction studies are
excluded from the default run. Inspect existing sealed artifacts before
launching an expensive regeneration.

## Evidence and provenance

A dissertation-facing numerical result requires an explicit config,
deterministic generator, machine-readable data and metadata, an identified
consumer and a regression check when it supports a claim. Calibration data and
held-out validation data remain separate.

The inverse architecture is described in
[Reconstruction architecture](docs/reconstruction_architecture.md). The frozen
ORCA-Fusion synthetic study, run identifiers and credibility controls are
recorded in the
[sealed reconstruction evidence](docs/reconstruction_orca_v4_evidence_2026_07_21.md).
Those artifacts must not be used as quantitative evidence for the current
signed-`kappa_F`, `NA=0.130`, Ultra quiet screening scenario.

## Experimental hand-off

The laboratory plan measures the quantities that enter the simulations as
provisional inputs: the installed optical path and polarisation response,
magnification and spatial transfer, camera calibration, delivered pulse,
effective Faraday response and net repeated-exposure disturbance. The task order
and acceptance criteria are in the
[experimental measurement plan](docs/experimental_measurement_plan.md).

## Citation and release

Citation and archive metadata are intentionally deferred until the repository
creator and release version are supplied. A Zenodo DOI has not been issued. No
repository licence is currently present, so reuse rights must not be assumed.

Before an archival public release:

1. supply the creator record and intended release version;
2. generate and review the citation and archive metadata;
3. select and add an explicit licence;
4. verify a clean release commit and version tag;
5. run the tests and maintained regeneration workflow;
6. inspect `results/reproducibility_manifest.json` and the retained-result index;
7. archive the tag and add the DOI only after the archive service mints it.
