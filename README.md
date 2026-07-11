# Continuous non-destructive imaging simulations for an ultracold 166Er Bose-Einstein condensate

This repository contains an MSc research simulator for continuous
non-destructive imaging of an ultracold `166Er` Bose-Einstein condensate.

The codebase preserves the historical notebook implementation while adding a
tested, reproducible helper package and dissertation-facing output workflow.
The physical model currently covers:

- Thomas-Fermi condensate construction;
- scalar dispersive phase shifts;
- PCI, DGI, and Faraday imaging;
- deterministic camera binning and stochastic camera noise;
- deterministic multi-shot imaging sequences;
- representative dissertation figure / plot reproduction;
- calibration-ready absorption / RAI observable extraction.

The original notebook remains the scientific reference implementation:

```text
1 calculations revised 2  multishot  6  extended.ipynb
```

The helper package in `src/non_destructive_image/` preserves
notebook-equivalent behaviour. It does not replace the notebook with a new
physics model.

## Repository Status

This is a Version 1 notebook-aligned simulation and dissertation-figure
repository. It is reproducibility-focused and is being prepared for later
citation through an archived Zenodo release.

Current status:

- migrated Version 1 simulator core is implemented and tested;
- notebook-aligned recovery figures are generated and documented;
- representative Version 1 Faraday optimisation outputs are generated;
- a dissertation detuning trade-off physics plot is generated;
- finite-phase and finite-Faraday-rotation approximation validity is audited;
- absorption / RAI calibration-readiness helpers exist;
- release metadata and reproducibility documentation are being prepared.

Current limitations:

- outputs are not experimentally calibrated yet;
- `kappa_F = 1.0` remains a Version 1 phenomenological placeholder;
- Faraday optimisation outputs are representative Version 1 results;
- real absorption / RAI calibration is future work;
- no final calibrated experimental operating-point prediction is claimed.

## Quick Start

Windows PowerShell:

```powershell
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
$env:PYTHONPATH="src;."
$env:PYTHONUTF8="1"
pytest -q
```

Validate that exported notebook sections remain in sync:

```powershell
python scripts\validate_notebook_sections.py
```

## Reproduce Dissertation Figures And Plots

Run the approved reproduction entry point:

```powershell
$env:PYTHONPATH="src;."
$env:PYTHONUTF8="1"
python scripts\run_all_dissertation_figures.py
```

The script runs the approved notebook-aligned recovery and representative
result-generation scripts in a controlled order. It writes:

```text
results/reproducibility_manifest.json
```

The manifest records the git commit, timestamp, scripts run, configs used,
expected outputs, and pending outputs not yet merged into `main`.

## Figure And Plot Outputs

Primary output locations:

- `results/notebook_aligned_recovery/` - canonical notebook-aligned recovery
  figures and stage metadata.
- `results/dissertation_plots_v1/` - approved dissertation physics plots.
- `results/linear_approximation_audit/` - finite-phase and finite-rotation
  approximation-validity audit outputs.
- `results/faraday_optimisation_v1/` - representative Version 1 deterministic
  Faraday optimisation tables and plots.

See `docs/figure_index.md` for the current index of approved figures, plots,
generating scripts, configs, metadata, and intended dissertation use.

## Parameter Provenance Rule

No figure without parameter provenance.

Every dissertation figure or plot should have metadata recording:

- config files used;
- physical parameters;
- optical parameters;
- camera and noise parameters where relevant;
- sweep parameters where relevant;
- units;
- normalisation rule;
- calibration status.

Notebook-aligned recovery outputs should use `configs/notebook_v1_defaults.json`
as the physical/default parameter source. Plot-specific ranges, display
settings, and output paths should live in separate plot configs.

## Validation And Reproducibility

Recommended validation sequence:

```powershell
$env:PYTHONPATH="src;."
$env:PYTHONUTF8="1"
pytest -q
python scripts\validate_notebook_sections.py
python scripts\run_all_dissertation_figures.py
```

For details, see `docs/reproducibility.md`.

## Repository Structure

- `src/non_destructive_image/` - simulator helper package.
- `tests/` - regression and helper tests.
- `notebook_sections/` - exported notebook sections for audit and validation.
- `scripts/` - validation, baseline, bundle, recovery, and result-generation
  scripts.
- `configs/` - notebook defaults and result-generation inputs.
- `results/` - generated dissertation outputs and reproducibility manifest.
- `regression/` - regression baselines used by tests.
- `docs/` - architecture, migration, optimisation, calibration,
  reproducibility, and dissertation-facing documentation.

## Citation

For now, cite the GitHub repository. After the first archived Zenodo release,
cite the Zenodo DOI.

Placeholder citation files are provided:

- `CITATION.cff`
- `.zenodo.json`

DOI will be added after the first archived Zenodo release.

## License

No repository license file is currently present. Choose and add a license before
public archival or reuse. Do not assume reuse rights until a license has been
selected.
