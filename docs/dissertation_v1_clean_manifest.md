# Dissertation V1 Clean Branch Manifest

## Purpose

`dissertation-v1-clean` is a temporary, usable dissertation working branch. It
starts from validated `main` commit
`319d96f5c122ed2da843cd7eede82bdd2299436d` and keeps the material needed to:

- run and test the simulator;
- reproduce notebook-aligned numerical stages;
- regenerate dissertation results and figures;
- audit thesis-facing numerical claims;
- continue absorption/RAI calibration work;
- write the dissertation from a compact documentation surface.

This branch does not rewrite history or alter `main`.

## Retained Working Surface

- `src/` - maintained simulator package;
- `tests/` and `regression/baseline/` - validation and fixed baselines;
- original notebook and `notebook_sections/` - historical Version 1 reference;
- `configs/` - canonical notebook, plot, result, and numerical contracts;
- `scripts/` - validation, canonical recovery, audit, and result generation;
- `results/notebook_aligned_recovery/` - canonical staged recovery outputs;
- `results/dissertation_plots_v1/` and `results/faraday_optimisation_v1/`;
- `results/thesis_numerical_consistency_v1/` and supporting validity audits;
- thesis-facing physics, architecture, calibration, result, and publication
  documents indexed by `docs/README.md`.

The large PCI/DGI and Faraday `.npz` files remain because regression tests use
them directly. The original notebook remains because it records the historical
computational pipeline.

## Removed From The Active Surface

- archived milestone reports and migration-time status notes;
- superseded notebook-pipeline, unit, parameter, and figure-recovery reviews;
- orphaned deep-discussion result tables superseded by the full physics audit;
- the obsolete pre-canonical single condensate recovery workflow;
- the tracked generated code-bundle zip.

These deletions remove duplicate process history, not scientific capability.
The deleted records remain recoverable from Git history. A portable zip can be
regenerated with `python scripts/create_code_bundle.py`; generated zips are
ignored by Git on this branch.

## Scientific Status

The branch preserves the Version 1 notebook-aligned simulator and its tested
numerical behavior. Results remain representative and uncalibrated. In
particular, `kappa_F = 1.0` is a placeholder, no real absorption/RAI fit has
been applied, and no microscopic multi-level Faraday model is claimed.

## Required Validation

From Windows PowerShell:

```powershell
$env:PYTHONPATH = "src;."
$env:PYTHONUTF8 = "1"
pytest -q
python scripts\validate_notebook_sections.py
python scripts\run_all_dissertation_figures.py --dry-run
```

The branch should only be used after all three commands pass and the working
tree is clean.

Validated on 14 July 2026:

- `pytest -q`: `127 passed`;
- notebook section validation: passed;
- dissertation generation dry run: all canonical steps resolved successfully.
