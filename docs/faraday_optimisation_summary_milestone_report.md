# Faraday Optimisation Summary Milestone Report

## Scope

This milestone adds a lightweight deterministic reporting helper for summarising
Faraday single-variable sweep outputs. It reads the dictionary returned by an
existing sweep helper and extracts the selected metric, best point, evaluated
point count, and small min/max diagnostics.

The helper does not write files, plot results, depend on pandas, run stochastic
noise averaging, or perform new optimisation searches.

## Helper Added

- `summarise_faraday_sweep(...)`

The helper accepts:

- the result dictionary from a deterministic Faraday sweep helper;
- a metric key, defaulting to `signal_per_scattered_photon`;
- an optional parameter key.

It returns a plain summary dictionary containing:

- `metric_key`
- `parameter_key`
- `best_index`
- `best_parameter_value`
- `best_metric_value`
- `num_points`
- `metric_min`
- `metric_max`
- min/max diagnostics for selected quantities already present in the sweep

## Behaviour Preserved

The helper only summarises existing sweep output. It does not alter:

- `evaluate_faraday_operating_point(...)` behaviour
- detuning sweep behaviour
- intensity sweep behaviour
- exposure-time sweep behaviour
- the phenomenological `kappa_f = 1.0` convention
- the Faraday physics model
- imaging, camera, or multi-shot helpers
- atomic or light-atom helpers
- notebook sections
- regression baselines
- deliverable bundle files

## Tests Added

Added `tests/regression/test_faraday_optimisation_summary.py`, covering:

- best index matching the maximum selected metric
- single-point sweep handling
- clear error handling for a missing metric
- clear error handling for an empty sweep
- deterministic repeated summaries

## Validation

Validation commands:

```powershell
$env:PYTHONPATH="src;."
$env:PYTHONUTF8="1"
pytest -q
python scripts\validate_notebook_sections.py
```

Results:

- `pytest -q`: 66 passed
- `python scripts\validate_notebook_sections.py`: passed

## Confirmation

No notebook logic, physics equations, imaging helpers, camera helpers,
multi-shot helpers, atomic/light-atom helpers, baselines, or deliverable bundle
files were changed.
