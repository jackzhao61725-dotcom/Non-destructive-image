# Faraday Intensity Sweep Milestone Report

## Scope

This milestone adds a small deterministic one-dimensional probe-intensity
scaling helper for Faraday operating-point evaluation. It uses the existing
`probe_power_mw` input from `evaluate_faraday_operating_point(...)` as the
swept intensity-like control parameter.

The helper remains intentionally narrow. It does not add a two-dimensional
detuning-intensity sweep, exposure-time sweep, plotting, stochastic noise
averaging, microscopic Faraday physics, or multi-parameter optimisation.

## Helper Added

- `sweep_faraday_intensity(...)`

The helper accepts a non-empty one-dimensional list or array of positive
`probe_power_mw` values and returns:

- `probe_power_mw`
- all scalar objective quantities returned by
  `evaluate_faraday_operating_point(...)`, stacked into arrays
- `objective_key`
- `best_index`
- `best_probe_power_mw`
- `best_objective_value`

The default objective is `signal_per_scattered_photon`.

## Behaviour Preserved

The helper delegates each operating-point calculation to
`evaluate_faraday_operating_point(...)`. It does not alter:

- `evaluate_faraday_operating_point(...)` behaviour
- `sweep_faraday_detuning(...)` behaviour
- the phenomenological `kappa_f = 1.0` convention
- the Faraday physics model
- imaging, camera, or multi-shot helpers
- atomic or light-atom helpers
- notebook sections
- regression baselines
- deliverable bundle files

## Tests Added

Added `tests/regression/test_faraday_intensity_sweep.py`, covering:

- output length matching probe-power input length
- finite output arrays
- best index matching the maximum selected objective
- deterministic repeated evaluation
- single probe-power sweep behaviour
- clear error handling for empty input
- clear error handling for non-positive probe-power values

## Validation

Validation commands:

```powershell
$env:PYTHONPATH="src;."
$env:PYTHONUTF8="1"
pytest -q
python scripts\validate_notebook_sections.py
```

Results:

- `pytest -q`: 54 passed
- `python scripts\validate_notebook_sections.py`: passed

## Confirmation

No notebook logic, physics equations, imaging helpers, camera helpers,
multi-shot helpers, atomic/light-atom helpers, baselines, or deliverable bundle
files were changed.
