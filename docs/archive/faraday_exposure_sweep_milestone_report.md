# Faraday Exposure-Time Sweep Milestone Report

## Scope

This milestone adds a small deterministic one-dimensional exposure-time scaling
helper for Faraday operating-point evaluation. It uses the existing
`pulse_duration_s` input from `evaluate_faraday_operating_point(...)` as the
swept exposure-time control parameter.

The helper remains intentionally narrow. It does not add two-dimensional or
three-dimensional sweeps, plotting, stochastic noise averaging, microscopic
Faraday physics, or multi-parameter optimisation.

## Helper Added

- `sweep_faraday_exposure_time(...)`

The helper accepts a non-empty one-dimensional list or array of positive
`pulse_duration_s` values and returns:

- `pulse_duration_s`
- all scalar objective quantities returned by
  `evaluate_faraday_operating_point(...)`, stacked into arrays
- `objective_key`
- `best_index`
- `best_pulse_duration_s`
- `best_objective_value`

The default objective is `signal_per_scattered_photon`.

## Behaviour Preserved

The helper delegates each operating-point calculation to
`evaluate_faraday_operating_point(...)`. It does not alter:

- `evaluate_faraday_operating_point(...)` behaviour
- `sweep_faraday_detuning(...)` behaviour
- `sweep_faraday_intensity(...)` behaviour
- the phenomenological `kappa_f = 1.0` convention
- the Faraday physics model
- imaging, camera, or multi-shot helpers
- atomic or light-atom helpers
- notebook sections
- regression baselines
- deliverable bundle files

## Tests Added

Added `tests/regression/test_faraday_exposure_sweep.py`, covering:

- output length matching exposure-time input length
- finite output arrays
- best index matching the maximum selected objective
- deterministic repeated evaluation
- single exposure-time sweep behaviour
- clear error handling for empty input
- clear error handling for non-positive exposure-time values

## Validation

Validation commands:

```powershell
$env:PYTHONPATH="src;."
$env:PYTHONUTF8="1"
pytest -q
python scripts\validate_notebook_sections.py
```

Results:

- `pytest -q`: 61 passed
- `python scripts\validate_notebook_sections.py`: passed

## Confirmation

No notebook logic, physics equations, imaging helpers, camera helpers,
multi-shot helpers, atomic/light-atom helpers, baselines, or deliverable bundle
files were changed.
