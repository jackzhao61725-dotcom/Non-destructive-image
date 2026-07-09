# Faraday Detuning Sweep Milestone Report

## Scope

This milestone adds a small deterministic detuning-sweep helper for evaluating a list of Faraday imaging operating points. The helper reuses `evaluate_faraday_operating_point(...)` for each detuning and collects the returned objective quantities into arrays.

The implementation is intentionally limited to a one-dimensional detuning sweep. It does not add intensity sweeps, exposure-time sweeps, plotting, stochastic noise averaging, microscopic Faraday physics, or multi-parameter optimisation.

## Helper Added

- `sweep_faraday_detuning(...)`

The helper accepts a non-empty one-dimensional list or array of detuning values and returns:

- `detuning_hz`
- all scalar objective quantities returned by `evaluate_faraday_operating_point(...)`, stacked into arrays
- `objective_key`
- `best_index`
- `best_detuning_hz`
- `best_objective_value`

The default objective is `signal_per_scattered_photon`.

## Behaviour Preserved

The sweep helper delegates the operating-point calculation to the existing deterministic objective helper. It does not alter:

- `evaluate_faraday_operating_point(...)` behaviour
- the phenomenological `kappa_f = 1.0` convention
- light-atom helper formulas
- imaging, camera, or multi-shot helper behaviour
- notebook sections
- regression baselines

## Tests Added

Added `tests/regression/test_faraday_detuning_sweep.py`, covering:

- output length matching detuning input length
- finite output arrays
- best index matching the maximum selected objective
- deterministic repeated evaluation
- single-detuning sweep behaviour
- clear error handling for an empty sweep

## Validation

Validation commands:

```powershell
$env:PYTHONPATH="src;."
$env:PYTHONUTF8="1"
pytest -q
python scripts\validate_notebook_sections.py
```

Results:

- `pytest -q`: 47 passed
- `python scripts\validate_notebook_sections.py`: passed

## Confirmation

No notebook logic, physics equations, imaging helpers, camera helpers, multi-shot helpers, atomic/light-atom helpers, baselines, or deliverable bundle files were changed.
