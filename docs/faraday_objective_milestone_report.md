# Faraday Objective Milestone Report

## Objective

Add a small deterministic helper for evaluating one Faraday imaging operating
point.

## Helper Added

`evaluate_faraday_operating_point(...)` was added to
`src/non_destructive_image/analysis.py` and exported from
`src/non_destructive_image/__init__.py`.

The helper combines existing notebook-equivalent light-atom helpers:

- `faraday_rotation_angle(...)`
- `scattered_photons_per_atom(...)`
- `reabsorption_fraction(...)`

It returns a dictionary containing:

- `faraday_signal_rad`
- `faraday_signal_scale`
- `scattered_photons_per_atom`
- `reabsorption_fraction`
- `destructiveness_metric`
- `estimated_per_frame_snr`
- `signal_per_scattered_photon`
- `information_per_scattered_photon`
- `signal_to_destruction`

## Scope

This is a single-operating-point deterministic objective helper.

It does not implement:

- a parameter scan;
- plotting;
- stochastic noise averaging;
- microscopic Faraday physics;
- calibration of `kappa_F`.

The current phenomenological Faraday convention remains unchanged:

```text
theta_F = kappa_F * phi_peak
kappa_F = 1.0
```

## Tests Added

`tests/regression/test_faraday_objective.py` checks:

- required output keys are present;
- outputs are finite;
- scattered photons are non-negative;
- destructiveness and signal-ratio metrics are non-negative;
- `signal_per_scattered_photon` matches the expected ratio;
- lower probe power improves signal per scattered photon in a controlled case;
- identical inputs produce identical outputs.

## Scope Confirmation

Validation results:

```text
pytest -q: 42 passed
python scripts\validate_notebook_sections.py: passed
```

No notebook sections were changed.

No existing imaging helpers were changed.

No camera helpers were changed.

No multi-shot helpers were changed.

No Atomic or Light-Atom helpers were changed.

No baselines were changed.

No deliverable zip was regenerated.

`docs/optimisation_readiness.md` was not changed.
