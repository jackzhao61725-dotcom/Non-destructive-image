# Multi-Shot Sequence Core Milestone Report

## Objective

Migrate only the deterministic multi-shot sequence core from the notebook into
a small helper module.

## Helpers Added

`src/non_destructive_image/multishot.py` adds:

- `simulate_multishot_sequence(...)`
- `accumulate_snr(...)`

The helpers are exported from `src/non_destructive_image/__init__.py`.

## Migrated Logic

`simulate_multishot_sequence(...)` preserves the notebook sequence bookkeeping:

```text
clean-loss model:
N0_now = N0_0 * exp(-eta_coll * N_gamma * shot)

heating model:
N0_now = N_tot_sc * (1 - (T / Tc)^3)
T_next = (T^4 + dE / A_E)^0.25

stop condition:
stop when loss_fraction >= loss_fraction_limit
```

The helper returns per-frame arrays for:

- `shot`
- `N0`
- `condensate_fraction`
- `loss_fraction`
- `frac`
- `T`
- `phi`
- `snr`

`accumulate_snr(...)` preserves the notebook RMS accumulated-SNR convention:

```text
sqrt(nancumsum(where(isfinite(snr), snr**2, 0)))
```

## Callback Boundary

The helper accepts callbacks for phase and deterministic SNR:

```text
phase_from_n0(...)
snr_from_phi(...)
```

This keeps imaging, camera noise, noisy frame rendering, detuning sweeps, and
plotting outside this migration.

## Tests Added

`tests/regression/test_multishot_sequence_core.py` checks:

- accumulated SNR matches the notebook RMS convention;
- clean-loss sequence values match hand-computed expectations;
- heating sequence values match hand-computed expectations;
- output arrays have consistent lengths;
- repeated calls with the same inputs are deterministic.

## Scope Confirmation

Validation results:

```text
pytest -q: 37 passed
python scripts\validate_notebook_sections.py: passed
```

No notebook sections were changed.

No imaging helpers were changed.

No camera helpers were changed.

No Atomic or Light-Atom helpers were changed.

No baseline `.npz` files were changed.

No noisy camera or frame-rendering code was migrated.

No Faraday dual-port frame sequence code was changed.

No detuning sweep, operating-map, plotting, SNR-map, or optimisation code was
changed.

No deliverable zip was regenerated.
