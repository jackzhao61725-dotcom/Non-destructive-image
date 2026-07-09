# Faraday Orchestration Milestone Report

## Objective

Add a thin `simulate_faraday_image(...)` helper that reproduces the existing
Faraday imaging baseline without changing notebook physics or lower-layer
helpers.

## Migrated Orchestration

The helper accepts an existing `theta_F` map and pupil. It preserves the current
notebook-equivalent Faraday conventions:

```text
sigma_plus_object_field = exp(+1j * theta_f_map)
sigma_minus_object_field = exp(-1j * theta_f_map)
Pp = 1 + ifft2(fft2(sigma_plus_object_field - 1) * pupil)
Pm = 1 + ifft2(fft2(sigma_minus_object_field - 1) * pupil)
Ex = (Pp + Pm) / 2
Ey = 1j * (Pp - Pm) / 2
I_dark = abs(Ey)**2
I_u = abs(Ex + Ey)**2 / 2
I_v = abs(Ex - Ey)**2 / 2
S = (I_v - I_u) / (I_v + I_u)
```

The phenomenological convention `theta_F = kappa_F * phi_peak` with
`kappa_F = 1.0` remains unchanged. The helper does not introduce a microscopic
Faraday model.

## Baseline Used

Regression testing compares helper outputs against:

```text
regression/baseline/imaging/faraday_imaging_baseline_v1.npz
```

The focused test checks:

- `theta_f_map_rad`
- `sigma_plus_object_field`
- `sigma_minus_object_field`
- `sigma_plus_field`
- `sigma_minus_field`
- `output_ex_field`
- `output_ey_field`
- `dark_field_intensity`
- `dual_port_u_intensity`
- `dual_port_v_intensity`
- `dual_port_signal`

## Scope Confirmation

Validation results:

```text
pytest -q: 25 passed
python scripts\validate_notebook_sections.py: passed
```

No notebook sections were changed.

No baseline generator scripts were changed.

No existing `.npz` baselines were changed.

No Atomic or Light-Atom helper APIs were changed.

No PCI, DGI, camera, noise, or multi-shot behaviour was changed.
