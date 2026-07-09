# Faraday Imaging Baseline Milestone Report

## Objective

Generate deterministic Faraday imaging baseline arrays before introducing any
Faraday orchestration helper.

This milestone is baseline generation only. It does not redesign the Faraday
model, introduce a microscopic circular-transition model, migrate notebook
logic, or add optimisation.

## Arrays Saved

The baseline file is:

```text
regression/baseline/imaging/faraday_imaging_baseline_v1.npz
```

It contains:

- `grid_axis_m`
- `grid_a_m`
- `grid_b_m`
- `spatial_frequency_axis_m_inv`
- `spatial_frequency_x_m_inv`
- `spatial_frequency_y_m_inv`
- `radii_m`
- `column_density_m2`
- `scalar_phase_peak_rad`
- `theta_f_peak_rad`
- `thomas_fermi_profile`
- `scalar_phase_map_rad`
- `theta_f_map_rad`
- `pupil`
- `sigma_plus_object_field`
- `sigma_minus_object_field`
- `sigma_plus_scattered_field`
- `sigma_minus_scattered_field`
- `sigma_plus_propagated_scattered_field`
- `sigma_minus_propagated_scattered_field`
- `sigma_plus_field`
- `sigma_minus_field`
- `output_ex_field`
- `output_ey_field`
- `dark_field_intensity`
- `dual_port_u_intensity`
- `dual_port_v_intensity`
- `dual_port_signal`
- `metadata_json`

## Parameters Used

The baseline uses the current notebook section 17.2 reference Faraday setup:

- imaging axis: `x`
- transverse plane: `y,z`
- detuning: `1.5e9 Hz`
- grid size: `1024 x 1024`
- field of view: `100e-6 m`
- wavelength: `401e-9 m`
- probe diameter: `24e-3 m`
- first imaging lens focal length: `150e-3 m`
- numerical aperture: `(D_probe / 2) / f1`
- `kappa_F = 1.0`

Atomic and optical parameters match the current reference notebook values for
the `166Er` 401 nm transition and BEC operating point.

## Generation Method

The generator is:

```text
scripts/generate_faraday_imaging_baseline.py
```

It directly mirrors the current notebook-equivalent Faraday equations from
notebook section 17.2 / cell 51:

```text
theta_F_peak = kappa_F * phi_peak
sigma_plus_object_field = exp(+1j * theta_F_map)
sigma_minus_object_field = exp(-1j * theta_F_map)
Pp = 1 + ifft2(fft2(sigma_plus_object_field - 1) * pupil)
Pm = 1 + ifft2(fft2(sigma_minus_object_field - 1) * pupil)
Ex = (Pp + Pm) / 2
Ey = 1j * (Pp - Pm) / 2
I_dark = abs(Ey)**2
I_u = abs(Ex + Ey)**2 / 2
I_v = abs(Ex - Ey)**2 / 2
S = (I_v - I_u) / (I_v + I_u)
```

This is a representative deterministic notebook-equivalent array baseline, not
a full notebook execution dump. No stochastic camera noise is included.

## Regeneration

From the repository root:

```bash
python scripts/generate_faraday_imaging_baseline.py
```

This rewrites:

```text
regression/baseline/imaging/faraday_imaging_baseline_v1.npz
```

## Tests Added

The regression test is:

```text
tests/regression/test_faraday_imaging_baseline.py
```

It checks:

- the baseline file exists,
- required keys are present,
- array shapes are as expected,
- arrays contain finite values,
- intensity-like arrays are non-negative,
- `theta_f_map_rad` is finite,
- metadata is present and records `kappa_F = 1.0`,
- regenerating the baseline produces arrays equal to the stored `.npz` values.

## Validation Results

Validation was run with:

```powershell
$env:PYTHONPATH="src;."
$env:PYTHONUTF8="1"
pytest -q
python scripts\validate_notebook_sections.py
python scripts\create_code_bundle.py
```

Results:

```text
pytest -q: 24 passed
notebook section validation: passed
code bundle generation: passed and remained reproducible
```

Generated baseline hash:

```text
regression/baseline/imaging/faraday_imaging_baseline_v1.npz:
  fe7977423fc05461b247ea6dc841d3e91e7895a15a2e97ec810533a9afa888d7
```

## Scope Confirmation

`kappa_F` remains `1.0`.

No microscopic Faraday model was introduced.

No Clebsch-Gordan, vector-polarizability, or detuning-dependent `kappa_F`
calibration was introduced.

No physics equations were changed.

No helper APIs were changed.

No notebook logic was changed.

No simulator behaviour was changed.

No `simulate_faraday_image(...)` helper was introduced.

Camera, noise, shot-noise, SNR, optimisation, and multi-shot migration were not
touched.
