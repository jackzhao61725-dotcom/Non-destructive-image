# Notebook-Aligned Recovery Status

## Current Recovery Path

The canonical notebook-aligned recovery path now covers:

```text
parameters -> Thomas-Fermi condensate -> projected profile -> column-density map -> scalar phase map -> PCI image -> DGI image -> ideal Faraday outputs
-> deterministic camera image
```

This is a recovery of the historical Version 1 notebook computation. It does
not validate the model against experiment and does not introduce calibrated
physics.

## Closed Stages

### Condensate

Script:

```text
scripts/recover_notebook_condensate_stage.py
```

Status:

- closed for tested deterministic quantities;
- grid shape: `1024 x 1024`;
- peak location: centre pixel `[512, 512]`;
- helper/notebook-expression differences for radii, column density, profile,
  and projected column-density map: `0.0`.

### Scalar Phase

Script:

```text
scripts/recover_notebook_phase_stage.py
```

Recovered convention:

```text
delta = 2 * Delta_Hz * 2*pi / Gamma
phi_peak = sigma0 * n_col_peak * delta / (2 * (1 + delta**2))
phase_map = phi_peak * projected_profile
```

Status:

- closed for tested deterministic quantities;
- detuning: `1.5e9 Hz`;
- dimensionless detuning: `101.69491525423727`;
- peak scalar phase: `0.20294165287929014 rad`;
- phase-map helper/notebook-expression max absolute difference: `0.0`.

### PCI

Script:

```text
scripts/recover_notebook_pci_stage.py
```

Recovered convention:

```text
object_field = np.exp(1j * phase_map)
scattered_field = object_field - 1
propagated_scattered_field = np.fft.ifft2(np.fft.fft2(scattered_field) * pupil)
reference_field = t_p * np.exp(1j * theta)
pci_image_intensity = np.abs(reference_field + propagated_scattered_field) ** 2
```

Status:

- closed for tested deterministic quantities;
- `t_p = 0.95`;
- `theta = pi/2`;
- numerical aperture: `0.08`;
- pupil nonzero pixels: `1245`;
- PCI image centre value: `1.1585677695076864`;
- PCI image max absolute difference against `simulate_pci_image(...)`: `0.0`.

### DGI

Script:

```text
scripts/recover_notebook_dgi_stage.py
```

Recovered convention:

```text
object_field = np.exp(1j * phase_map)
scattered_field = object_field - 1
propagated_scattered_field = np.fft.ifft2(np.fft.fft2(scattered_field) * pupil)
reference_field = 10 ** (-OD / 2)
dgi_image_intensity = np.abs(reference_field + propagated_scattered_field) ** 2
```

Status:

- closed for tested deterministic quantities;
- `OD = 4.0`;
- residual reference field: `0.01`;
- residual reference intensity: `0.0001`;
- DGI image centre value: `0.01595651906354294`;
- DGI image max absolute difference against `simulate_dgi_image(...)`: `0.0`.

### Faraday

Script:

```text
scripts/recover_notebook_faraday_stage.py
```

Recovered notebook references:

- cell 39: `kappa_F = 1.0` placeholder;
- cell 43: `theta_F_peak(...)`;
- cell 51: `sim_faraday_fields(...)` and `faraday_maps(...)`;
- cell 85: ideal dark-field Faraday display logic;
- cell 89: ideal dual-port Faraday display logic.

Recovered convention:

```text
theta_f_map = kappa_F * scalar_phase_map
sigma_plus_object_field = np.exp(+1j * theta_f_map)
sigma_minus_object_field = np.exp(-1j * theta_f_map)
sigma_plus_field = 1 + np.fft.ifft2(np.fft.fft2(sigma_plus_object_field - 1) * pupil)
sigma_minus_field = 1 + np.fft.ifft2(np.fft.fft2(sigma_minus_object_field - 1) * pupil)
Ex = (sigma_plus_field + sigma_minus_field) / 2
Ey = 1j * (sigma_plus_field - sigma_minus_field) / 2
I_dark = np.abs(Ey) ** 2
I_u = np.abs(Ex + Ey) ** 2 / 2
I_v = np.abs(Ex - Ey) ** 2 / 2
S = (I_v - I_u) / (I_v + I_u)
```

Status:

- closed for tested deterministic quantities;
- `kappa_F = 1.0` remains a Version 1 phenomenological placeholder;
- no microscopic Faraday model was introduced;
- no experimental calibration was applied;
- peak Faraday rotation: `0.20294165287929014 rad`;
- dark-field centre value: `0.015956443479481053`;
- dual-port signal centre value: `0.25116900542198173`;
- all tested Faraday helper/notebook-expression max absolute differences
  against `simulate_faraday_image(...)`: `0.0`.

Existing Faraday baseline compatibility:

- `regression/baseline/imaging/faraday_imaging_baseline_v1.npz` was read but
  not modified;
- it records `kappa_F = 1.0` and no microscopic Faraday model;
- canonical recovery differs from that older baseline at small nonzero levels
  for some arrays, for example `theta_f_map_rad` max absolute difference
  `< 3e-5` and dark-field intensity max absolute difference `< 5e-6`;
- the canonical notebook recovery and current helper output are the primary
  closure criteria for this stage.

### Deterministic Camera

Script:

```text
scripts/recover_notebook_camera_stage.py
```

Recovered notebook reference:

- cell 20: `N_phot_pix(...)` and `to_camera(...)`.

Recovered deterministic convention:

```text
Mag = f2 / f1
pix_obj = pix_cam / Mag
N_phot_pix = intensity_at_atoms(P_mW) * pix_obj**2 * tau_s * QE / E_phot
nb = (Ngrid // 15) * 15
binned = Iratio[:nb, :nb].reshape(nb//15, 15, nb//15, 15).mean(axis=(1, 3))
deterministic_counts = binned * N_phot_pix
camera_image = deterministic_counts / N_phot_pix
```

Status:

- closed for tested deterministic quantities;
- input ideal image: recovered PCI intensity image;
- probe power: `2.0 mW`;
- exposure used by the notebook camera helper default: `100 us`;
- quantum efficiency: `0.40`;
- read noise: `7.0 e-`, recorded but not applied in deterministic recovery;
- binning: `15 x 15`;
- trimmed high-resolution shape: `1020 x 1020`;
- camera image shape: `68 x 68`;
- detected photons per camera pixel: `1532.3236613066642`;
- deterministic camera centre value: `1.1339167724448815`;
- deterministic count centre value: `1737.5275003697766`;
- all tested deterministic camera helper/notebook-expression max absolute
  differences against `simulate_camera_image(...)`: `0.0`.

## Generated Recovery Outputs

Current recovery outputs are grouped under:

```text
results/notebook_aligned_recovery/
```

Stage directories:

- `condensate_stage/`;
- `phase_stage/`;
- `pci_stage/`;
- `dgi_stage/`;
- `faraday_stage/`;
- `camera_stage/`.

The Faraday stage generates only:

- `faraday_dark_field_stage.svg`;
- `faraday_dual_port_signal_stage.svg`;
- JSON comparison/summary/metadata files;
- central lineouts CSV.

It does not generate camera, noisy-frame, or multi-shot figures.

The camera stage generates only:

- `camera_deterministic_stage.svg`;
- JSON comparison/summary/metadata files;
- central lineouts CSV.

It does not generate Poisson photon noise, Gaussian read-noise frames, or
multi-shot sequences.

## Regression Tests

Focused regression tests:

```text
tests/regression/test_notebook_condensate_recovery.py
tests/regression/test_notebook_phase_recovery.py
tests/regression/test_notebook_pci_recovery.py
tests/regression/test_notebook_dgi_recovery.py
tests/regression/test_notebook_faraday_recovery.py
tests/regression/test_notebook_camera_recovery.py
```

These tests check stable numerical outputs and helper/notebook-expression
agreement. They do not test SVG pixel appearance.

## What Remains Uncertain

The recovery still has not locked down end-to-end notebook-aligned recipes for:

- noisy frame rendering;
- multi-shot frame rendering;
- optimisation and SNR operating maps.

Display styling is documented only where it is needed to avoid plotting the
wrong physical quantity. The numerical arrays remain the primary recovery
target.

## Recommended Next Step

The next candidate recovery should be the stochastic camera/noise stage:

```text
deterministic camera image -> noisy frame
```

That recovery should keep the existing deterministic camera stage fixed and
compare the notebook Poisson/read-noise recipe before any multi-shot frame
sequence is attempted.
