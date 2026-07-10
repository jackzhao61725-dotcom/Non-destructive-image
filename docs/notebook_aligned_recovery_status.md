# Notebook-Aligned Recovery Status

## Current Recovery Path

The canonical notebook-aligned recovery path now covers:

```text
parameters -> Thomas-Fermi condensate -> projected profile -> column-density map -> scalar phase map -> PCI image -> DGI image -> ideal Faraday outputs
-> deterministic camera image -> stochastic camera frame -> deterministic multishot sequence -> noisy PCI multishot filmstrip
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

### Condensate Three-View Projection Extension

Script:

```text
scripts/generate_condensate_three_view.py
```

Status:

- notebook-aligned condensate-model extension;
- not an exact notebook figure recovery unless the historical notebook is later
  shown to contain the same three-view layout;
- not experimental calibration;
- not a final calibrated prediction;
- uses the same canonical condensate parameters, grid, field of view, and unit
  conventions as `scripts/recover_notebook_condensate_stage.py`;
- plots absolute Thomas-Fermi column-density distributions in `m^-2`, not
  normalised density;
- labels 2D maps as `n_col(y,z)`, `n_col(x,z)`, or `n_col(x,y)`;
- reserves `n_tilde_x`, `n_tilde_y`, and `n_tilde_z` for peak scalar values
  matching the thesis parameter table.

Axis conventions:

```text
integrate along x -> display y-z plane
integrate along y -> display x-z plane
integrate along z -> display x-y plane
```

Generated outputs:

```text
results/notebook_aligned_recovery/condensate_three_view/condensate_three_view.svg
results/notebook_aligned_recovery/condensate_three_view/condensate_three_view_summary.json
results/notebook_aligned_recovery/condensate_three_view/central_lineouts.csv
results/notebook_aligned_recovery/condensate_three_view/metadata.json
```

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

### Stochastic Camera

Script:

```text
scripts/recover_notebook_noisy_camera_stage.py
```

Recovered notebook references:

- cell 0: `rng = np.random.default_rng(7)`;
- cell 20: `to_camera(...)`.

Recovered stochastic convention:

```text
counts = rng.poisson(np.clip(binned, 0, None) * N_phot_pix)
       + rng.normal(0, read_e, binned.shape)
noisy_camera_image = counts / N_phot_pix
```

Status:

- closed for the tested explicit-seed stochastic recipe;
- input: recovered deterministic PCI camera image / binned PCI ideal image;
- RNG seed used for recovery: `7`;
- notebook RNG policy: global `np.random.default_rng(7)`;
- recovery RNG policy: explicit `np.random.default_rng(seed)`;
- exact replay is available for the isolated first PCI camera call under clean
  notebook execution order;
- arbitrary interactive notebook-global RNG state is not claimed to be
  reproducible;
- photons per camera pixel: `1532.3236613066642`;
- read noise: `7.0 e-`;
- noisy camera image shape: `68 x 68`;
- noisy camera image centre value: `1.1547265453246764`;
- noisy camera image mean: `0.9047881903189995`;
- deterministic camera image mean: `0.9043924044803204`;
- residual standard deviation: `0.02517625607875563`;
- residual standard deviation divided by mean expected per-pixel image noise:
  `1.0184981488308726`;
- all tested seeded stochastic helper/notebook-expression max absolute
  differences against `simulate_noisy_camera_image(...)`: `0.0`.

### Deterministic Multishot

Script:

```text
scripts/recover_notebook_multishot_stage.py
```

Recovered notebook references:

- cell 40: `run_sequence(...)` and `accumulate(...)`;
- cell 42: deterministic sequence-evolution display.

Recovered deterministic convention:

```text
Ng = N_scatt(Delta_Hz, P_mW, tau_s)
dE = Ng * (1 + reabs) * E_rec
A_E = 3 * (zeta4 / zeta3) * kB / Tc**3
T_next = (T**4 + dE / A_E)**0.25
N0_heating = N_tot_sc * (1 - (T / Tc)**3)
N0_clean_loss = N0_0 * exp(-eta_coll * Ng * s)
phi_s = phi_peak(Delta_Hz, ncol_axis(N0_s))
accumulated_snr = sqrt(nancumsum(where(isfinite(snr), snr**2, 0)))
```

Status:

- closed for tested deterministic sequence quantities;
- operating point: `Delta = 1.5 GHz`, `P = 3.5 mW`, `tau = 40 us`;
- imaging axis: `0` / `x`;
- stop condition: `30%` condensate loss;
- max-shot safety cap: `400`;
- `eta_coll = 1.3`;
- reabsorption fraction at this operating point: `0.029708652968257532`;
- scattered photons per atom per shot: `0.009274742967987243`;
- self-consistent total atom number: `116166.77870472142`;
- critical temperature: `216.8262079928517 nK`;
- heating sequence length: `15` frames, last shot `14`;
- clean-loss sequence length: `31` frames, last shot `30`;
- heating final condensate number: `17369.583984817014`;
- heating final loss fraction: `0.30521664060731946`;
- heating final accumulated SNR: `27.281738255565195`;
- clean-loss final condensate number: `17412.021337141854`;
- clean-loss final loss fraction: `0.3035191465143259`;
- clean-loss final accumulated SNR: `39.020132203959946`;
- all tested deterministic multishot helper/notebook-expression max absolute
  differences against `simulate_multishot_sequence(...)`: `0.0`.

### Noisy PCI Multishot Filmstrip

Script:

```text
scripts/recover_notebook_noisy_multishot_filmstrip.py
```

Recovered notebook references:

- cell 44: related start/middle/end camera-frame variant;
- cell 93: Step 14 noisy PCI filmstrip.

Recovered cell 93 convention:

```text
st_show = [0, 5, 10, 14]
st_Nd_run = N_phot_pix(SEQ_P_mW, SEQ_tau)
phi_s = seq_h['phi'][s]
I_s, _ = sim_image(SEQ_axis, phi_s, 'PCI')
b = I_s[:_nb2, :_nb2].reshape(_nb2//15, 15, _nb2//15, 15).mean(axis=(1, 3))
frame = (rng.poisson(clip(b, 0, None) * st_Nd_run)
         + rng.normal(0, read_e, b.shape)) / st_Nd_run
```

Status:

- closed for the tested explicit-seed noisy PCI filmstrip recipe;
- primary recovered notebook cell: `93`;
- related notebook cell 44 frame selection: `[0, 7, 14]`;
- recovered Step 14 selected frames: `[0, 5, 10, 14]`;
- profile convention: cell 93 reuses the fresh-cloud PCI profile and applies
  each frame's phase value;
- photons per camera pixel at the run operating point: `1072.6265629146646`;
- read noise: `7.0 e-`;
- noisy frame shape: `68 x 68`;
- frame 0 noisy mean: `0.9048341400329197`;
- frame 5 noisy mean: `0.9045503168872816`;
- frame 10 noisy mean: `0.9044529944801176`;
- frame 14 noisy mean: `0.9043769748075294`;
- seeded helper replay max absolute differences for binned images, noisy
  counts, and noisy frames: `0.0`.

RNG status:

- notebook RNG policy: global `np.random.default_rng(7)`;
- recovery RNG policy: explicit `np.random.default_rng(seed)`;
- exact full-notebook global RNG reproduction is not claimed, because cell 93
  appears after earlier stochastic camera and detuning-strip draws;
- the recovered filmstrip is exactly reproducible under the explicit-seed
  replay used by the recovery script.

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
- `camera_stage/`;
- `noisy_camera_stage/`;
- `multishot_stage/`;
- `noisy_multishot_filmstrip/`.

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

The noisy camera stage generates only:

- `noisy_camera_stage.svg`;
- JSON comparison/summary/metadata files;
- noise statistics CSV.

It does not generate multi-shot frame sequences.

The multishot stage generates only:

- `multishot_sequence_stage.svg`;
- JSON comparison/summary/metadata files;
- deterministic sequence CSV.

It does not generate noisy multishot filmstrips or operating maps.

The noisy multishot filmstrip stage generates only:

- `noisy_multishot_pci_filmstrip.svg`;
- JSON comparison/summary/metadata files;
- selected-frame and frame-statistics CSV files.

It does not generate Faraday camera panels, dual-port flicker, operating maps,
or shot-noise maps.

## Regression Tests

Focused regression tests:

```text
tests/regression/test_notebook_condensate_recovery.py
tests/regression/test_notebook_phase_recovery.py
tests/regression/test_notebook_pci_recovery.py
tests/regression/test_notebook_dgi_recovery.py
tests/regression/test_notebook_faraday_recovery.py
tests/regression/test_notebook_camera_recovery.py
tests/regression/test_notebook_noisy_camera_recovery.py
tests/regression/test_notebook_multishot_recovery.py
tests/regression/test_notebook_noisy_multishot_filmstrip.py
```

These tests check stable numerical outputs and helper/notebook-expression
agreement. They do not test SVG pixel appearance.

## What Remains Uncertain

The recovery still has not locked down end-to-end notebook-aligned recipes for:

- Faraday camera-level reference panels;
- Faraday dual-port noisy frame and flicker robustness;
- optimisation and SNR operating maps.

Display styling is documented only where it is needed to avoid plotting the
wrong physical quantity. The numerical arrays remain the primary recovery
target.

## Recommended Next Step

The next candidate recovery should be the Faraday camera-level reference panel:

```text
ideal Faraday dark-field and dual-port outputs -> camera-level Faraday reference panel
```

That recovery should keep the ideal Faraday outputs and camera/noise recipes
fixed. Dual-port flicker and operating maps should remain separate later
recovery targets.
