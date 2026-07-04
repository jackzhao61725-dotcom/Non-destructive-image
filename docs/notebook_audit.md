# Notebook Audit Report

Notebook audited: `1 calculations revised 2  multishot  6  extended.ipynb`

This audit is Step 1 of the refactoring workflow. No physics, numerical algorithms, parameter values, or notebook code have been changed. The purpose is to document the current structure, identify duplicated calculations and helper functions, and propose a safe refactoring path for Version 1.0.

A machine-readable companion inventory is available at `docs/notebook_audit.json`; it records the same cell/function audit in JSON for future refactor scripts and regression-check tooling.

## 1. Current notebook structure

The notebook currently has 96 cells: 51 markdown cells and 45 code cells. It is scientifically complete, but the section numbering and responsibilities have grown organically. The main sections are:

| Current cells | Current section | Main responsibility |
|---|---|---|
| 0-2 | Title, imports, constants | Imports, plotting defaults, random number generator, fundamental constants |
| 3-8 | Apparatus, condensate parameters, atomic constants, Thomas-Fermi profile | Defines Er, optical arm, BEC parameters, derived constants, TF radii, peak density, column densities |
| 9-14 | Dispersive phase and scattering | Detuning conversion, phase shift, residual OD, scattering rate, loss and heating bounds |
| 15-22 | Phase-to-camera image formation and shot-noise SNR | PCI/DGI Fourier optics, camera binning/noise, lineouts, SNR helper functions |
| 23-30 | Shot budget, spatial resolution, operating point summary | Threshold sweeps, resolution estimates, summary table/report |
| 31-38 | Pulse duration and operating-point optimisation | Realistic destruction models, pulse-duration sweep, operating maps, recommendation |
| 39-46 | Multi-shot run evolution and detuning sweep | Condensate evolution, heating/loss sequences, fading camera frames, detuning trade-off maps |
| 47-55 | Faraday imaging | Faraday rotation model, dark-field and dual-port Faraday simulations, Faraday SNR, sequence comparison |
| 56-63 | End-to-end signal explanation | Pedagogical five-stage explanation of shared imaging pipeline |
| 64-93 | Step-by-step imaging illustrations | Detailed visual walkthrough for phase map, PCI, DGI, Faraday, camera noise, and run filmstrip |
| 94-95 | Closing separators | Markdown separators only |

## 2. Helper functions currently present

### Atomic model / condensate state

- `_frac_residual(Ntot)` solves for the self-consistent total atom number that gives the requested condensate number.
- `_tf_profile(Ra, Rb)` builds a 2D Thomas-Fermi column profile on the imaging grid.
- `tf_state(Nc)` recalculates Thomas-Fermi chemical potential, peak density, radii, and column density for a changed condensate number.

### Light-atom interaction

- `delta_of(Delta_Hz)` converts Hz detuning to the dimensionless linewidth-normalised detuning.
- `phi_peak(Delta_Hz, n_col_peak)` computes the scalar dispersive phase shift.
- `od_resonant_equiv(Delta_Hz, n_col_peak)` computes the residual optical depth.
- `intensity_at_atoms(P_mW)` converts probe power to intensity at the cloud.
- `N_scatt(Delta_Hz, P_mW, tau_s=None)` computes scattered photons per atom per shot.
- `theta_F_peak(Delta_Hz, n_col_peak)` computes the Faraday rotation angle using the same dispersive lineshape and the placeholder `kappa_F` factor.

### Imaging

- `sim_image(axis, phi_peak_val, mode='PCI', OD=4.0)` simulates PCI, DGI, and clear-field propagation with the current static cloud.
- `blur_for_axis(axis, phi_test=0.1)` estimates NA blur for each imaging axis.
- `sim_image_state(axis, phi_val, R_state, mode='PCI', OD=4.0)` repeats the same Fourier-optics model as `sim_image`, but accepts updated Thomas-Fermi radii for multi-shot evolution.
- `sim_faraday_fields(axis, theta_F_val)` propagates the two circular-polarisation components and recombines them into linear components.
- `faraday_maps(Delta_Hz, axis=0)` returns dark-field and dual-port Faraday intensity maps.

### Camera and noise

- `N_phot_pix(P_mW, tau_s=None, QE=None)` computes detected photons per camera pixel.
- `to_camera(Iratio, P_mW, QE=None)` bins the simulated image to camera pixels and adds Poisson photon noise plus Gaussian read noise.
- `camera_from_image(Iratio, P_mW, tau_s)` duplicates `to_camera` with explicit pulse duration.
- `_bin_center(I)` bins an image and extracts the central camera pixel.

### SNR / destruction / analysis

- `Nmax_loss`, `Nmax_heat`, `Nmax_cleanloss`, and `Nmax_heating` compute shot budgets under different destruction assumptions.
- `reabs_frac(Delta_Hz)` estimates reabsorption probability.
- `I_full(phi)`, `regime_of(phi)`, `SNR_shot_ideal`, `SNR_pixel`, `SNR_reselem_sim`, `SNR_pixel_phi`, `SNRres_fast`, `snr_abs_pix`, `snr_pci_pix`, `usable_snr`, and `SNR_faraday_sim` implement related SNR calculations.
- `_pci_block(Delta_Hz, axis=0, half=1)` caches the signal block used for fast PCI SNR calculations.
- `usable_curve(D, n_images)`, `make_strip(Delta_GHz, n_frames)`, `run_sequence(...)`, and `accumulate(seq)` implement multi-shot analysis and visualisation support.
- `operating_point_report(...)` prints the selected operating-point summary.

## 3. Repeated code and duplicated calculations

### 3.1 Camera binning and noise

The 15-by-15 grid-cell camera binning appears in several places:

- `to_camera(...)`
- `camera_from_image(...)`
- `_bin_center(...)`
- `SNR_reselem_sim(...)`
- `_pci_block(...)`
- `SNR_faraday_sim(...)`

Recommended consolidation: introduce a shared `bin_to_camera_pixels(image)` helper and keep the random-noise addition in one `add_camera_noise(...)` / `to_camera(...)` function. This should be a mechanical extraction only.

### 3.2 Thomas-Fermi profile generation

The 2D profile expression

```python
np.maximum(0, 1 - x**2/Ra**2 - y**2/Rb**2)**1.5
```

appears in `_tf_profile(...)`, `sim_image_state(...)`, and later explanatory cells. Recommended consolidation: keep `_tf_profile` as the single static-grid implementation and add an explicit radii-aware variant, e.g. `thomas_fermi_column_profile(radii, axis)`.

### 3.3 Fourier-optics propagation

The pattern

```python
np.fft.ifft2(np.fft.fft2(object_field - reference_field) * pupil)
```

appears in PCI/DGI simulation, state-dependent PCI frames, and Faraday circular-polarisation propagation. Recommended consolidation: extract a single `propagate_scattered_field(object_field)` helper that preserves the exact FFT convention and pupil multiplication.

### 3.4 Photon-number calculations

`N_phot_pix(...)` is the central photon-number helper, but some code relies on its default pulse duration while other code passes `tau_s` explicitly. The current default uses a local `tau_set = 100e-6` when `tau_s is None`, which differs from the global `tau = 40e-6`. This may be intentional legacy behaviour, so it should not be changed during the audit. It should be documented before any later refactor.

### 3.5 Peak phase / column-density recalculation

The combination `phi_peak(Delta_Hz, n_col[axis])` is repeated throughout the notebook. Multi-shot code repeats the pattern with `tf_state(N0_now)` followed by `phi_peak(...)`. Recommended consolidation: add a small helper such as `peak_phase_for_axis(detuning_hz, column_density_vector, axis)` only if it improves readability without changing call semantics.

### 3.6 SNR variants

SNR calculations are spread across ideal PCI, realistic PCI per pixel, realistic PCI per resolution element, absorption, DGI, dark-field Faraday, and dual-port Faraday. These should remain separate scientifically, but they should be grouped under an Analysis/SNR section with docstrings explaining inputs, outputs, and validity regimes.

### 3.7 Plot setup and figure styling

Many plotting cells repeat axis limits, extent construction, centre-line extraction, colourbar setup, labels, and lineout plotting. Recommended consolidation: only extract non-physics plotting conveniences, such as `set_cloud_axes(ax)`, `camera_extent_um()`, and `central_lineout(image)`. Avoid refactoring plot data generation until image identity has been regression-checked.

## 4. Plotting sections

The notebook contains many plotting and explanatory visualisation sections. Current major plotting groups are:

- Table 1 phase shift vs detuning.
- PCI/DGI ideal phase-to-intensity curves.
- PCI/DGI camera frames and lineouts.
- Shot-budget vs detuning threshold family.
- Pulse-duration SNR / shot-budget sweep.
- Operating map over power and pulse duration.
- Multi-shot condensate evolution, accumulated SNR, and fading PCI camera frames.
- Detuning-sweep trade-off maps and frame strips.
- Faraday rotation table, dark-field and dual-port Faraday frames, Faraday SNR comparison, and Faraday sequence comparison.
- Pedagogical step-by-step pipeline figures for shared phase formation, PCI, DGI, dark-field Faraday, dual-port Faraday, and filmstrip evolution.

These plotting cells are valuable for the MSc narrative but currently mix explanatory plotting with reusable simulation functions. In Version 1.0, plotting should move after the reusable function definitions inside the relevant logical section, not into a separate package-style plotting API unless needed.

## 5. Proposed notebook structure

The requested logical structure maps cleanly to the current content as follows:

```text
00 Imports
   - imports, plotting defaults, random seed, constants

01 Parameters
   - Er transition constants
   - imaging-arm parameters
   - BEC operating point
   - detection and loss-model controls

02 Atomic Model
   - Thomas-Fermi initial-state calculations
   - self-consistent thermal pedestal calculation
   - build_condensate / tf_state-style helper
   - column_density-style helper
   - update_after_scattering-style helper or multi-shot state update note

03 Light-Atom Interaction
   - detuning conversion
   - scalar phase shift
   - residual OD
   - scattering rate
   - reabsorption and destruction budget helpers
   - Faraday rotation angle helper

04 PCI
   - shared Fourier propagation helper
   - PCI image simulation
   - PCI SNR helpers directly tied to the image model
   - PCI camera examples

05 DGI
   - DGI image simulation through the same propagation helper
   - DGI read-noise floor / camera examples

06 Faraday
   - Faraday polarisation model
   - dark-field and dual-port maps
   - Faraday SNR comparison
   - Faraday multi-shot readout comparison

07 Camera
   - photon-per-pixel helper
   - binning helper
   - shot noise and read noise helper
   - digitisation/camera normalisation helper

08 Shot Noise
   - ideal and realistic shot-noise SNR functions
   - regime labels and per-resolution-element SNR
   - accumulated-SNR invariance demonstration

09 Multi-shot Simulation
   - heating/loss models
   - run_sequence
   - accumulate
   - state-dependent images

10 Analysis
   - operating-point summary
   - pulse-duration sweep
   - operating map
   - detuning sweep
   - spatial resolution
   - pedagogical pipeline figures and narrative plots
```

Note: section order 07/08 can be debated because shot noise is part of the camera model. To match the requested structure exactly, camera binning/digitisation should be in 07 and SNR/noise analysis should be in 08.

## 6. Refactoring plan

### Step 1 — Audit only

- Create this audit report.
- Do not edit notebook code or physics.
- Confirm that the current notebook remains the reference implementation.

### Step 2 — Rename variables conservatively

- Rename only local variables where intent is unambiguous.
- Prefer names such as `phase_map`, `field_object`, `field_fourier`, `field_image`, `camera_image`, `column_density_peak`, and `detuning_hz`.
- Avoid renaming globally reused scientific parameters (`N0`, `R`, `Gamma`, `tau`, `n_col`, etc.) until there is a regression check, because they appear throughout the notebook.
- Execute the notebook after each small rename batch.

### Step 3 — Extract duplicated helpers

Mechanical extractions only:

1. `bin_to_camera_pixels(image, bin_size=15)`
2. `add_camera_noise(binned_image, photons_per_pixel, rng, read_noise_e)`
3. `camera_extent_um(fov)`
4. `thomas_fermi_profile_2d(radius_a, radius_b)` or `column_profile_for_axis(axis, radii)`
5. `propagate_scattered_field(scattered_field, pupil)`

Each extraction should initially be used in one existing call site, then rolled out to duplicates after verifying outputs.

### Step 4 — Reorganise notebook sections

- Move cells into the requested 00-10 structure without changing their internal calculations.
- Keep definitions before first use.
- Place demonstrations and plots immediately after the functions they exercise unless a later Analysis section is clearer.
- Preserve random seed placement so noisy figures remain reproducible.

### Step 5 — Improve markdown and documentation

- Add a short markdown block at each major section covering purpose, physics, inputs, and outputs.
- Add docstrings to every reusable function.
- Mark placeholder physics explicitly, especially `kappa_F`, without changing its value or behaviour.
- Document legacy conventions such as the `N_phot_pix` default pulse duration before changing anything.

### Step 6 — Regression checks for future refactor PRs

Before and after each refactor step, compare:

- Key printed scalar values: `sigma0`, `Isat`, `T_rec`, `mu/kB`, `R`, `N_check`, `phi_peak(1.5e9, n_col[0])`, and `N_scatt(1.5e9, 2.0)`.
- Representative image arrays from `sim_image` for PCI and DGI at the reference detuning.
- Representative Faraday arrays from `faraday_maps(1.5e9, axis=0)`.
- Multi-shot sequence arrays from `run_sequence` for the current reference settings.
- Saved figures where deterministic; for noisy plots, preserve the random seed and compare summary statistics if exact image comparison is fragile.

## 7. Risks and cautions

- The notebook is the reference implementation; refactoring must be behaviour-preserving.
- `N_phot_pix` has a surprising default pulse duration (`100e-6`) rather than the global `tau`; document and preserve it until explicitly reviewed.
- The same 15-cell binning convention is embedded in many functions; extraction must preserve array shape and central-pixel indexing.
- Faraday currently uses `kappa_F = 1.0` as a placeholder/bound. Do not replace it with atomic-structure calculations in Version 1.0.
- Some variables are short but standard physics notation. Rename only where it improves maintainability without obscuring the original equations.

## 8. Audit method and exact notebook inventory

The audit was generated by reading the notebook JSON directly and extracting cell type, first non-empty source line, top-level function definitions, and approximate plotting calls. The notebook itself was not edited.

### 8.1 Exact cell inventory

| Cell | Type | Heading / first line | Functions defined | Plot calls |
|---:|---|---|---|---:|
| 0 | markdown | `# Phase-contrast and dark-ground imaging on the Oxford Er apparatus — design calculations` | — | 0 |
| 1 | markdown | `## 1. Imports and fundamental constants` | — | 0 |
| 2 | code | `import numpy as np` | — | 1 |
| 3 | markdown | `## 2. Apparatus and condensate parameters` | — | 0 |
| 4 | code | `# ---- 166Er and the 401 nm transition (K24 3.2.2) ----` | — | 0 |
| 5 | markdown | `## 3. Atomic and optical constants` | — | 0 |
| 6 | code | `sigma0 = 3*lam**2 / (2*np.pi)` | — | 0 |
| 7 | markdown | `## 4. Condensate profile (Thomas-Fermi)` | — | 0 |
| 8 | code | `omega = 2*np.pi*trap_Hz` | `_frac_residual` | 0 |
| 9 | markdown | `## 5. Dispersive phase shift` | — | 0 |
| 10 | code | `def delta_of(Delta_Hz):` | `delta_of`, `phi_peak`, `od_resonant_equiv` | 0 |
| 11 | markdown | `### Table 1 — Peak phase shift at representative detunings` | — | 0 |
| 12 | code | `detunings_GHz = [0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0]` | — | 3 |
| 13 | markdown | `## 6. Photon scattering and per-shot loss` | — | 0 |
| 14 | code | `def intensity_at_atoms(P_mW):` | `intensity_at_atoms`, `N_scatt`, `Nmax_loss`, `Nmax_heat` | 0 |
| 15 | markdown | `## 7. Image formation — phase to camera intensity` | — | 0 |
| 16 | code | `phi = np.linspace(0, 2.6, 600)` | — | 7 |
| 17 | markdown | `### 7.3 Fourier-optics simulation through the real arm` | — | 0 |
| 18 | code | `Ngrid, FOV = 1024, 100e-6` | `_tf_profile`, `sim_image`, `blur_for_axis` | 0 |
| 19 | markdown | `### 7.4 Simulated camera frames and lineouts` | — | 0 |
| 20 | code | `Mag = f2/f1; pix_obj = pix_cam/Mag` | `N_phot_pix`, `to_camera` | 11 |
| 21 | markdown | `## 8. Photon shot-noise SNR` | — | 0 |
| 22 | code | `def I_full(phi):` | `I_full`, `regime_of`, `SNR_shot_ideal`, `SNR_pixel`, `SNR_reselem_sim` | 0 |
| 23 | markdown | `## 9. Maximum non-destructive shots — threshold family [revised]` | — | 0 |
| 24 | code | `# Nmax_loss and Nmax_heat are defined in section 6.3 (the loss model); applied here.` | — | 0 |
| 25 | markdown | `### Figure — shot budget vs detuning for the candidate thresholds` | — | 0 |
| 26 | code | `Dg = np.linspace(0.5, 4.0, 400)` | — | 3 |
| 27 | markdown | `## 10. Spatial resolution` | — | 0 |
| 28 | code | `NA_e = (D_probe/2)/f1` | — | 0 |
| 29 | markdown | `## 11. Operating-point summary` | — | 0 |
| 30 | code | `# ---- operating-point report (folded from former section 13) ----` | `operating_point_report` | 0 |
| 31 | markdown | `## 12. Pulse duration and the $(\Delta, P, \tau)$ operating point [new]` | — | 0 |
| 32 | code | `# ===== CELL D: pulse duration & realistic destruction models =====` | `reabs_frac`, `Nmax_cleanloss`, `Nmax_heating` | 0 |
| 33 | markdown | `### 12.1 Pulse-duration sweep at the reference detuning` | — | 0 |
| 34 | code | `# ===== CELL F: tau sweep — SNR rises, DGI clears the read-noise floor...` | — | 8 |
| 35 | markdown | `### 12.2 $(\Delta, P, \tau)$ optimisation` | — | 0 |
| 36 | code | `# ===== CELL H: (Delta, P, tau) optimisation — operating map...` | `_pci_block`, `SNRres_fast` | 10 |
| 37 | markdown | `### 12.3 Recommendation and revised headline numbers` | — | 0 |
| 38 | code | empty cell | — | 0 |
| 39 | markdown | `---` | — | 0 |
| 40 | code | `# ============================================================================` | `tf_state`, `SNR_pixel_phi`, `run_sequence`, `accumulate` | 0 |
| 41 | markdown | `## 14. Evolution of the run` | — | 0 |
| 42 | code | `# ---- 14. Evolution of the run: condensate, temperature, phase, SNR ----------` | — | 12 |
| 43 | markdown | `## 15. The cloud fading on camera` | — | 0 |
| 44 | code | `# ---- 15. The cloud fading on camera: PCI frames across the run...` | `sim_image_state`, `camera_from_image` | 7 |
| 45 | markdown | `## 16. Detuning sweep: the SNR–destruction trade-off` | — | 0 |
| 46 | code | `# ============================================================================` | `snr_abs_pix`, `snr_pci_pix`, `usable_snr`, `usable_curve`, `make_strip` | 5 |
| 47 | markdown | `---` | — | 0 |
| 48 | markdown | `### The physical picture: why unequal σ± phases *rotate* the polarization [added]` | — | 0 |
| 49 | code | `# ============================================================================` | `theta_F_peak` | 3 |
| 50 | markdown | `### 17.2 Two detection schemes, and how the imaging simulation differs from PCI/DGI` | — | 0 |
| 51 | code | `# ============================================================================` | `sim_faraday_fields`, `faraday_maps` | 12 |
| 52 | markdown | `### 17.3 SNR, and why the destruction budget carries over unchanged` | — | 0 |
| 53 | code | `# ============================================================================` | `SNR_faraday_sim`, nested `_bin` | 0 |
| 54 | markdown | `### 17.4 Faraday imaging over the multi-shot run` | — | 0 |
| 55 | code | `# ---- re-image the SAME Sec-13 run (seq_h) with dual-port Faraday...` | `_bin_center` | 6 |
| 56-63 | mixed | Section 18 explanatory five-stage signal pipeline | — | 27 total |
| 64-95 | mixed | Section 19 step-by-step imaging walkthrough | — | 105 total |

### 8.2 Function-to-target-section mapping

| Target section | Existing functions to move or keep nearby |
|---|---|
| 02 Atomic Model | `_frac_residual`, `_tf_profile`, `tf_state` |
| 03 Light-Atom Interaction | `delta_of`, `phi_peak`, `od_resonant_equiv`, `intensity_at_atoms`, `N_scatt`, `reabs_frac`, `theta_F_peak` |
| 04 PCI | `sim_image` PCI branch, `blur_for_axis`, `_pci_block`, PCI-specific SNR calls |
| 05 DGI | `sim_image` DGI branch and DGI-specific camera/read-noise examples |
| 06 Faraday | `sim_faraday_fields`, `faraday_maps`, `SNR_faraday_sim` |
| 07 Camera | `N_phot_pix`, `to_camera`, `camera_from_image`, `_bin_center` |
| 08 Shot Noise | `I_full`, `regime_of`, `SNR_shot_ideal`, `SNR_pixel`, `SNR_reselem_sim`, `SNR_pixel_phi`, `SNRres_fast`, `snr_abs_pix`, `snr_pci_pix`, `usable_snr` |
| 09 Multi-shot Simulation | `Nmax_loss`, `Nmax_heat`, `Nmax_cleanloss`, `Nmax_heating`, `run_sequence`, `accumulate`, `sim_image_state`, `usable_curve`, `make_strip` |
| 10 Analysis | `operating_point_report` plus plotting-only cells for tables, maps, lineouts, sequence plots, and pedagogical figures |

## 9. Deliverables checklist

- Notebook audit report: complete in this file.
- Proposed notebook structure: provided in Section 5.
- Refactoring plan: provided in Section 6.
- Notebook code changes: none.
- Physics, numerical algorithms, parameters, and figures: unchanged.
