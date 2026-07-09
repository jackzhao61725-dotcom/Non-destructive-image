# Code-Notebook Pipeline Correspondence Audit

## Scope

This audit checks whether the current migrated code can reproduce the original
notebook computational pipeline when the helpers are used in the correct order.
It is not a plotting task and does not propose new physics.

The main finding is that many low-level notebook quantities are faithfully
migrated and tested, but the repository still lacks one canonical
notebook-aligned orchestration recipe for the full chain:

```text
parameters -> condensate -> projected column density -> phase/rotation
-> PCI/DGI/Faraday image -> camera image/noise -> sequence/SNR/destruction
```

This explains why tests can pass while broad generated figures still diverge
visually from the notebook: the helpers are mostly correct, but the figure
script can call them with different parameters, omit notebook display/camera
stages, or plot the wrong intermediate quantity.

## Current Public Helper Surface

The current package exports helpers for:

- atomic model: `build_thomas_fermi_state(...)`, `recoil_quantities(...)`;
- profiles: `thomas_fermi_profile_2d(...)`;
- light-atom interaction: `dimensionless_detuning(...)`,
  `scalar_phase_shift(...)`, `residual_optical_depth(...)`,
  `intensity_at_atoms(...)`, `scattered_photons_per_atom(...)`,
  `faraday_rotation_angle(...)`, `reabsorption_fraction(...)`;
- imaging: `propagate_scattered_field(...)`, `simulate_fourier_image(...)`,
  `simulate_pci_image(...)`, `simulate_dgi_image(...)`,
  `simulate_faraday_image(...)`;
- camera: `bin_to_camera_pixels(...)`, `add_camera_noise(...)`,
  `normalize_camera_counts(...)`, `simulate_camera_image(...)`,
  `simulate_noisy_camera_image(...)`;
- multi-shot: `simulate_multishot_sequence(...)`, `accumulate_snr(...)`;
- analysis/calibration: deterministic Faraday objective/sweeps and
  absorption-image observable helpers.

The helper surface is broad enough to reproduce many notebook arrays, but it is
not itself a complete notebook recipe.

## Stage-by-Stage Correspondence

### 1. Parameters and Constants

Notebook location:

- original notebook cells 1-6;
- `notebook_sections/00_imports.py`;
- `notebook_sections/01_parameters.py`.

What the stage does:

- defines constants, erbium transition parameters, probe geometry, condensate
  parameters, phase-plate parameters, camera parameters, and loss/scattering
  controls.

Key variables:

- `hbar`, `h`, `c`, `kB`, `amu`, `a0`;
- `m`, `lam`, `k`, `Gamma`;
- `P_probe_mW`, `D_probe`, `f1`, `f2`, `pix_cam`, `tau`;
- `N0`, `a_s`, `trap_Hz`, `T_cloud`;
- `t_p`, `theta`, `QE_cam`, `read_e`, `use_peak_intensity`, `eta_coll`;
- `sigma0`, `Isat`, `E_rec`, `T_rec`, `v_rec`, `E_phot`.

Current helper mapping:

- no single parameter object exists in `src/`;
- helpers accept constants and parameters explicitly;
- scripts/configs must supply notebook defaults.

Output identity:

- mostly scalars in SI units;
- some notebook display units are derived later, for example `cm^-3`,
  `cm^-2`, `um`, and `nK`.

Coverage:

- scalar notebook-output baseline exists in
  `regression/baseline/notebook_outputs.json`;
- there is no canonical source-of-truth config that contains all notebook
  defaults in one place for end-to-end reproduction.

Classification:

- **B. Helpers exist downstream, but no canonical notebook-aligned workflow
  config exists.**

### 2. Condensate / Thomas-Fermi State

Notebook location:

- original notebook cell 8;
- `notebook_sections/02_atomic_model.py`;
- Stage 18.1 in `notebook_sections/10_analysis.py`.

What the stage does:

- computes Thomas-Fermi condensate quantities from `N0`, `a_s`, `trap_Hz`,
  `m`, `hbar`, and `kB`.

Key variables:

- `omega`, `omega_bar`, `a_ho`, `mu`, `T_mu`;
- `n_peak`;
- `R`;
- `n_col`;
- `N_check`;
- later thermal bookkeeping: `N_tot_sc`, `Tc_sc`.

Current helper mapping:

- `build_thomas_fermi_state(...)` maps directly to `omega`, `omega_bar`,
  `a_ho`, `mu`, `T_mu`, `n_peak`, `R`, `n_col`, and `N_check`;
- `recoil_quantities(...)` maps to `E_rec`, `T_rec`, and `v_rec`.

Output identity:

- scalar and vector physical quantities;
- SI units internally;
- displayed as `nK`, `cm^-3`, `um`, and `cm^-2` in notebook figures.

Tests/baselines:

- notebook-output baseline covers derived scalar values;
- the notebook-aligned Stage 18.1 recovery compares direct notebook formulas
  with helper output for radii, column densities, profile, and column-density
  map;
- no general end-to-end test currently asserts that one canonical config
  recreates all condensate arrays used by later imaging sections.

Classification:

- **A for the helper formulas; B for canonical workflow orchestration.**

### 3. Projected Profile and Column Density

Notebook location:

- `_tf_profile(...)` in original notebook cell 18;
- Stage 18.1 and Step 19.1 in `notebook_sections/10_analysis.py`.

What the stage does:

- selects an imaging axis;
- selects the transverse plane;
- builds the 2D projected Thomas-Fermi profile:

```text
profile = max(0, 1 - a^2/Ra^2 - b^2/Rb^2)^1.5
```

- forms the absolute column-density map:

```text
column_density_map = n_col[axis] * profile
```

Current helper mapping:

- `thomas_fermi_profile_2d(...)` reproduces `_tf_profile(...)`;
- no helper takes a full `ThomasFermiState` plus `axis` and returns
  `gax`, grids, profile, and column-density map.

Output identity:

- `profile`: dimensionless 2D map, shape normally `(1024, 1024)`;
- `column_density_map`: absolute column density in `m^-2`;
- notebook plots convert it to `cm^-2` with `*1e-4`.

Tests/baselines:

- PCI/DGI and Faraday imaging baselines include `thomas_fermi_profile`,
  `grid_axis_m`, and `column_density_m2`;
- Stage 18.1 recovery records exact agreement between notebook direct formulas
  and helpers for the recovered figure.

Classification:

- **A for the profile helper and recovered Stage 18.1 script; B for reusable
  canonical workflow.**

### 4. Scalar Phase and Scattering

Notebook location:

- original notebook cells 10 and 14;
- `notebook_sections/03_light_atom_interaction.py`.

What the stage does:

- computes dimensionless detuning;
- computes scalar phase shift and residual OD;
- computes probe intensity at atoms;
- computes scattered photons per atom per pulse;
- computes clean-loss/heating/reabsorption quantities later used by multi-shot
  logic.

Key variables:

- `delta_of(Delta_Hz)`;
- `phi_peak(Delta_Hz, n_col_peak)`;
- `od_resonant_equiv(...)`;
- `intensity_at_atoms(P_mW)`;
- `N_scatt(...)`;
- `reabs_frac(...)`.

Current helper mapping:

- `dimensionless_detuning(...)` -> `delta_of`;
- `scalar_phase_shift(...)` -> `phi_peak`;
- `residual_optical_depth(...)` -> `od_resonant_equiv`;
- `intensity_at_atoms(...)` -> notebook intensity convention;
- `scattered_photons_per_atom(...)` -> `N_scatt`;
- `reabsorption_fraction(...)` -> `reabs_frac`.

Output identity:

- scalar values, not images;
- phase is radians;
- residual OD is dimensionless;
- scattered photons per atom is dimensionless per shot;
- reabsorption is a scalar averaged over principal axes.

Tests/baselines:

- helper regression status and notebook-output baseline cover these formulas;
- imaging baselines consume the scalar phase map;
- no single pipeline test currently verifies
  `condensate -> n_col[axis] -> phi_peak -> phase_map` as a named canonical
  recipe.

Classification:

- **A for formulas; B for canonical stage-to-stage recipe.**

### 5. PCI Imaging

Notebook location:

- original notebook cell 18, `sim_image(..., mode='PCI')`;
- notebook sections 7.3-7.4, 18.3-18.4, and 19.2;
- `notebook_sections/04_pci.py` and `notebook_sections/10_analysis.py`.

What the stage does:

- creates `object_field = exp(1j*phase_map)`;
- creates `scattered_field = object_field - 1`;
- propagates scattered field with `ifft2(fft2(scattered)*pupil)`;
- adds phase-plate reference `t_p*exp(1j*theta)`;
- returns `I_pci = abs(E)^2`.

Current helper mapping:

- `propagate_scattered_field(...)` maps to the FFT/pupil operation;
- `simulate_fourier_image(...)` maps to shared scattered-field propagation and
  reference recombination;
- `simulate_pci_image(...)` maps to PCI reference construction and intensity.

Output identity:

- ideal image intensity array in `I/I0`;
- usually `(1024, 1024)` for notebook baseline grid;
- display limits and camera binning are separate from the helper output.

Tests/baselines:

- `pci_dgi_imaging_baseline_v1.npz` includes phase map, object field,
  scattered field, propagated scattered field, PCI reference, and PCI image;
- `test_pci_orchestration.py` compares helper outputs to baseline arrays.

Classification:

- **A for ideal PCI image core; B/C for notebook-style camera/display workflows.**

### 6. DGI Imaging

Notebook location:

- original notebook cell 18, `sim_image(..., mode='DGI')`;
- sections 7.3-7.4, 18.3-18.4, and 19.3;
- `notebook_sections/05_dgi.py` and `notebook_sections/10_analysis.py`.

What the stage does:

- uses the same scalar phase map and propagated scattered field as PCI;
- replaces the phase-plate reference with attenuated carrier
  `10**(-OD/2)`;
- returns `I_dgi = abs(E)^2`.

Current helper mapping:

- `simulate_dgi_image(...)` maps to DGI orchestration above
  `simulate_fourier_image(...)`.

Output identity:

- ideal DGI image intensity in `I/I0`;
- small absolute values near a dark floor;
- display limits are important and are not part of the helper.

Tests/baselines:

- `pci_dgi_imaging_baseline_v1.npz` includes DGI reference and DGI image;
- `test_dgi_orchestration.py` compares helper output to baseline arrays.

Classification:

- **A for ideal DGI image core; B/C for notebook-style camera/display workflows.**

### 7. Faraday Imaging

Notebook location:

- original notebook cells 49, 51, 53, 55, and sections 17-19;
- `notebook_sections/06_faraday.py`;
- `notebook_sections/10_analysis.py`.

What the stage does:

- computes `theta_F = kappa_F * phi_peak`;
- current `kappa_F = 1.0` is a placeholder calibration convention;
- propagates sigma+ and sigma- fields with opposite phase maps;
- recombines into `Ex` and `Ey`;
- computes dark-field intensity and dual-port outputs.

Current helper mapping:

- `faraday_rotation_angle(...)` -> `theta_F_peak`;
- `simulate_faraday_image(...)` -> circular propagation, `Ex/Ey`
  recombination, `I_dark`, `I_u`, `I_v`, and `S`.

Output identity:

- `theta_f_map_rad`: 2D rotation map in radians;
- circular/object/propagated fields: complex maps;
- dark-field and port outputs: intensity maps;
- dual-port signal: dimensionless normalised difference.

Tests/baselines:

- `faraday_imaging_baseline_v1.npz` includes intermediate fields and outputs;
- `test_faraday_orchestration.py` compares helper output to baseline arrays.

Classification:

- **A for ideal Faraday field/image core; C for camera-level dual-port frame
  workflows and multi-shot Faraday sequence.**

### 8. Camera and Noise

Notebook location:

- original notebook cell 20, `to_camera(...)`;
- section 15 fading camera frames;
- section 19 camera steps;
- `notebook_sections/07_camera.py`.

What the stage does:

- bins a high-resolution image to camera pixels using 15x15 averaging;
- converts binned normalised intensity into photons/electrons;
- applies Poisson photon noise and Gaussian read noise;
- divides by `Nd` to return normalised `I/I0`-like camera units.

Current helper mapping:

- `bin_to_camera_pixels(...)` -> notebook reshape/mean binning;
- `add_camera_noise(...)` -> Poisson plus read noise;
- `normalize_camera_counts(...)` -> `counts/Nd`;
- `simulate_camera_image(...)` -> deterministic binning/normalisation;
- `simulate_noisy_camera_image(...)` -> stochastic camera recipe with explicit
  RNG.

Output identity:

- binned deterministic image: normalised image units;
- noisy counts: electrons if not normalised;
- noisy image: normalised camera image if `normalize=True`.

Tests/baselines:

- camera tests verify shape, binning equivalence, deterministic
  normalisation, explicit RNG reproducibility, and direct helper composition;
- no canonical end-to-end test asserts notebook-specific values for
  `N_phot_pix`, `pix_obj`, `Mag`, or a specific camera frame from the notebook.

Classification:

- **A for camera helper recipe; B/C for notebook-exact camera workflow because
  parameter defaults and RNG state are not centralised.**

### 9. SNR, Destructiveness, and Optimisation

Notebook location:

- original notebook sections 8, 9, 12, 16, and 17.3;
- `notebook_sections/08_shot_noise.py`;
- `notebook_sections/10_analysis.py`.

What the stage does:

- computes ideal and realistic SNR estimates;
- computes shot budgets and heating/loss constraints;
- builds operating maps over detuning, power, and exposure time;
- compares PCI/DGI/Faraday readout efficiency under a shared destruction
  budget.

Current helper mapping:

- `evaluate_faraday_operating_point(...)`;
- `sweep_faraday_detuning(...)`;
- `sweep_faraday_intensity(...)`;
- `sweep_faraday_exposure_time(...)`;
- `summarise_faraday_sweep(...)`.

Output identity:

- deterministic scalar metrics and arrays from one-dimensional sweeps;
- not notebook display maps;
- not stochastic/noise-averaged objectives;
- not full notebook operating maps.

Tests/baselines:

- regression tests cover deterministic objective keys, finite outputs,
  one-dimensional sweep behaviour, best-index selection, and summary helper;
- no baseline covers notebook full 2D operating maps or absorption/dispersive
  crossover plots.

Classification:

- **C. Partially migrated. Notebook analysis logic remains local.**

### 10. Multi-Shot / Continuous Imaging

Notebook location:

- original notebook cell 40, `run_sequence(...)`;
- sections 14-16 and 19.6;
- `notebook_sections/09_multishot_simulation.py`;
- `notebook_sections/10_analysis.py`.

What the stage does:

- steps frame by frame;
- updates atom number or temperature;
- recomputes Thomas-Fermi state and phase for each frame;
- computes deterministic SNR;
- optionally renders noisy frame sequences and filmstrips.

Current helper mapping:

- `simulate_multishot_sequence(...)` maps to deterministic bookkeeping;
- `accumulate_snr(...)` maps to RMS accumulated SNR.

Output identity:

- arrays for shot index, `N0`, condensate fraction/loss, temperature, phase,
  SNR;
- no noisy image frames by itself;
- no Faraday dual-port sequence by itself.

Tests/baselines:

- tests verify clean-loss hand calculation, heating update, output lengths,
  deterministic repeated calls, and accumulated-SNR convention;
- no canonical end-to-end test links
  `N0 -> TF state -> n_col -> phase -> imaging -> camera -> SNR` across frames.

Classification:

- **A for deterministic sequence core; C for image/frame sequence workflows.**

## Missing Canonical Orchestration Paths

The table below classifies each path.

| Path | Status | Reason |
| --- | --- | --- |
| parameters -> condensate | B | Helper exists, but no single canonical notebook-default parameter object exists. |
| condensate -> column density | B | Helper and recovery script exist; no reusable recipe from `ThomasFermiState + axis + grid`. |
| column density -> phase shift | B | Scalar helper exists, but no canonical phase-map workflow binding axis/grid/profile/default detuning. |
| phase shift -> PCI image | A/B | Core helper and baseline exist; workflow needs canonical parameter/grid recipe for notebook figures. |
| phase shift -> DGI image | A/B | Core helper and baseline exist; workflow needs notebook display/camera conventions. |
| phase/Faraday rotation -> Faraday image | A/C | Core helper and baseline exist; dual-port camera and sequence logic remain notebook-local. |
| ideal image -> camera image | B/C | Helpers exist; notebook `N_phot_pix`, global RNG, and mode-specific camera recipes are not centralised. |
| camera/image quantities -> multi-shot sequence | C | Deterministic sequence core exists, but full frame regeneration and SNR callbacks remain external. |

## Parameter and Default Mismatch Risks

A helper call can be mathematically correct and still not reproduce a notebook
figure if any of these differ:

- atom number: notebook uses `N0 = 2.5e4`;
- grid: notebook imaging grid uses `Ngrid = 1024`, `FOV = 100e-6`;
- imaging axis: notebook reference figures commonly use `axis = 0` and a `y,z`
  transverse plane;
- units: helpers use SI, notebook display often uses `um`, `cm^-3`, `cm^-2`,
  `nK`, or normalised `I/I0`;
- detuning: many reference figures use `1.5 GHz`, but not all;
- probe power: PCI/DGI/Faraday examples use different powers in different
  notebook sections;
- exposure time: notebook defaults and later sections differ;
- `kappa_F`: Faraday uses `1.0` as a placeholder, not a calibrated value;
- camera parameters: `QE_cam`, `read_e`, magnification, `pix_obj`, and
  `N_phot_pix` must match;
- normalisation: plotting a dimensionless profile is not the same as plotting
  absolute column density or camera-normalised intensity;
- stochastic assumptions: notebook RNG state depends on prior draws, while
  helpers require explicit RNG input.

## Output Identity Rules

For future recovery work, each output must state what it is:

- `ThomasFermiState`: scalar/vector SI physical quantities;
- `profile`: dimensionless 2D projected TF profile;
- `column_density_map`: `m^-2` physical map, displayed as `cm^-2`;
- `phase_map`: radians;
- `theta_f_map`: radians;
- `object_field`, `scattered_field`, propagated fields: complex fields;
- PCI/DGI/Faraday intensity images: ideal normalised intensity maps;
- binned image: deterministic camera-pixel image;
- noisy counts: electron counts;
- noisy normalised image: camera-normalised image after stochastic noise;
- multi-shot sequence arrays: deterministic per-frame state values.

Display transforms are separate from physical output:

- axis unit conversion;
- image zoom;
- colour limits;
- colour map choice;
- lineout direction;
- normalisation for visual comparison;
- stochastic sampling.

The removed broad figure gallery blurred these layers.

## Tests and Baselines

Current tests verify:

- scalar notebook baseline status;
- atomic/light-atom helper behaviour;
- PCI/DGI baseline existence, shapes, keys, and reproducibility;
- PCI and DGI orchestration against baseline arrays;
- Faraday baseline existence, shapes, keys, and reproducibility;
- Faraday orchestration against baseline arrays;
- deterministic camera pipeline;
- stochastic camera helper with explicit RNG;
- deterministic multi-shot core;
- Faraday optimisation objective/sweeps/summary;
- absorption calibration helper basics.

Current baselines:

- `regression/baseline/notebook_outputs.json`;
- `regression/baseline/imaging/pci_dgi_imaging_baseline_v1.npz`;
- `regression/baseline/imaging/faraday_imaging_baseline_v1.npz`.

What is missing:

- canonical notebook-default config test;
- end-to-end pipeline test from parameters to Stage 18.1 arrays;
- phase-map recovery test from condensate to `phase_map`;
- PCI/DGI/Faraday recipe tests that include canonical grid, axis, and display
  metadata;
- camera test that reproduces notebook `N_phot_pix` and one fixed-seed frame
  from a canonical recipe;
- multi-shot test that recomputes `N0 -> TF state -> phase -> SNR` across
  frames with notebook parameters;
- tests that compare canonical intermediate arrays, not only individual helper
  calls.

## Diagnosis of Failed Figure Generation

The failed broad gallery can diverge even with all tests passing because:

1. It bypassed the notebook's staged pipeline and began from independently
   chosen plotting recipes.
2. It likely used helper defaults or new config values that did not match the
   notebook's figure-specific defaults.
3. It plotted helper outputs without the notebook's display transforms:
   units, zoom, colour limits, lineout direction, or background reference.
4. It combined multiple stages before proving earlier intermediate arrays
   matched.
5. It treated the helper output as the figure, even where the notebook figure
   included camera binning, noise, ratio construction, or explanatory phasors.
6. It did not preserve stochastic notebook state for noisy camera views.
7. It did not use comparison reports for intermediate arrays before plotting.

In short: the helpers can be correct, but the workflow can still be wrong.

## Recommended Minimal Repair Plan

Do not add broad plotting yet. The minimal repair path is:

1. Define a canonical notebook-default config.
   - Include constants, apparatus, condensate, grid, phase-plate, camera,
     Faraday, and multi-shot defaults.
   - Keep all values explicit.

2. Define canonical notebook-aligned workflow recipes.
   - Start with scripts/configs rather than new `src/` APIs.
   - Recipes should produce named intermediate arrays and comparison reports.

3. Recover one stage at a time.
   - Stage 18.1 condensate/density is already recovered.
   - Next: Step 19.1 scalar phase map.
   - Then: one PCI ideal image.
   - Then: camera binning/noise for PCI.
   - Only after that: DGI, Faraday, and multi-shot figures.

4. Add comparison reports before figures.
   - Shape, min, max, mean, sum, centre value, peak index;
   - selected lineout samples;
   - max absolute and relative differences against notebook-direct logic.

5. Add focused regression tests after recipes stabilise.
   - Test canonical intermediate arrays first.
   - Avoid testing only SVG output.

6. Defer new `src/` orchestration functions.
   - Scripts/configs are enough while recovering notebook alignment.
   - Add small `src/` orchestration only if repeated recipes become stable and
     clearly belong in the public simulator API.

## Are Code Changes Required?

Not immediately.

The next repair step should be scripts/configs/documentation plus comparison
reports, not simulator-code changes. Current helpers are largely sufficient for
the next one or two recovery milestones, provided they are called in the exact
notebook order with exact notebook parameters and display metadata.

Small `src/` orchestration functions may become useful later, but adding them
now risks freezing a workflow before the notebook correspondence has been
proven end to end.

## Bottom Line

The current migrated code can reproduce many notebook computational stages when
used correctly, but the repository does not yet provide a single canonical
notebook-aligned pipeline recipe. That missing orchestration layer is the main
gap, not an obvious failure of the individual helpers.
