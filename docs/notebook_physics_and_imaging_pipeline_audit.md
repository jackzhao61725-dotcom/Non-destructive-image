# Full Notebook Physics and Imaging Pipeline Audit

## Scope

This audit re-reads the original notebook as the most complete historical
Version 1 computational prototype. It does not treat the notebook as final
calibrated theory, but it does treat the notebook's computational flow as the
reference that any notebook-aligned figure recovery must understand before
plotting begins.

The key conclusion is that figure recovery must start from the condensate
object, not from image plotting. Most visible notebook figures are downstream
views of a specific Thomas-Fermi condensate, projected along a chosen imaging
axis, converted into phase or rotation maps, propagated through a finite-NA
Fourier system, then optionally binned and noised by a camera model.

## 1. Condensate Generation / Physical Object

The notebook defines physical constants first, then apparatus and condensate
parameters, and only then builds the condensate object. The first physical
object is not an image. It is a Thomas-Fermi Bose-Einstein condensate specified
by atom number, scattering length, trap frequencies, mass, and fundamental
constants.

The main inputs are defined in the early notebook cells and exported in
`notebook_sections/01_parameters.py`:

- atomic mass: `m = 166 * amu`;
- wavelength and wavevector: `lam = 401e-9`, `k = 2*pi/lam`;
- natural linewidth: `Gamma = 2*pi*29.5e6`;
- probe geometry and imaging arm: `P_probe_mW`, `D_probe`, `f1`, `f2`,
  `pix_cam`, `tau`;
- condensate inputs: `N0 = 2.5e4`, `a_s = 72*a0`,
  `trap_Hz = [293, 14, 233] Hz`, `T_cloud = 200 nK`;
- phase-plate and camera parameters: `t_p`, `theta`, `QE_cam`, `read_e`;
- loss/scattering controls: `use_peak_intensity`, `eta_coll`.

The Thomas-Fermi quantities are then computed in the notebook's section 4 and
exported through `notebook_sections/02_atomic_model.py`:

- trap angular frequencies:
  `omega = 2*pi*trap_Hz`;
- geometric mean frequency:
  `omega_bar = (omega.prod())**(1/3)`;
- harmonic oscillator length:
  `a_ho = sqrt(hbar/(m*omega_bar))`;
- chemical potential:
  `mu = 0.5*(15*N0*a_s/a_ho)**(2/5) * hbar*omega_bar`;
- chemical potential temperature:
  `T_mu = mu/kB`;
- peak 3D density:
  `n_peak = mu*m/(4*pi*hbar**2*a_s)`;
- Thomas-Fermi radii:
  `R = sqrt(2*mu/(m*omega**2))`;
- peak column densities along the three principal imaging axes:
  `n_col = (4/3)*n_peak*R`;
- atom-number consistency check:
  `N_check = (8*pi/15)*n_peak*R.prod()`.

The vector `R` is the physical cloud size in the three trap axes. The vector
`n_col` is critical: later scalar phase shifts and Faraday rotations are usually
computed from one component of `n_col`, selected by the imaging axis. For the
main across-cigar imaging path, the notebook uses `axis = 0`, so the image plane
is the transverse `y,z` plane and the peak column density is `n_col[0]`.

The notebook also computes a self-consistent total atom number and critical
temperature for the thermal fraction:

- solve `Ntot*f - N0 = 0`;
- `Tc = 0.94*hbar*omega_bar/kB * Ntot**(1/3)`;
- condensate fraction `f = 1 - (T_cloud/Tc)**3`.

Those quantities do not replace the fresh condensate object used for imaging,
but they are later used in heating and multi-shot bookkeeping.

## 2. Density, Profile, and Normalisation Logic

The notebook uses both 3D density and 2D projected density, but at different
stages.

The 3D central density logic is used to derive and summarise the condensate:

- `n_peak` is an absolute peak 3D density in SI units, often printed in
  `cm^-3`;
- 1D density cuts in later explanatory figures use
  `n_peak * max(0, 1 - (x/R_i)**2)`.

The imaging calculation uses projected 2D column density. The repeated notebook
profile helper is:

```text
profile(a,b) = max(0, 1 - a^2/Ra^2 - b^2/Rb^2)^1.5
```

This is the normalised 2D Thomas-Fermi column-density shape in the plane
transverse to the imaging axis. The absolute column-density map is:

```text
column_density_map = n_col[axis] * profile
```

The exponent `1.5` is important: a helper-generated figure using a 2D
Thomas-Fermi exponent, a Gaussian, a flattened profile, or a different axis
will not visually match the notebook.

The grid for the main Fourier imaging sections is:

- `Ngrid = 1024`;
- `FOV = 100e-6 m`;
- `dgrid = FOV/Ngrid`;
- `gax = (arange(Ngrid)-Ngrid//2)*dgrid`;
- `GA, GB = meshgrid(gax, gax)`.

The plot extent is usually:

```text
ext = [-FOV/2*1e6, FOV/2*1e6, -FOV/2*1e6, FOV/2*1e6]
```

Most image display panels then zoom to `xlim(-45,45)` and `ylim(-12,12)`.
For across-cigar imaging, the horizontal plotted axis is usually `y (um)` and
the vertical plotted axis is `z (um)`.

Some notebook figures show absolute density or column density with units:

- 3D density line cuts use `cm^-3`;
- column density maps use `cm^-2`;
- phase maps use radians.

Other notebook figures show normalised image intensity, camera-normalised
`I/I0`, or normalised signal `S`. A recovery script must preserve which type is
being displayed. A plot of the profile alone is not equivalent to a notebook
camera image.

## 3. Light-Atom Interaction

The notebook scalar phase shift comes from the selected peak column density.
The dimensionless detuning is:

```text
delta = 2 * Delta_Hz * 2*pi / Gamma
```

The scalar phase shift is:

```text
phi_peak = sigma0 * n_col_peak * delta / (2*(1 + delta^2))
```

The residual optical depth is:

```text
OD_residual = sigma0 * n_col_peak / (1 + delta^2)
```

The phase map used by PCI and DGI is:

```text
phase_map = phi_peak * profile
```

Here `phi_peak` is a scalar and `phase_map` is a 2D array. The scalar is derived
from the condensate's selected peak column density and the detuning. The profile
is derived from the condensate radii and the imaging axis.

Scattering and destructiveness use probe intensity and pulse duration, not the
2D density profile:

```text
I_atoms = 2 * (P_mW*1e-3)/(pi*(D_probe/2)^2)
N_gamma = (Gamma/2) * s/(1 + s + delta^2) * tau
```

where `s = I_atoms/Isat`. The factor of two is used when
`use_peak_intensity=True`, because the atoms are taken to sit at the Gaussian
beam centre. Changing this convention changes destructiveness and camera photon
counts.

Reabsorption later uses the three principal column densities:

```text
OD_i = sigma0*n_col[i]/(1 + delta^2)
reabs_frac = mean(1 - exp(-OD_i))
```

This is a scalar correction for heating/destruction, not an image map.

## 4. PCI Pipeline

The notebook's PCI path is built around a pure scalar phase object. For a chosen
imaging axis:

```text
profile = max(0, 1 - GA^2/Ra^2 - GB^2/Rb^2)^1.5
object_field = exp(1j * phi_peak * profile)
scattered = object_field - 1
propagated_scattered = ifft2(fft2(scattered) * pupil)
reference = t_p * exp(1j*theta)
E_pci = reference + propagated_scattered
I_pci = abs(E_pci)^2
```

The Fourier pupil is:

```text
NA = (D_probe/2)/f1
fx = fftfreq(Ngrid, dgrid)
pupil = sqrt(FX^2 + FY^2) <= NA/lam
```

The notebook does not propagate the full carrier through the FFT. It propagates
the scattered component and adds a mode-specific carrier/reference field after
the finite-NA propagation. This is a common source of mismatch if a recovery
script propagates `exp(i phi)` as a whole and treats the carrier numerically.

The final PCI image array is a normalised intensity `I/I0`. The no-atom
phase-plate background is `bg_plate = t_p**2`. Many PCI plots use fixed display
limits near `vmin=0.85`, `vmax=1.15`, and lineouts include a horizontal
`t_p^2` background reference.

Current helper mapping:

- `simulate_pci_image(...)` correctly forms `object_field`, the PCI reference
  field, the propagated scattered field, and the intensity.
- It is regression-tested against
  `regression/baseline/imaging/pci_dgi_imaging_baseline_v1.npz`.

## 5. DGI Pipeline

DGI uses the same scalar phase map, the same scattered field, and the same
finite-NA propagation as PCI. The contrast mechanism differs only in the
reference/carrier:

```text
reference = 10**(-OD/2)
E_dgi = reference + propagated_scattered
I_dgi = abs(E_dgi)^2
```

In the notebook, the default dark-stop optical depth is `OD = 4.0`, giving a
small reference amplitude. DGI therefore appears as a faint bright signal on a
dark floor. It is quadratic-like in signal strength and is often displayed with
limits such as `vmin=-0.01`, `vmax=0.05`.

A recovery plot that normalises DGI to its own maximum can look visually clear
but will not match the notebook's physical framing. Notebook DGI panels are
normally in `I/I0` units with small absolute contrast.

Current helper mapping:

- `simulate_dgi_image(...)` correctly mirrors the DGI orchestration above
  `simulate_fourier_image(...)`.
- It is regression-tested against the same PCI/DGI baseline.

## 6. Faraday Pipeline

The notebook's current Faraday model is phenomenological:

```text
theta_F_peak = kappa_F * phi_peak
kappa_F = 1.0
```

The notebook explicitly does not implement a microscopic circular-transition
model. The scalar phase lineshape is reused and scaled by `kappa_F`. This is a
placeholder calibration convention, not a fitted experimental result.

The Faraday image path differs from PCI/DGI because it propagates two circular
components with opposite phase:

```text
Pp = 1 + ifft2(fft2(exp(+1j*theta_F*profile) - 1) * pupil)
Pm = 1 + ifft2(fft2(exp(-1j*theta_F*profile) - 1) * pupil)
Ex = (Pp + Pm)/2
Ey = 1j*(Pp - Pm)/2
```

The dark-field output is:

```text
I_dark = abs(Ey)^2
```

The dual-port outputs are:

```text
I_u = abs(Ex + Ey)^2 / 2
I_v = abs(Ex - Ey)^2 / 2
S = (I_v - I_u)/(I_v + I_u)
```

The notebook often displays:

- dark-field Faraday as `I_dark/I0` with dark-field colour limits;
- dual-port Faraday as the normalised difference `S`, with a diverging colour
  map and fixed limits around `[-0.5, 0.5]`;
- noisy dual-port `S_map` computed after applying camera noise separately to
  the two ports.

Current helper mapping:

- `faraday_rotation_angle(...)` preserves `theta_F = kappa_F * phi_peak`;
- `simulate_faraday_image(...)` preserves the opposite circular phase maps,
  finite-NA propagation, `Ex/Ey` recombination, dark-field intensity, dual-port
  intensities, and dual-port signal.
- Faraday arrays are regression-tested against
  `regression/baseline/imaging/faraday_imaging_baseline_v1.npz`.

Likely mismatch risk:

- helper-based plots may show `theta_F_map`, `I_dark`, and `S` side by side
  without notebook display limits, camera binning, zoom, or noisy-port
  construction;
- dual-port camera figures must noise `I_u` and `I_v` first, then form `S`;
  noising an already formed `S` map is not notebook-equivalent.

## 7. Camera and Noise

The notebook converts ideal image arrays into camera-level arrays using a fixed
binning recipe and stochastic photon/read noise.

The magnification and object pixel size are:

```text
Mag = f2/f1
pix_obj = pix_cam/Mag
```

Detected photons per camera pixel at `I/I0 = 1` are:

```text
N_phot_pix = intensity_at_atoms(P_mW) * pix_obj^2 * tau * QE / E_phot
```

The notebook binning recipe is:

```text
nb = (Ngrid//15)*15
binned = Iratio[:nb,:nb].reshape(nb//15,15,nb//15,15).mean(axis=(1,3))
```

The stochastic camera recipe is:

```text
counts = poisson(clip(binned,0,None)*Nd) + normal(0, read_e, binned.shape)
camera_image = counts/Nd
```

The deterministic part is the high-resolution image and the binned noiseless
image. The stochastic part is the Poisson draw and Gaussian read noise. The
notebook initializes `rng = np.random.default_rng(12345)` near the start, so
figure appearance can depend on prior random draws and cell execution order.

Current helper mapping:

- `bin_to_camera_pixels(...)` matches the notebook binning;
- `simulate_camera_image(...)` covers deterministic binning and normalisation;
- `simulate_noisy_camera_image(...)` applies the stochastic recipe with an
  explicit RNG, which is better for testing than notebook-global state.

For notebook-aligned recovery, stochastic figures should either use saved
notebook arrays or reproduce the notebook RNG state carefully. Otherwise the
visual noise pattern will differ even if the physics is correct.

## 8. Shot Noise, SNR, and Destructiveness

The notebook uses several SNR layers:

- ideal PCI shot-noise SNR:
  `2*t_p*phi*sqrt(Nphot)`;
- realistic PCI pixel SNR:
  contrast relative to `bg_plate`, with photon noise and read noise;
- realistic resolution-element SNR:
  sum over a small block of binned pixels;
- DGI and Faraday SNRs read from simulated image intensities or dual-port
  error propagation.

The notebook repeatedly stresses that signal extraction and destructiveness are
not the same object. Scattered photons per atom are computed from detuning,
power, and exposure. Image SNR is computed from the phase/rotation image,
camera photon budget, and read noise. Destructiveness limits the maximum number
of useful frames.

Two destruction models appear:

- clean loss:
  `N0(s) = N0*exp(-eta_coll*N_gamma*s)`;
- heating:
  `T = (T**4 + dE/A_E)**0.25`, then
  `N0_now = N_tot_sc*(1 - (T/Tc)**3)`.

The notebook's optimisation and trade-off plots are therefore not just
brightness plots. They relate per-frame signal or SNR to a finite destruction
budget.

## 9. Multi-Shot / Continuous Imaging

The notebook multi-shot engine is defined in section 13. It starts from the
fresh Thomas-Fermi condensate and updates the condensate state frame by frame.

For each frame, the notebook:

1. computes the current condensate atom number;
2. recomputes Thomas-Fermi state quantities through `tf_state(N0_now)`;
3. selects the current peak column density;
4. computes the current peak phase;
5. computes the current deterministic SNR;
6. stops when the condensate-loss fraction reaches the chosen limit.

The fresh-cloud state is therefore not fixed across the whole run. The cloud
shrinks and the peak phase changes as `N0` changes. The section 15 fading-camera
figure explicitly recomputes the radii and phase for selected frames and then
passes them through the same PCI imaging and camera recipe.

The notebook has both deterministic sequence plots and stochastic frame plots:

- deterministic: condensate fraction, temperature, phase, per-shot SNR,
  accumulated SNR;
- stochastic: noisy PCI frames, noisy Faraday ports, noisy filmstrips.

Current helper mapping:

- `simulate_multishot_sequence(...)` captures the deterministic sequence core
  and accepts callbacks for phase and SNR;
- it does not migrate noisy frame rendering, Faraday dual-port frame sequence,
  detuning maps, operating maps, or filmstrips.

## 10. Notebook Figure Families

The figure families below should be treated as downstream views of the
condensate pipeline.

### Condensate / Density Figures

Notebook locations:

- section 4 prints derived TF quantities;
- section 18.1 plots 3D density line cuts and a column-density map.

Plotted quantities:

- 3D central density cuts `n(r)` in `cm^-3`;
- column density map `n_col[axis] * profile` in `cm^-2`.

Upstream variables:

- `N0`, `a_s`, `trap_Hz`, `m`, `mu`, `n_peak`, `R`, `n_col`;
- imaging axis determines which `n_col` and which transverse radii appear.

Recovery status:

- can be recovered from current helpers, but the first recovery target must
  match notebook variables and units before any image figures are attempted.

### Phase / Rotation Figures

Notebook locations:

- sections 5, 17.1, 18.2, and 19.1;
- phase map and lineout in step 1.

Plotted quantities:

- scalar phase map `phi_peak * profile`;
- Faraday rotation map `theta_F_peak * profile`;
- lineouts along the notebook axes.

Recovery status:

- can be recovered from current helpers if axis, grid, units, extent, and
  colour limits match the notebook.

### Fourier / Wavelet / Pupil Diagnostic Figures

Notebook locations:

- section 19.1, steps 2-5.

Plotted quantities:

- `object_field = exp(i phi_map)`;
- scattered wavelet `w = object_field - 1`;
- Fourier power before/after pupil;
- clear image without phase plate or dark stop.

Recovery status:

- mostly notebook-section logic, not current high-level helper outputs;
- scientifically valuable for explanation, but should not be regenerated until
  the phase map and FFT arrays are numerically matched.

### PCI Figures

Notebook locations:

- sections 7.3-7.4;
- section 18.3-18.4;
- section 19.2, steps 6-7;
- section 15 for fading PCI frames.

Plotted quantities:

- ideal PCI image `I_pci`;
- binned noiseless camera image;
- one noisy camera frame;
- lineouts against `t_p^2` background.

Recovery status:

- core image helper is available and tested;
- notebook-aligned visual recovery still needs exact axis selection, grid,
  camera binning, colour limits, zoom window, and RNG handling.

### DGI Figures

Notebook locations:

- sections 7.3-7.4;
- section 18.3-18.4;
- section 19.3, steps 8-9.

Plotted quantities:

- DGI ideal image `I_dgi`;
- binned noiseless and noisy camera images;
- dark-floor lineouts.

Recovery status:

- core helper is available and tested;
- visual recovery must preserve small absolute `I/I0` scale and notebook
  display limits.

### Faraday Figures

Notebook locations:

- sections 17.1-17.4;
- sections 18.2-18.4;
- section 19.4-19.5, steps 10-13.

Plotted quantities:

- `theta_F` table and maps;
- dark-field Faraday image;
- dual-port `I_u`, `I_v`, and `S`;
- noisy dual-port frames;
- extracted angle robustness under common-mode probe flicker;
- Faraday SNR over the same multi-shot sequence.

Recovery status:

- core helper is available and tested;
- camera-level dual-port figures require notebook-style port-specific camera
  noise and ratio construction.

### SNR, Destruction, and Optimisation Figures

Notebook locations:

- sections 8, 9, 12, 16;
- section 17.3 for mode comparison;
- section 18.4 summary table.

Plotted quantities:

- SNR per pixel or resolution element;
- shot budget vs detuning;
- pulse-duration sweep;
- operating map in power/exposure space;
- detuning sweep with absorption/dispersive crossover;
- mode comparison table.

Recovery status:

- some deterministic optimisation helpers exist;
- many notebook plots remain notebook-local because they combine camera
  binning, read-noise floors, heating budgets, and image-specific SNR
  estimates.

### Multi-Shot / Filmstrip Figures

Notebook locations:

- sections 14-16 and section 19.6.

Plotted quantities:

- condensate fraction and temperature vs frame;
- phase and SNR vs frame;
- accumulated SNR;
- selected noisy frames or filmstrips.

Recovery status:

- deterministic sequence core is migrated;
- noisy frame sequences and Faraday dual-port sequence rendering remain
  notebook-local.

## 11. Mapping to Current Modular Code

### Faithfully Migrated and Tested

- `build_thomas_fermi_state(...)`: notebook TF algebra for `omega`,
  `omega_bar`, `a_ho`, `mu`, `n_peak`, `R`, `n_col`, and `N_check`.
- `thomas_fermi_profile_2d(...)`: notebook projected profile
  `max(0, 1-a^2/Ra^2-b^2/Rb^2)^1.5`.
- `scalar_phase_shift(...)`, `residual_optical_depth(...)`,
  `scattered_photons_per_atom(...)`, `reabsorption_fraction(...)`: notebook
  light-atom scalar formulas.
- `faraday_rotation_angle(...)`: notebook phenomenological
  `theta_F = kappa_F*phi_peak`.
- `simulate_fourier_image(...)`: notebook scattered-field FFT/pupil convention.
- `simulate_pci_image(...)` and `simulate_dgi_image(...)`: thin
  notebook-equivalent orchestration.
- `simulate_faraday_image(...)`: notebook Faraday circular-component
  propagation and recombination.
- `simulate_camera_image(...)` and `simulate_noisy_camera_image(...)`:
  notebook binning/noise recipe with explicit RNG.
- `simulate_multishot_sequence(...)` and `accumulate_snr(...)`:
  deterministic multi-shot bookkeeping.

### Partially Migrated

- Figure-level camera framing, colour limits, axis labels, zoom windows, and
  notebook-specific lineout choices are not fully captured by helpers.
- Multi-shot helpers do not regenerate shrinking-cloud noisy frames by
  themselves; they require callbacks and additional notebook logic.
- Faraday camera-level dual-port ratio figures require separate port noising.
- Optimisation helpers cover deterministic scalar one-dimensional sweeps, not
  the full notebook operating maps or absorption/dispersive crossover figure.

### Notebook-Local

- explanatory phasor diagrams;
- Fourier-plane diagnostic plots;
- clear-image "invisible phase object" demonstration;
- noisy PCI/DGI/Faraday filmstrips;
- dual-port flicker-cancellation demonstration;
- absorption/dispersive crossover map;
- full two-dimensional operating maps;
- notebook narrative figure styling and axis/colour choices.

## 12. Why Previous Helper-Based Figures Diverged

The removed broad gallery likely diverged for several reasons:

1. It started from helper outputs and independent plotting assumptions rather
   than reconstructing the notebook's physical object and display chain first.
2. It used representative configuration values that did not necessarily match
   the notebook's canonical figure-specific choices, such as `N0 = 2.5e4`,
   `Ngrid = 1024`, `FOV = 100 um`, `axis = 0`, and the notebook's zoom window.
3. It did not preserve every notebook display convention: fixed colour limits,
   zoom range, unit conversions, lineout direction, background references, and
   camera binning are visually important.
4. It compressed multiple figure families into a broad gallery at once, making
   it difficult to identify whether mismatch came from condensate generation,
   phase, FFT/pupil, camera, or plotting.
5. It treated existing helpers as sufficient for finished figures, but helpers
   intentionally preserve core arrays, not notebook-specific explanatory
   plotting logic.
6. Stochastic camera figures cannot visually match the notebook unless the RNG
   state and order of random draws are reproduced or fixed output arrays are
   compared.

## 13. Recommended Recovery Strategy

Do not restart with many figures. Recover the pipeline in strict dependency
order.

1. Recover condensate and density logic.
2. Recover scalar phase and Faraday rotation maps from the same condensate.
3. Recover one imaging mode, preferably PCI, because it exposes the phase-plate
   background and finite-NA blur clearly.
4. Recover DGI only after PCI arrays match.
5. Recover Faraday after scalar phase and PCI/DGI propagation are verified.
6. Recover deterministic camera binning before stochastic camera frames.
7. Recover multi-shot deterministic evolution before noisy filmstrips.

### First Recovery Milestone

Recover only one or two notebook-aligned figures:

1. Section 18.1 Stage 1 condensate figure:
   - 3D TF density line cuts;
   - across-cigar column density map.
2. Section 19.1 Step 1 phase-map figure:
   - scalar phase map;
   - centre lineouts along `y` and `z`.

These are the right first targets because all later figures depend on their
arrays.

### Required Numerical Checks

Before visual inspection, compare arrays numerically:

- `R` from notebook formulas vs helper state radii;
- `n_col` from notebook formulas vs helper state column density;
- `gax`, `GA`, `GB`, and `profile` shape and centre values;
- `column_density_map = n_col[0]*profile`;
- `phi_peak(1.5e9, n_col[0])`;
- `phase_map = phi_peak*profile`;
- lineouts `phase_map[Ngrid//2, :]` and `phase_map[:, Ngrid//2]`;
- min, max, centre value, integrated sum, and selected lineout points.

Suggested tolerances:

- exact shape match: `(1024, 1024)` for the main notebook imaging grid;
- relative tolerance around `1e-10` for deterministic formula arrays if using
  the same grid and constants;
- looser visual-only tolerance is not acceptable until the numerical arrays
  match.

### Required Manual Visual Checks

For the first two figures, inspect:

- axis labels and units: `um`, `cm^-3`, `cm^-2`, or `rad` as appropriate;
- correct imaging plane: across-cigar `y,z` for `axis=0`;
- zoom limits near notebook convention: `xlim(-45,45)`, `ylim(-12,12)` where
  used;
- colourbar units and limits;
- lineout direction and background/reference lines;
- whether the cloud's aspect ratio matches the notebook's cigar geometry.

## 14. Risks

- Treating the notebook as final calibrated theory would be wrong; it is a
  historical prototype with placeholder Faraday calibration.
- Treating helpers as complete figure recipes would also be wrong; they are
  tested computational kernels and thin orchestration layers.
- Changing profile exponent, imaging axis, grid size, FOV, display window, or
  colour limits can produce visually different but superficially plausible
  figures.
- Stochastic figures are especially fragile because the notebook RNG state is
  global and cell-order dependent.
- Faraday dual-port figures can be wrong if the ratio is formed before camera
  noise instead of after separately noising the two ports.
- Multi-shot frame figures can be wrong if the condensate object is not
  regenerated as `N0` and `R` evolve.

## 15. Bottom Line

The notebook pipeline is:

```text
physical constants and apparatus
-> Thomas-Fermi condensate state
-> projected column-density profile for a chosen imaging axis
-> scalar phase or Faraday rotation map
-> finite-NA scattered-field propagation
-> PCI, DGI, or Faraday contrast mechanism
-> deterministic camera binning
-> stochastic photon/read noise where needed
-> SNR/destructiveness calculation
-> multi-shot update of condensate number, temperature, phase, and SNR
```

Any future figure recovery should prove that this chain is matched numerically
from the condensate upward before attempting dissertation-quality plotting.
