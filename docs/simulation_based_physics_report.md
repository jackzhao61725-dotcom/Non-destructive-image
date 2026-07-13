# Simulation-Based Physics Report

## Version, Purpose, and Scientific Boundary

This report consolidates the physical model, numerical conventions, validated
results, parameter provenance, approximation limits, and future development
path of the `Non-destructive-image` repository.

It is written to support an MSc dissertation, but it is not a substitute for a
literature review or an experimental calibration report. Its authority is
implementation traceability:

1. the historical notebook defines the Version 1 computational route;
2. exported notebook sections expose that route for audit;
3. the helper package preserves migrated equations and numerical conventions;
4. regression tests and stored outputs check numerical behaviour;
5. explicit configs define which parameter set belongs to each result.

The notebook is therefore a historical computational reference, not a claim of
final calibrated theory. The current simulation is representative and
uncalibrated. In particular, the Faraday coupling factor `kappa_F = 1.0` is a
phenomenological placeholder rather than a measured or microscopic Er result.

## 1. Research Question

The central question is not simply which imaging method produces the brightest
image. It is:

> How much statistically useful information about an ultracold condensate can
> be acquired for a fixed physical disturbance of the same atoms?

For continuous imaging, signal and damage are coupled. Increasing probe power
or exposure time raises the number of detected photons, but also increases
spontaneous scattering. Increasing detuning suppresses scattering faster than
it suppresses a dispersive phase shift, but weaker image contrast eventually
meets camera read noise, finite numerical aperture, technical noise, and a
finite experimental time budget. Repeated frames can be combined, but the
condensate itself evolves after each pulse.

The simulator therefore separates six physical layers:

```text
condensate state
  -> light-atom interaction
  -> coherent image formation
  -> camera detection and noise
  -> repeated-measurement evolution
  -> information-versus-destruction analysis
```

This separation is scientifically useful because each layer introduces a
different assumption and a different possible source of systematic error.

## 2. Model Hierarchy

### 2.1 Exact within the implemented Version 1 model

The following calculations use the full expressions implemented by the code,
not their weak-signal expansions:

- the two-level dispersive lineshape `delta/(1+delta^2)`;
- the scattering denominator `1+s+delta^2`;
- complex object fields `exp(i phi)`;
- FFT, pupil multiplication, and inverse FFT without additional shifts;
- PCI and DGI field recombination;
- circular Faraday fields `exp(+/- i theta_F)`;
- linear-basis recombination and dual-port intensity ratios;
- Poisson photon noise and Gaussian read noise;
- framewise RMS SNR accumulation `sqrt(sum_i SNR_i^2)`.

"Exact" here means exact relative to the chosen Version 1 equations. It does
not imply that those equations include all Er atomic structure or experimental
systematics.

### 2.2 Analytical approximations used for interpretation

The following are scaling arguments or idealised reference models:

- `phi ~ 1/Delta` and scattering `N_gamma ~ 1/Delta^2` at large detuning;
- PCI signal `Delta I ~ 2 t_p phi` at weak phase;
- ideal dark-field signal `I ~ phi^2` or `theta_F^2`;
- dual-port signal `S ~ 2 theta_F` at small rotation;
- identical-frame accumulation `SNR_total = SNR_shot sqrt(N_max)`;
- continuous clean loss as an optimistic destruction bound.

Quantitative image and evolving-sequence results retain the fuller numerical
expressions.

### 2.3 Physics not yet included

- multi-level Er scalar/vector/tensor polarisability;
- calibrated circular-transition strengths;
- vector-polarisability derivation of `kappa_F`;
- coherent multiple scattering;
- full radiation transport for reabsorbed photons;
- dipolar, LHY, droplet, supersolid, or mixture state dynamics;
- technical laser noise, drift, reference-image noise, and pixel covariance;
- experimental parameter posterior distributions.

## 3. Canonical Version 1 Parameter Set

The canonical notebook-aligned source is
`configs/notebook_v1_defaults.json`. All internal quantities use SI units.

### 3.1 Atomic and optical constants

| Quantity | Symbol | Value | Status |
| --- | --- | ---: | --- |
| Isotope | | `166Er` | Version 1 species |
| Atomic mass | `m` | `2.7564948562e-25 kg` | `166 u` convention |
| Transition wavelength | `lambda` | `401 nm` | notebook value |
| Natural linewidth | `Gamma` | `1.8535396656e8 rad s^-1` | `Gamma/2pi = 29.5 MHz` |
| Resonant cross-section | `sigma_0` | `7.6786733412e-14 m^2` | `3 lambda^2/(2 pi)` |
| Saturation intensity | `I_sat` | `597.9627 W m^-2` | `59.7963 mW cm^-2` |
| Photon energy | `E_ph` | `4.9537303140e-19 J` | `h c/lambda` |
| Recoil energy | `E_rec` | `4.9526328170e-30 J` | `(hbar k)^2/(2m)` |
| Recoil temperature | `T_rec` | `358.718 nK` | `E_rec/k_B` |
| Recoil velocity | `v_rec` | `5.9945 mm s^-1` | `hbar k/m` |

The cross-section and saturation intensity are effective two-level quantities.
They are appropriate for reproducing the notebook model, but a calibrated Er
analysis must revisit transition strengths, polarisation, Zeeman populations,
and multilevel corrections.

### 3.2 Condensate inputs

| Quantity | Value |
| --- | ---: |
| Condensate atom number `N0` | `2.5e4` |
| Scattering length `a_s` | `72 a0 = 3.8100759192e-9 m` |
| Trap frequencies `(f_x,f_y,f_z)` | `(293, 14, 233) Hz` |
| Initial temperature `T0` | `200 nK` |
| Self-consistent total atoms | `116166.78` |
| Self-consistent critical temperature | `216.826 nK` |
| Initial condensate fraction | approximately `0.215` |

The self-consistent total atom number is not an independently measured input.
It is the value required by the ideal trapped-gas condensate-fraction relation
to reconcile `N0=25000` with `T0=200 nK`.

### 3.3 Optical and camera geometry

| Quantity | Value |
| --- | ---: |
| Probe `1/e^2` diameter | `24 mm` |
| First/second focal length | `150/300 mm` |
| Magnification | `2` |
| Numerical aperture | `0.080` |
| Rayleigh scale at NA 0.080 | `3.0576 um` |
| Comparison Rayleigh scale at NA 0.40 | `0.6115 um` |
| Camera pixel | `2.93 um` at camera |
| Object-plane pixel | `1.465 um` |
| Simulation grid | `1024 x 1024` |
| Field of view | `100 um x 100 um` |
| Camera binning | `15 x 15` simulation samples |
| Binned frame | `68 x 68` pixels after truncation |
| Quantum efficiency | `0.40` |
| Read noise | `7 e- rms` per binned pixel |

Probe power and exposure are deliberately not universal defaults. They must be
stated for every camera or destruction result.

## 4. Thomas-Fermi Condensate Physics

### 4.1 Trapped interacting condensate

For harmonic angular frequencies `omega_i = 2 pi f_i`, define

```math
\bar{\omega} = (\omega_x\omega_y\omega_z)^{1/3},
\qquad
a_{\mathrm{ho}} = \sqrt{\frac{\hbar}{m\bar{\omega}}}.
```

The Version 1 Thomas-Fermi chemical potential is

```math
\mu = \frac{\hbar\bar{\omega}}{2}
\left(\frac{15N_0a_s}{a_{\mathrm{ho}}}\right)^{2/5}.
```

The peak density and Thomas-Fermi radii are

```math
n_0 = \frac{\mu m}{4\pi\hbar^2a_s},
\qquad
R_i = \sqrt{\frac{2\mu}{m\omega_i^2}}.
```

The three-dimensional profile is

```math
n(x,y,z)=n_0\max\left(0,
1-\frac{x^2}{R_x^2}-\frac{y^2}{R_y^2}-\frac{z^2}{R_z^2}
\right).
```

The numerical atom-number check is

```math
N_0 = \frac{8\pi}{15}n_0R_xR_yR_z.
```

### 4.2 Derived canonical state

| Derived quantity | Value |
| --- | ---: |
| `bar(omega)` | `618.914 rad s^-1` |
| `bar(omega)/2pi` | `98.503 Hz` |
| Harmonic oscillator length | `0.7862 um` |
| Chemical potential | `6.5681763e-31 J` |
| `mu/k_B` | `47.573 nK` |
| Peak density | `3.4002134e20 m^-3` |
| Radii `(R_x,R_y,R_z)` | `(1.1858, 24.8171, 1.4912) um` |
| Atom-number reconstruction | `25000` |

The cloud is therefore cigar-shaped along `y`. The small transverse radii are
below the 3.06 um Rayleigh scale of the existing imaging arm, while the long
axis is resolved.

### 4.3 Column densities and imaging axes

Integrating along principal axis `i` gives the peak column density

```math
\tilde n_i = \frac{4}{3}n_0R_i.
```

The full x-integrated distribution in the displayed `y,z` plane is

```math
n_{\mathrm{col}}(y,z)=\tilde n_x
\max\left(0,1-\frac{y^2}{R_y^2}-\frac{z^2}{R_z^2}\right)^{3/2}.
```

Equivalent expressions apply for integration along `y` and `z`. The notation
`n_col(y,z)` denotes a two-dimensional distribution; `tilde n_x` denotes only
its peak scalar.

| Line of sight | Display plane | Peak column density |
| --- | --- | ---: |
| `x` | `y,z` | `5.3759625e14 m^-2` |
| `y` | `x,z` | `1.1251121e16 m^-2` |
| `z` | `x,y` | `6.7603305e14 m^-2` |

The long-axis line of sight gives much larger phase, but its image plane
contains two unresolved transverse radii. It is more naturally an integrated
atom-number probe than a structural image.

### 4.4 Scaling during depletion

Within the Thomas-Fermi model,

```math
\mu\propto N_0^{2/5},\quad
R_i\propto N_0^{1/5},\quad
n_0\propto N_0^{2/5},\quad
\tilde n_i\propto N_0^{3/5}.
```

Thus repeated imaging changes both signal amplitude and cloud size. A
multishot calculation that holds the first image fixed cannot reproduce the
full evolving sequence.

## 5. Dispersive Light-Atom Interaction

### 5.1 Detuning convention

The code accepts ordinary-frequency detuning magnitude `Delta_Hz` and linewidth
`Gamma` in radians per second. The dimensionless detuning is

```math
\delta = \frac{2(2\pi\Delta_{\mathrm{Hz}})}{\Gamma}.
```

Thesis tables should label detuning as `|Delta|/2pi` in GHz. At
`|Delta|/2pi = 1.5 GHz`, `delta = 101.694915`, normally reported as
`delta approximately 101.7`.

### 5.2 Scalar phase

For peak column density `tilde n`, the Version 1 scalar phase is

```math
\phi = \frac{\sigma_0\tilde n}{2}
\frac{\delta}{1+\delta^2}.
```

The spatial phase map is `phi_peak` multiplied by the projected Thomas-Fermi
profile. In the far-detuned limit,

```math
\phi \sim \frac{\sigma_0\tilde n}{2\delta}\propto\Delta^{-1}.
```

At 1.5 GHz:

| Imaging axis | Peak phase | Interpretation |
| --- | ---: | --- |
| `x`, across cigar | `0.20294 rad` | weak finite phase; suitable for exact PCI calculation |
| `y`, along cigar | `4.247 rad` | phase wrapped; a single-valued linear PCI inversion is invalid |
| `z`, across cigar | approximately `0.255 rad` | weak finite phase |

At 13 GHz along `y`, the peak phase is approximately `0.490 rad`, close to the
notebook's practical linear-PCI boundary.

### 5.3 Residual optical depth

The residual optical depth is

```math
\mathrm{OD}(\Delta)=\frac{\sigma_0\tilde n}{1+\delta^2},
```

which scales as `Delta^-2`. Absorptive attenuation therefore falls faster than
the dispersive phase. This is the basic reason off-resonant imaging can be less
destructive, but "less destructive" is not equivalent to "non-destructive":
spontaneous scattering remains finite.

### 5.4 Probe intensity and saturation

For a Gaussian beam of `1/e^2` diameter `D`, the Version 1 atoms-at-centre
intensity is

```math
I_{\mathrm{peak}} = \frac{2P}{\pi(D/2)^2}.
```

The factor of two relative to the area-averaged intensity is intentional. The
cloud is much smaller than the 24 mm probe diameter. The saturation parameter
is `s=I_peak/I_sat`.

### 5.5 Scattering

The scattering rate and photons scattered per atom per pulse are

```math
R_{\mathrm{sc}}=\frac{\Gamma}{2}
\frac{s}{1+s+\delta^2},
\qquad
N_\gamma=R_{\mathrm{sc}}\tau.
```

For the canonical continuous-imaging operating point
`P=3.5 mW`, `tau=40 us`, `|Delta|/2pi=1.5 GHz`, the code gives

```text
s = 0.0258769
N_gamma = 0.00927474 photons per atom per pulse.
```

In the far-detuned, low-saturation regime,

```math
N_\gamma\propto\frac{P\tau}{\Delta^2}.
```

At 1.5 GHz, replacing exact phase and scattering denominators by their
far-detuned limits produces relative errors around `1e-4`; the asymptotic
language is accurate for trends, though the code retains exact denominators.

## 6. Coherent Image Formation

### 6.1 Projected object field

For scalar imaging,

```math
E_{\mathrm{obj}}(x,y)=\exp[i\phi(x,y)],
\qquad
E_{\mathrm{sc}}=E_{\mathrm{obj}}-1.
```

Subtracting one separates the unscattered carrier from the spatially varying
field. This separation is central to PCI and DGI: both propagate the same
scattered field but treat the carrier differently.

### 6.2 Fourier convention and finite aperture

The numerical pupil is

```math
P(f_x,f_y)=
\begin{cases}
1, & \sqrt{f_x^2+f_y^2}\le\mathrm{NA}/\lambda,\\
0, & \mathrm{otherwise}.
\end{cases}
```

The propagated scattered field is

```math
E_{\mathrm{sc}}' = \mathrm{IFFT2}\{\mathrm{FFT2}(E_{\mathrm{sc}})P\}.
```

No extra FFT shift, padding, or alternative normalisation is introduced. The
same convention is used by the notebook and helper package.

### 6.3 Phase-contrast imaging (PCI)

The phase plate changes the carrier to

```math
E_{\mathrm{ref}}^{\mathrm{PCI}}=t_p e^{i\theta},
\qquad t_p=0.95,\quad\theta=\pi/2.
```

The image is

```math
I_{\mathrm{PCI}}=\left|t_pe^{i\theta}+E_{\mathrm{sc}}'\right|^2.
```

The no-atom reference is `t_p^2=0.9025`, not unity. For weak phase and ideal
resolution, the contrast is approximately linear:

```math
\Delta I_{\mathrm{PCI}}\approx2t_p\phi.
```

PCI is therefore a homodyne-like measurement: the bright carrier acts as a
local reference for the weak scattered field. The exact field expression is
used for images.

At the canonical 1.5 GHz across-cigar phase, the recovered ideal PCI image has
centre intensity `1.15857`, mean `0.90437`, and shape `1024 x 1024`.

### 6.4 Dark-ground imaging (DGI)

The Version 1 finite stop is represented by residual carrier amplitude

```math
E_{\mathrm{ref}}^{\mathrm{DGI}}=10^{-\mathrm{OD}_{\mathrm{stop}}/2},
\qquad \mathrm{OD}_{\mathrm{stop}}=4.
```

Thus the reference amplitude is `0.01` and its intensity is `1e-4`. The image is

```math
I_{\mathrm{DGI}}=\left|10^{-\mathrm{OD}_{\mathrm{stop}}/2}
+E_{\mathrm{sc}}'\right|^2.
```

For an ideal opaque stop and a spatially uniform point,

```math
I_{\mathrm{DGI}}=|e^{i\phi}-1|^2=4\sin^2(\phi/2)
\approx\phi^2.
```

DGI rejects the bright carrier, giving a dark background, but the signal is
quadratic at weak phase and can become read-noise limited. The finite OD=4
leakage is included in full Fourier calculations but not in the idealised
analytical Fig. 3.2 DGI scaling curve.

The canonical recovered DGI image centre is `0.0159565`.

### 6.5 Phenomenological Faraday rotation

Version 1 defines

```math
\theta_F(\mathbf r)=\kappa_F\phi(\mathbf r),
\qquad \kappa_F=1.0.
```

This is a phenomenological mapping. It does not derive separate circular
transition phases from Clebsch-Gordan coefficients or Er vector
polarisability. Consequently, Faraday amplitudes are not calibrated physical
predictions.

The circular object fields are

```math
E_+=e^{+i\theta_F},\qquad E_-=e^{-i\theta_F}.
```

Each scattered part is propagated through the same pupil, after which the full
circular fields are reconstructed and recombined:

```math
E_x=\frac{E_+'+E_-'}{2},
\qquad
E_y=\frac{i(E_+'-E_-')}{2}.
```

The common scalar phase is neglected in this pure-rotation model.

### 6.6 Dark-field Faraday

A crossed analyser transmits the rotated component:

```math
I_{\mathrm{dark}}=|E_y|^2.
```

For an ideal uniform small rotation this behaves as
`sin^2(theta_F) approximately theta_F^2`. The exact propagated field is used by
the simulation. At the canonical point, the ideal dark-field centre intensity
is `0.0159564`.

### 6.7 Dual-port Faraday

The two analyser ports are

```math
I_u=\frac{|E_x+E_y|^2}{2},
\qquad
I_v=\frac{|E_x-E_y|^2}{2}.
```

The normalised difference is

```math
S=\frac{I_v-I_u}{I_v+I_u}.
```

For an ideal uniform field, `S=sin(2 theta_F) approximately 2 theta_F`. The two
bright ports use nearly all detected light and reject common multiplicative
intensity fluctuations in the ratio, although unequal gain, offset,
polarisation leakage, and uncorrelated camera noise still matter.

At the canonical ideal-image point, the dual-port signal centre is `0.251169`.

## 7. Finite-Signal Validity

At the canonical across-cigar operating point,

```text
max |phi| = max |theta_F| = 0.202942 rad.
```

The simulation uses exact fields, but the errors of common explanatory
approximations are:

| Approximation at peak | Relative error |
| --- | ---: |
| `exp(i phi)` versus `1+i phi` field | `2.057%` |
| `sin(theta_F)` versus `theta_F` | `0.690%` |
| `sin^2(theta_F)` versus `theta_F^2` | `1.384%` |
| `sin(2 theta_F)` versus `2 theta_F` | `2.799%` |

Therefore:

- weak-phase language is suitable for explaining scaling;
- the canonical point is weak but not infinitesimal;
- quantitative PCI, DGI, and Faraday results should retain exact fields;
- along-cigar imaging at 1.5 GHz is phase wrapped and must not use a linear
  inversion.

## 8. Camera and Noise Model

### 8.1 Photon scale

The object-plane pixel is `p_obj=p_cam/M`. For normalised intensity `I/I0=1`,
the expected detected photon/electron count is

```math
N_{\mathrm{ph,pix}}=
\frac{I_{\mathrm{peak}}p_{\mathrm{obj}}^2\tau\,\mathrm{QE}}
{hc/\lambda}.
```

The result depends explicitly on probe power, exposure, quantum efficiency,
and pixel area. It must never be inferred from an unlabeled default.

### 8.2 Binning

The `1024 x 1024` high-resolution image is truncated to `1020 x 1020` and
averaged in `15 x 15` blocks, producing `68 x 68` camera pixels. Averaging
represents the ideal intensity ratio over a camera pixel; expected photon
counts are then obtained from the photon scale.

### 8.3 Stochastic detection

For binned normalised image `I_b`, the code generates

```math
C=\mathrm{Poisson}(\max(I_b,0)N_{\mathrm{ph,pix}})
+\mathcal N(0,\sigma_{\mathrm{read}}),
```

and optionally returns `C/N_ph,pix`. The RNG is passed explicitly as a
`numpy.random.Generator`; helpers do not create hidden global seeds.

The notebook used `default_rng(7)` globally. Recovery scripts can exactly
replay an isolated, documented call order with seed 7, but cannot claim exact
reproduction of an arbitrary interactively executed notebook RNG state.

### 8.4 Distinct camera parameter contexts

The notebook-aligned deterministic PCI camera recovery uses `P=2 mW` and
`tau=100 us`, giving `1532.324` photons per unit-intensity camera pixel. The
Faraday camera panel uses `P=5 mW` and `tau=100 us`, giving `3830.809` photons
per pixel. The continuous-imaging route uses `P=3.5 mW` and `tau=40 us`, giving
`1072.627` photons per pixel.

These are separate examples. Their outputs must not be compared without
renormalising the photon budget and stating the observable.

## 9. Signal-to-Noise Definitions

There is no single universal SNR in this repository. A valid SNR statement must
specify the estimator and spatial support.

### 9.1 Analytical peak-pixel SNR

The idealised Fig. 3.2 PCI model uses

```math
S_{\mathrm{PCI}}=2t_p|\phi|N_{\mathrm{ph}},
\qquad
\sigma_{\mathrm{PCI}}=
\sqrt{t_p^2N_{\mathrm{ph}}+\sigma_{\mathrm{read}}^2}.
```

The ideal DGI model uses

```math
I_{\mathrm{DGI}}=4\sin^2(\phi/2),
\quad
S_{\mathrm{DGI}}=I_{\mathrm{DGI}}N_{\mathrm{ph}},
\quad
\sigma_{\mathrm{DGI}}=
\sqrt{S_{\mathrm{DGI}}+\sigma_{\mathrm{read}}^2}.
```

These are analytical object-centre pixel quantities before NA/PSF blur and
without spatial summation.

### 9.2 Resolution-element SNR

A resolution-element calculation sums signal over a defined pixel block and
adds photon and read variances over the same pixels. It is not obtained by
arbitrarily multiplying a peak-pixel SNR. The notebook has a specific
`3 x 3`-pixel block example, but it is not the same observable as the later
fixed ROI.

### 9.3 Matched-ROI SNR

The full evolving analysis uses a fixed 228-pixel spatial ROI, shared by PCI
and DGI. For pixel `j`,

```math
s_j=(I_{\mathrm{atoms},j}-I_{\mathrm{ref},j})N_{\mathrm{ph}},
```

with diagonal variance

```math
v_j=I_{\mathrm{atoms},j}N_{\mathrm{ph}}
+\sigma_{\mathrm{read}}^2.
```

The matched-template statistic is

```math
\mathrm{SNR}^2=\sum_{j\in\mathrm{ROI}}\frac{s_j^2}{v_j}.
```

This estimator includes distributed image information and provides a fairer
same-ROI mode comparison than one peak pixel. It still assumes independent
pixels and a noiseless known reference image.

## 10. Destructiveness Models

### 10.1 Continuous optimistic clean loss

The clean-loss model assumes

```math
N_0(q)=N_0(0)e^{-\eta N_\gamma q},
```

where `q` is a continuous pulse index and `eta=1.3` represents effective
collisional secondaries. For tolerated condensate loss `f`,

```math
N_{\max}^{\mathrm{clean}}=
\frac{-\ln(1-f)}{\eta N_\gamma}.
```

This is a fractional continuous threshold, not automatically an integer number
of acquired frames. It is retained as an optimistic analytical upper bound.

### 10.2 Heating and condensate melting

Because `E_rec/k_B approximately 359 nK` is below a microkelvin-scale trap
depth, a recoiling atom may remain trapped and thermalise. The Version 1
heating model uses

```math
A_E=3\frac{\zeta(4)}{\zeta(3)}\frac{k_B}{T_c^3},
```

```math
\Delta E=N_\gamma(1+r)E_{\mathrm{rec}},
\qquad
T_{q+1}=\left(T_q^4+\frac{\Delta E}{A_E}\right)^{1/4},
```

and

```math
N_0(q)=N_{\mathrm{tot}}\left[1-\left(\frac{T_q}{T_c}\right)^3\right].
```

The threshold is condensate depletion, not necessarily total atom loss.

### 10.3 Reabsorption

The current approximation is

```math
r=\frac{1}{3}\sum_{i=x,y,z}
\left[1-e^{-\mathrm{OD}_i(\Delta)}\right].
```

At the canonical 1.5 GHz point, `r=0.0297087`. This is a three-axis average,
not a full angular radiation-transport calculation. The full evolving script
holds this initial-density value fixed during a sequence, so density-updated
reabsorption remains future work.

### 10.4 Frame-count semantics

Three quantities must remain distinct:

1. continuous threshold, such as `29.582` clean-loss pulses;
2. strict integer accepted frames, such as `29` clean or `13` heating-aware;
3. a stored sequence length that may include frame 0 and/or a threshold-crossing
   state.

The quantitative full-model convention accepts a frame only when its post-pulse
state remains within 30% condensate depletion. The first crossing state is
excluded.

## 11. Accumulated SNR and the Invariance Argument

For statistically independent frames with framewise SNR `SNR_i`, the repository
uses

```math
\mathrm{SNR}_{\mathrm{acc}}=
\sqrt{\sum_i\mathrm{SNR}_i^2}.
```

If every frame is identical, this reduces to

```math
\mathrm{SNR}_{\mathrm{acc}}=
\mathrm{SNR}_{\mathrm{shot}}\sqrt{N}.
```

In the far-detuned photon-limited PCI scaling model,

```math
\mathrm{SNR}_{\mathrm{shot}}
\propto\frac{\sqrt{P\tau}}{|\Delta|},
\qquad
N_{\max}\propto\frac{\Delta^2}{P\tau}.
```

Therefore the factors cancel at fixed destruction budget. Ideal DGI also gives
an approximately constant accumulated shot-noise SNR because its per-frame
signal scales quadratically while its photon noise scales with the square root
of that signal.

This cancellation is a controlled scaling result, not a universal theorem. It
breaks when:

- read noise dominates weak DGI frames;
- finite stop leakage contributes background photons;
- the cloud depletes and shrinks;
- heating changes the accepted frame count;
- integer stopping creates steps;
- the estimator or ROI changes;
- technical noise is correlated between frames.

## 12. Verified Accumulated-SNR Results

### 12.1 Fig. 3.2 idealised scaling set

Parameter set:

```text
|Delta|/2pi scan: 0.75 to 5.0 GHz
P: 3.5 mW
tau: 40 us
axis: x, across cigar
QE: 0.40
read noise: 7 e- rms for shot+read curves
observable: analytical peak object-space camera pixel
image model: pre-NA/PSF, no spatial sum
N_max: continuous 30% optimistic clean loss
accumulation: identical frames
```

At 1.5 GHz:

| Mode | Noise | SNR/frame | Continuous `N_max` | Accumulated SNR |
| --- | --- | ---: | ---: | ---: |
| PCI | shot only | `13.2931` | `29.5820` | `72.3001` |
| PCI | shot + read | `12.9689` | `29.5820` | `70.5370` |
| DGI | shot only | `6.63514` | `29.5820` | `36.0881` |
| DGI | shot + read | `4.56457` | `29.5820` | `24.8264` |

Across 0.75-5 GHz, DGI shot-only accumulated SNR spans
`35.8975-36.1460`, a relative range of `0.689%`. It is therefore
approximately detuning-independent, not exactly constant.

The PCI/DGI shot-noise ratio at 1.5 GHz is `2.0034`. This is consistent with
the chosen homodyne-like PCI and ideal dark-field observables. It is not a
universal factor of two for arbitrary apertures, stops, ROIs, or noise models.

### 12.2 Full evolving matched-ROI set

The optical point remains `P=3.5 mW`, `tau=40 us`, axis `x`, `QE=0.4`, and
read noise `7 e-`, but the model now includes finite aperture, OD=4 DGI
leakage, camera binning, a fixed 228-pixel ROI, heating, reabsorption, state
regeneration, strict integer stopping, and framewise RMS accumulation.

At 1.5 GHz, 13 frames are accepted:

| Mode | Noise | Initial SNR | Final SNR | Full accumulated SNR |
| --- | --- | ---: | ---: | ---: |
| PCI | shot only | `37.5020` | `28.9887` | `120.3309` |
| PCI | shot + read | `36.7210` | `28.3657` | `117.7894` |
| DGI | shot only | `19.7160` | `14.9902` | `62.8219` |
| DGI | shot + read | `7.43766` | `4.66287` | `21.9222` |

The identical-initial-frame approximation overestimates these totals by about
12-22% at this operating point. Across the full scan, the overestimate ranges
from 8.2% to 27.7%, depending on mode and noise.

The `117.79` PCI value is not a replacement or correction factor for the
`70.54` Fig. 3.2 value. They answer different estimator questions.

### 12.3 Frame budgets

At `|Delta|/2pi=1.5 GHz`, 30% condensate-depletion threshold:

| Parameter set | Clean continuous | Heating continuous | Heating+reabs continuous | Strict realistic frames |
| --- | ---: | ---: | ---: | ---: |
| `P=2 mW`, `tau=15 us` | `138.049` | `66.109` | `64.202` | `64` |
| `P=2 mW`, `tau=40 us` | `51.768` | `24.791` | `24.076` | `24` |
| `P=3.5 mW`, `tau=40 us` | `29.582` | `14.166` | `13.758` | `13` |

This table explains why the old phrase "52/25/24 at 15 us" was wrong: those
values were computed with 40 us defaults.

### 12.4 Along-cigar reference

For axis `y`, `P=2 mW`, `tau=40 us`:

- peak phase at 1.5 GHz: `4.247 rad`, phase wrapped;
- peak phase at 13 GHz: `0.490 rad`;
- clean-loss continuous budget at 13 GHz: approximately `3888` pulses.

The large frame count is an optimistic clean-loss quantity, not a calibrated
camera or heating-aware prediction.

## 13. Corrected Legacy Numerical Claims

### 13.1 Deprecated `171.7`

The old accumulated SNR `171.7` combined:

- single-frame SNR evaluated with a 100 us camera default;
- `N_max` evaluated with the 40 us global pulse duration.

Using the same legacy ideal model consistently at 40 us gives `108.60`. Using
the canonical Fig. 3.2 definitions gives `72.30` shot only or `70.54` with read
noise. The old `171.7` must not be quoted as a valid 40 us result.

### 13.2 Corrected explicit-40-us operating-point SNRs

| `|Delta|/2pi` | Axis | P | Ideal QE=1 peak pixel | QE=0.4/read7 peak pixel with NA/PSF | QE=0.4/read7 resolution element |
| ---: | --- | ---: | ---: | ---: | ---: |
| `1.5 GHz` | `x` | `2 mW` | `15.094` | `5.781` | `11.692` |
| `13 GHz` | `y` | `2 mW` | `36.453` | `6.015` | `8.525` |

All values use `tau=40 us`. The older rows
`23.86/9.32/18.88` and `57.63/9.70/13.78` used the 100 us default and are not
valid 40 us results.

## 14. Faraday Optimisation Layer

The analysis layer evaluates deterministic scalar proxies:

- `faraday_signal_scale = |theta_F|`;
- scattered photons per atom;
- approximate reabsorption;
- a destructiveness metric `N_gamma(1+r)`;
- `signal_per_scattered_photon`;
- `signal_to_destruction`;
- optional `2 |theta_F| sqrt(N_ph)` SNR proxy.

One-dimensional detuning, power, and exposure sweeps call the same operating-
point function. They do not propagate Faraday images, average stochastic
noise, or optimise multiple variables jointly.

The current `configs/dissertation_results_v1.json` is a workflow placeholder,
not the canonical notebook cloud. It uses, among other values:

```text
N = 100000
peak column density = 2.0e14 m^-2
I_sat = 560 W m^-2
photons per camera pixel = 500
P = 2 mW
tau = 40 us
kappa_F = 1.0
```

These differ from the canonical `N0=25000`, x-axis peak column density
`5.376e14 m^-2`, and `I_sat=597.963 W m^-2`. Therefore the generated Faraday
optimisation outputs demonstrate software workflow and qualitative trade-offs;
they must not be mixed with canonical notebook-aligned numbers or presented as
final operating-point predictions.

The present metric `|theta_F|/N_gamma` tends to favour lower power/exposure and
larger detuning because it does not yet include all camera detectability and
finite acquisition constraints. A future objective should include a calibrated
measurement likelihood or information metric.

## 15. Absorption/RAI Calibration Readiness

The calibration module currently provides deterministic preprocessing:

```math
\mathrm{OD}=-\ln\left(
\frac{I_{\mathrm{atoms}}-I_{\mathrm{dark}}}
{I_{\mathrm{probe}}-I_{\mathrm{dark}}}
\right),
```

with clipping to avoid invalid logarithms. It can extract peak OD, integrated
OD, centre coordinates, and RMS widths.

These observables can eventually constrain:

- atom number or density scale;
- Thomas-Fermi radii or empirical widths;
- magnification and object-plane pixel size;
- camera offset, gain, and read noise;
- probe intensity and detuning offsets;
- optical-depth and cross-section scale;
- residual absorption and reabsorption corrections.

No file-format-specific laboratory loader, fit, uncertainty model, held-out
validation, or `kappa_F` calibration is currently implemented.

## 16. Numerical Provenance Register

Every generated conclusion should emit the following columns:

```text
quantity | value | |Delta|/2pi | P | tau | imaging axis |
normalisation | N_max model | QE/read | repository path
```

The most important current result contexts are:

| Quantity family | `|Delta|/2pi` | P | tau | Axis | Normalisation | `N_max` model | QE/read | Source |
| --- | ---: | ---: | ---: | --- | --- | --- | --- | --- |
| Condensate state | n/a | n/a | n/a | all | absolute SI density | n/a | n/a | `results/notebook_aligned_recovery/condensate_stage/` |
| Scalar phase map | `1.5 GHz` | n/a | n/a | `x` | absolute rad | n/a | n/a | `results/notebook_aligned_recovery/phase_stage/` |
| Ideal PCI/DGI images | `1.5 GHz` | n/a | n/a | `x` | `I/I0`, 1024 grid | n/a | n/a | `results/notebook_aligned_recovery/{pci,dgi}_stage/` |
| Ideal Faraday fields | `1.5 GHz` | n/a | n/a | `x` | dark `I/I0`; dual-port ratio | n/a | n/a | `results/notebook_aligned_recovery/faraday_stage/` |
| PCI camera example | `1.5 GHz` | `2 mW` | `100 us` | `x` | binned `I/I0` | n/a | `0.4/7 e-` | `results/notebook_aligned_recovery/camera_stage/` |
| Faraday camera panel | `1.5 GHz` | `5 mW` | `100 us` | `x` | camera `I/I0` and ratio | n/a | `0.4/7 e-` | `results/notebook_aligned_recovery/faraday_camera_panel/` |
| Fig. 3.2 scaling | scan | `3.5 mW` | `40 us` | `x` | analytical peak pixel | clean continuous 30% | `0.4/0 or 7 e-` | `results/dissertation_plots_v1/accumulated_snr_invariance/` |
| Full multishot | scan | `3.5 mW` | `40 us` | `x` | 228-pixel matched ROI | heating+reabs, strict frames | `0.4/0 or 7 e-` | `results/dissertation_plots_v1/full_multishot_accumulated_snr/` |
| Placeholder Faraday optimisation | configured sweep | configured | configured | scalar | signal/scattering proxy | no sequence | placeholder | `results/faraday_optimisation_v1/` |

If any entry is unknown, the number is not ready for a thesis table.

## 17. Reproducibility and Validation

Install the environment and set the local import path:

```powershell
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
$env:PYTHONPATH="src;."
$env:PYTHONUTF8="1"
```

Core validation:

```powershell
pytest -q
python scripts\validate_notebook_sections.py
```

At the 2026-07-13 closure review, the full suite completed with `127 passed`,
notebook-section validation passed, both accumulated-SNR generators passed,
the numerical consistency audit passed with the documented legacy
corrections, and the full-multishot SVG reproduced byte-for-byte across two
successive generations.

Numerical consistency and the two accumulated-SNR result families:

```powershell
python scripts\generate_accumulated_snr_invariance_plot.py --config configs\dissertation_plots_v1.json
python scripts\generate_full_multishot_accumulated_snr.py --config configs\dissertation_plots_v1.json
python scripts\audit_thesis_numerical_consistency.py --config configs\thesis_numerical_contract_v1.json
```

All approved recovery and result workflows:

```powershell
python scripts\run_all_dissertation_figures.py
```

Primary machine-readable sources:

- `configs/notebook_v1_defaults.json`;
- `configs/dissertation_plots_v1.json`;
- `configs/thesis_numerical_contract_v1.json`;
- `regression/baseline/notebook_outputs.json`;
- `regression/baseline/imaging/pci_dgi_imaging_baseline_v1.npz`;
- `regression/baseline/imaging/faraday_imaging_baseline_v1.npz`.

## 18. Scientific Risks and Error Budget

### 18.1 Atomic-physics risk: high

`kappa_F` is not calibrated, and the two-level model omits real Er multilevel
structure. This dominates the uncertainty of absolute Faraday predictions.

### 18.2 Experimental parameter risk: high

Atom number, widths, effective cross-section, probe intensity at the cloud,
magnification, QE, and read noise are notebook values or placeholders. Until
RAI/absorption data are fitted, absolute SNR and frame counts are representative.

### 18.3 Destruction-model risk: medium to high

Clean loss and heating bracket different physical pictures. The current
heating model is more defensible for a deep trap, but uses a simplified
equilibrium condensate relation and fixed initial-density reabsorption.

### 18.4 Imaging-model risk: medium

The coherent pupil model captures finite NA but omits aberrations, defocus,
phase-plate size/misalignment, DGI stop geometry, polariser extinction, camera
MTF, and pixel cross-talk.

### 18.5 Statistical-model risk: medium

Poisson and independent Gaussian read noise are implemented. Technical noise,
reference subtraction, common-mode correlations, drift, and atom-number
fluctuations are not.

### 18.6 Numerical risk: controlled for Version 1

Helper outputs are regression-tested against notebook-derived expressions and
array baselines. The larger remaining risk is selecting the wrong config,
normalisation, or estimator, not basic FFT or algebra drift. The thesis
parameter contract directly addresses this risk.

## 19. Dissertation-Safe Claims

The following claims are supported by the current simulation if their scope is
stated:

1. Far-detuned dispersive phase falls approximately as `1/Delta`, while
   scattering and residual absorption fall approximately as `1/Delta^2`.
2. PCI converts weak phase to a carrier-referenced linear intensity change,
   whereas ideal DGI produces a quadratic dark-field signal.
3. In an ideal photon-limited clean-loss model, accumulated SNR can become
   approximately independent of detuning, power, and pulse duration at fixed
   destruction budget.
4. Read noise breaks this cancellation most strongly for dark-field signals.
5. Heating, reabsorption, condensate depletion, finite aperture, integer
   stopping, and evolving images remove exact invariance.
6. A matched spatial estimator can produce different absolute SNR values from
   a peak-pixel estimator even at the same optical operating point.
7. Across-cigar imaging is weak phase but partially unresolved; along-cigar
   imaging has much larger phase and can be phase wrapped.
8. Dual-port Faraday is linear in small rotation and has common-mode rejection
   in principle, but the absolute rotation scale is uncalibrated.

Claims that are not yet supported:

- a final optimal detuning, power, or exposure for the experiment;
- an experimentally accurate number of non-destructive frames;
- an absolute Faraday rotation for `166Er` based on `kappa_F=1`;
- a universal factor-of-two PCI advantage;
- quantitative performance on droplets, supersolids, or mixtures;
- experimental superiority of one imaging mode without a common calibrated
  estimator and noise model.

## 20. Recommended Development Path

### Phase 1: close the numerical contract

- require explicit `tau`, power, axis, QE, read noise, and estimator arguments
  in every thesis generator;
- emit a synchronized parameter register for every numerical conclusion;
- keep canonical notebook-aligned and placeholder optimisation configs visibly
  separate;
- remove or quarantine deprecated mixed-default notebook tables from thesis use.

### Phase 2: calibrate the apparatus

- add laboratory image loaders with immutable raw-data provenance;
- fit dark offset, gain, QE or effective photon scale, read noise, magnification,
  and PSF;
- fit atom number, column-density scale, and cloud widths from absorption/RAI
  images;
- validate on held-out runs and report uncertainties.

### Phase 3: calibrate Faraday physics

- derive scalar, vector, and tensor polarisabilities for the relevant Er state,
  transition, detuning, and polarisation;
- map state populations and magnetic-field geometry to circular phase shifts;
- fit or validate `kappa_F` against simultaneous/paired destructive references;
- propagate calibration uncertainty into all Faraday objectives.

### Phase 4: improve repeated-measurement realism

- update reabsorption from each evolving density;
- include atom loss as well as heating where experimentally indicated;
- add probe/reference technical noise and temporal correlations;
- model polariser extinction, PBS imbalance, and camera-channel gain mismatch;
- validate framewise depletion and SNR against experiment.

### Phase 5: use an information-theoretic objective

- define the estimated parameter: atom number, centre, width, mode amplitude,
  or density perturbation;
- construct a calibrated likelihood for the camera data;
- compute Fisher information or estimator covariance per frame;
- optimise information per scattered photon or per allowed condensate loss;
- only then perform two- or three-dimensional operating-point optimisation.

### Phase 6: extend the state model

- add phenomenological droplet and supersolid profiles as additive providers;
- accept external GPE or experimentally reconstructed density maps;
- later add component-resolved mixture states and transition-specific optical
  responses;
- preserve the Thomas-Fermi Version 1 path as a regression reference.

## 21. Final Assessment

The repository now provides a coherent and well-tested Version 1 simulation
chain from a Thomas-Fermi `166Er` condensate to PCI, DGI, Faraday, camera, and
multishot observables. Its strongest dissertation contribution is not a claim
of final experimental prediction. It is the explicit separation of:

- dispersive signal generation;
- spontaneous-scattering destructiveness;
- mode-specific coherent image formation;
- finite-aperture and camera detection;
- state evolution under repeated probing;
- estimator-dependent information accumulation.

The corrected numerical contract resolves the principal internal consistency
risk: results generated with different exposure defaults, spatial
normalisations, or destruction models must no longer be presented as if they
were directly comparable. Within the documented Version 1 assumptions, the
physics and code are mutually consistent and reproducible. The next decisive
scientific gain will come from experimental calibration and a microscopic or
measured Faraday coupling, not from adding more uncalibrated plots.
