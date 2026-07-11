# Accumulated-SNR Invariance and Its Breakdown

## Purpose

This document develops the accumulated-SNR argument from first principles and
separates four levels that must not be conflated:

1. far-detuned single-frame scaling;
2. idealised clean-loss accumulation with identical frames;
3. an evolving heating-plus-reabsorption multishot sequence;
4. future experimentally calibrated performance.

The current dissertation figure belongs to level 2. It is a useful analytical
scaling figure, not a realistic usable-frame prediction.

## Part I: Reconstructing the Physical Argument

### 1. General measurement model

Let one frame produce an observable `x` with a no-atom or reference expectation
`x0`. Define the useful signal

```text
A(Delta) = E[x|atoms] - E[x|reference]
```

and total one-frame variance

```text
sigma^2(Delta) = sigma_ph^2 + sigma_read^2 + sigma_tech^2 + ... .
```

Then

```text
SNR_shot(Delta) = |A(Delta)| / sigma(Delta).
```

`A` is not universal. It may be a field quadrature, intensity contrast, one
pixel, a summed region of interest, a lineout amplitude, a matched-filter
output, or an inferred physical parameter. Different observables can share the
same detuning exponent but have different numerical prefactors.

For independent frames with frame-dependent signals and variances, the optimal
RMS combination is

```text
SNR_total^2 = sum_i SNR_i^2.
```

Only when every included frame has the same SNR does this reduce to

```text
SNR_total = SNR_shot*sqrt(N).
```

The destruction per frame is likewise model-dependent. The current light-atom
formula gives scattered photons per atom per pulse,

```text
N_gamma = (Gamma/2) * s/(1+s+delta^2) * tau,
delta = 2*Delta/Gamma,
```

with the repository's explicit Hz/rad-s conversion. `N_gamma` is not itself
atom loss, condensate loss, or temperature increase. A destruction model maps
it to those quantities.

### 2. Linear-mode scaling

For a far-detuned scalar phase object,

```text
phi(Delta) ~ C_phi/Delta.
```

For PCI near quadrature, the weak-phase intensity contrast is

```text
Delta I_PCI ~ 2*t_p*phi.
```

If `N_c` carrier photons are detected in the selected pixel or ROI, the useful
electron signal is proportional to `phi*N_c`.

#### Photon-shot-noise limit

The carrier supplies the local oscillator and also its shot noise:

```text
sigma_ph ~ sqrt(N_c),
SNR_PCI,shot ~ 2*|phi|*sqrt(N_c) ~ Delta^-1
```

at fixed probe power, exposure, optical throughput, QE, pixel, and ROI.

#### Shot plus read noise

With one read-noise contribution `sigma_r` per output pixel,

```text
SNR_PCI ~ 2*t_p*|phi|*N_c / sqrt(t_p^2*N_c + sigma_r^2).
```

At fixed `N_c`, the denominator is essentially detuning-independent, so the
same `Delta^-1` exponent remains. Read noise lowers the prefactor without
necessarily changing the detuning exponent. This is why PCI can remain nearly
flat after multiplication by `sqrt(N_max)` even when read noise is included.

#### Technical-noise limit

If fractional carrier noise dominates, `sigma_tech ~ epsilon*N_c`, then

```text
SNR_PCI,tech ~ |phi|/epsilon ~ Delta^-1.
```

The exponent can still look favourable, but repeated-frame accumulation may
not follow `sqrt(N)` because laser-intensity drift is often correlated. Thus a
single-frame exponent alone does not establish accumulated invariance.

Which photon number enters the denominator depends on the estimator. A raw PCI
pixel is dominated by carrier photons; a background-subtracted ROI includes
the shot variance of all contributing pixels; a balanced or matched estimator
can reject common-mode technical noise but cannot remove fundamental photon
shot noise.

### 3. Quadratic-mode scaling

For an ideal DGI opaque stop,

```text
I_DGI = |exp(i*phi)-1|^2 = 4*sin^2(phi/2) ~ phi^2 ~ Delta^-2.
```

#### Photon-shot-noise limit

If the dark port has no leakage background, the signal photons are also the
total detected photons. Therefore

```text
sigma_ph ~ sqrt(I_DGI*N_ph),
SNR_DGI,shot ~ sqrt(I_DGI*N_ph) ~ |phi|*sqrt(N_ph) ~ Delta^-1.
```

A quadratic intensity signal can therefore have a linear-in-phase SNR because
Poisson noise is the square root of the signal photon number.

#### Read-noise-dominated limit

When `I_DGI*N_ph << sigma_r^2`,

```text
SNR_DGI,read ~ I_DGI*N_ph/sigma_r ~ phi^2 ~ Delta^-2.
```

This steeper single-frame decline is the origin of the falling DGI accumulated
curve.

#### Mixed regime

The analytical interpolation is

```text
SNR_DGI = N_sig/sqrt(N_sig + sigma_r^2),
N_sig = 4*sin^2(phi/2)*N_ph.
```

Its local log slope moves continuously from `-1` toward `-2`. The audited
high-detuning slope is not exactly `-2` because the finite scan still crosses
the transition between photon and read noise.

A real DGI stop has finite leakage. Then the photon variance contains signal
plus leakage-background photons. At sufficiently weak phase, this background
can make even the nominal shot-noise-limited DGI SNR scale quadratically. The
current figure intentionally uses the ideal opaque-stop limit and therefore
does not show that crossover.

### 4. Clean-loss frame budget

In the far-detuned limit at fixed power and exposure,

```text
N_gamma ~ Delta^-2.
```

The notebook clean-loss survival model is

```text
N0(s)/N0(0) = exp(-eta_coll*N_gamma*s).
```

For a permitted fractional condensate loss `f`, the continuous threshold is

```text
N_max,cont = -log(1-f)/(eta_coll*N_gamma) ~ Delta^2.
```

This is exact for the chosen exponential survival law, but only asymptotically
quadratic because the exact scattering denominator is `1+s+delta^2`.

An additive approximation would use `N0/N0_initial ~ 1-eta*N_gamma*s` and
`N_max ~ f/(eta*N_gamma)`. A discrete fixed-probability model instead gives
`log(1-f)/log(1-p)`. These are close only for small per-pulse loss and are not
interchangeable definitions.

The continuous threshold, the integer number of allowed pulses, the index of
the first threshold-crossing state, and the stored sequence length are four
different quantities. The current figure uses the continuous threshold.

### 5. Ideal accumulated invariance

Combining `sqrt(N_max) ~ |Delta|` with the per-frame exponents gives:

| Mode and noise regime | Per-frame SNR | Accumulated scaling | Status |
| --- | --- | --- | --- |
| PCI, carrier shot noise | `Delta^-1` | `Delta^0` | asymptotically invariant |
| PCI, fixed read-noise term | `Delta^-1` | `Delta^0` | approximately invariant at fixed photons |
| PCI, correlated technical noise | often `Delta^-1` per frame | not generally `sqrt(N)` | no accumulated invariance guarantee |
| ideal DGI, signal shot noise | `Delta^-1` | `Delta^0` | asymptotically invariant |
| DGI, read-noise dominated | `Delta^-2` | `Delta^-1` | decreases with detuning |
| DGI with leakage-background shot noise | approaches `Delta^-2` at weak signal | approaches `Delta^-1` | ideal invariance eventually breaks |

No statement in this table is an exact all-detuning theorem. Each result is
conditional on a specified observable, photon budget, background, destruction
law, and noise covariance.

## Part II: Breakdown Mechanisms

### 6. Frame-to-frame evolution

The identical-frame shortcut fails because the probe changes the object being
measured. Clean loss lowers `N0`; heating lowers the condensate fraction even
without equivalent total-atom loss; Thomas-Fermi chemical potential, peak
density, radii, and column density change; phase and image contrast decline;
and the per-frame SNR falls.

The correct independent-frame accumulator is consequently

```text
SNR_total = sqrt(sum_i SNR_i^2).
```

Using the initial SNR for every frame overestimates the audited PCI sequence by
about 5.6%-10.6%. This bias would become larger for deeper depletion or a signal
with stronger density dependence.

Equal weighting is optimal only for equal frame variances and equal response.
For an estimated physical parameter, the correct weights should follow each
frame's derivative with respect to that parameter and its covariance.

### 7. Heating and reabsorption

Direct clean loss removes condensed atoms according to an imposed exponential
survival model. Recoil heating instead raises `T`; the total atom number can
remain nearly fixed while the condensate population falls through
`N0=N_total[1-(T/Tc)^3]`. Reabsorbed photons multiply the deposited recoil
energy and therefore accelerate condensate depletion.

Both clean loss and heating can retain an apparent exponent near `+2` because
both are ultimately driven by `N_gamma ~ Delta^-2`. They nevertheless have
different coefficients and stopping meanings. In the benchmark, heating plus
reabsorption permits roughly half as many pulses as clean loss. Similar
exponents do not imply equivalent physical budgets.

The current notebook evaluates the reabsorption factor from the initial cloud
and holds it fixed through the sequence. A future self-consistent model could
update it as density and geometry evolve, producing an additional departure
from simple scaling.

### 8. Finite aperture and estimator choice

Finite numerical aperture removes spatial frequencies and can reduce the
central contrast differently for narrow and extended clouds. Camera binning
changes both signal integration and the number of read-noise events. Pixel size
selects how much of the point-spread function enters one sample. A peak pixel,
central lineout, summed ROI, and matched filter are different estimators.

Effects that mainly change prefactors:

- constant throughput and QE;
- detuning-independent blur over a fixed cloud;
- fixed ROI efficiency;
- fixed phase-plate transmission.

Effects that can change the apparent exponent or mode ranking:

- density-dependent cloud size and blur;
- DGI leakage becoming dominant at high detuning;
- ROI size changing with cloud size;
- multiple read-noise contributions from a larger ROI;
- nonlinear PCI transfer at larger phase;
- estimator saturation or phase wrapping.

Thus peak-pixel PCI/DGI prefactors should not be interpreted as total
information efficiency.

### 9. Noise hierarchy

The useful hierarchy is:

1. **Photon shot noise only:** analytical Poisson limit.
2. **Shot plus read noise:** analytical detector variance used in the figure.
3. **Camera-level stochastic simulation:** explicit Poisson and Gaussian draws,
   binning, and normalisation.
4. **Technical probe noise:** intensity, pointing, frequency, and polarization
   fluctuations.
5. **Common-mode and temporal correlations:** drift that does not average as
   `sqrt(N)`.
6. **Atom-number and state-preparation noise:** fluctuations of the object.
7. **Model and calibration uncertainty:** uncertain density, magnification,
   cross section, read noise, gain, and destruction mapping.

Calling level 2 `realistic` is too strong. It is more precise to say `shot +
read noise`. Realistic performance requires at least levels 3-6 and experimental
calibration.

### 10. Mode-dependent ideal prefactors

Equal detuning exponents do not imply equal SNR. PCI mixes the scattered field
with a carrier local oscillator. In the current peak-pixel normalization this
gives approximately `2|phi|sqrt(N_ph)`. Ideal DGI counts dark-field photons and
gives approximately `|phi|sqrt(N_ph)`. The ratio near two is therefore expected
for these observables.

It is not fundamental. A different phase-plate transmission, finite stop,
aperture, ROI, or normalization changes the prefactor. A matched filter can
recover information spread over several pixels. A Fisher-information
comparison may assign equal or different efficiencies depending on which
photons, nuisance parameters, and detection losses are included. The current
factor of two is an estimator-dependent homodyne advantage, not a universal law
that PCI contains twice as much information as DGI.

## Part III: Thesis Claim Audit

The full claim table is stored in `claim_audit.csv`. The most important
corrections are:

- `accumulated SNR is independent of detuning` is too strong without
  `far-detuned`, `photon-limited`, `fixed clean-loss budget`, and
  `identical-frame` qualifiers;
- `the frames-versus-SNR curve is the same for every detuning` is contradicted
  at the several-percent level by the evolving multishot benchmark;
- `usable frames N_max (realistic)` is misleading where the plotted quantity is
  a continuous clean-loss threshold;
- `N_max*tau is fixed` is asymptotic under fixed power and destruction model,
  not an exact full-model invariant;
- `read noise breaks invariance` is mode- and observable-dependent: it strongly
  breaks an ideal quadratic dark-field mode, while for bright-background PCI it
  can mainly reduce the prefactor;
- `linear modes have a fundamental advantage` is unjustified. The demonstrated
  advantage is detector- and estimator-dependent.

The current figure metadata uses the qualified wording and is acceptable. The
stronger legacy notebook prose should not be copied verbatim into the thesis.

### Is accumulated SNR the best metric?

`SNR_shot*sqrt(N_max)` is useful for a compact scaling argument but is not the
best final performance metric. Better choices depend on the research question:

- **Main text:** explicit `sqrt(sum_i SNR_i^2)` from the evolving sequence,
  alongside condensate survival or temperature.
- **Mode comparison:** matched-filter or integrated-ROI SNR with identical
  optical and camera assumptions.
- **Parameter estimation:** total Fisher information or predicted uncertainty
  in phase, atom number, position, width, or another target parameter.
- **Destructiveness efficiency:** information per scattered photon and
  information per condensate atom depleted.
- **Appendix:** peak-pixel analytical scaling and alternative clean-loss bounds.

Fisher information is the most general theoretical metric because it includes
the derivative of the full image likelihood with respect to a specified
parameter. It is also more demanding: it requires a calibrated likelihood,
noise covariance, nuisance parameters, and an explicit scientific observable.

## Part IV: Layered Interpretation

### Level 1: Far-detuned single-frame scaling

Modelled: fixed cloud, fixed probe photons, weak phase, analytical noise.

Valid conclusion: phase falls as `Delta^-1`; scattering falls as `Delta^-2`;
linear and ideal dark-field shot-noise SNR fall approximately as `Delta^-1`.

Invalid conclusion: exact all-detuning equality or experimentally calibrated
mode performance.

Recommended wording: `In the far-detuned weak-phase limit, the leading
single-frame scalings are...`.

### Level 2: Clean-loss identical-frame accumulation

Modelled: continuous clean-loss threshold, identical independent frames,
analytical peak-pixel SNR.

Valid conclusion: accumulated photon-limited SNR is asymptotically independent
of detuning; quadratic read-noise-limited SNR decreases approximately as
`Delta^-1`.

Invalid conclusion: realistic usable frame number or exact multishot
invariance.

Recommended wording: `This idealised clean-loss construction exposes the
scaling cancellation between per-frame sensitivity and frame budget.`

### Level 3: Evolving heating-plus-reabsorption sequence

Modelled: temperature and condensate fraction evolve; phase and per-frame SNR
are recomputed; RMS SNR is accumulated.

Valid conclusion: the clean-loss cancellation remains a useful organizing
principle, but the accumulated result varies by several to ten percent and the
usable-frame coefficient is substantially lower.

Invalid conclusion: `SNR_initial*sqrt(N_max)` is the sequence result.

Recommended wording: `The full deterministic sequence retains only an
approximate remnant of the ideal scaling because the object and measurement
sensitivity evolve after every pulse.`

### Level 4: Experimentally calibrated performance

Modelled in future: measured cloud parameters, optical throughput, gain,
read-noise covariance, technical noise, drift, destruction response, and
possibly calibrated DGI leakage.

Valid conclusion now: none about final optimal detuning or absolute performance.

Recommended wording: `The present outputs are Version 1 predictions to be
regenerated after closed-loop calibration.`

## Publication-Ready Paragraphs

### Theory section

> In the far-detuned weak-phase limit, the dispersive phase decreases as
> `|Delta|^-1`, whereas the scattering probability per pulse decreases as
> `Delta^-2`. Under an idealised clean-loss budget this allows a pulse number
> proportional to `Delta^2`. Consequently, a photon-shot-noise-limited
> single-frame SNR proportional to `|Delta|^-1` yields an accumulated SNR that
> is asymptotically independent of detuning. This cancellation is a scaling
> argument rather than an exact property of the evolving experiment.

### Figure caption

> Accumulated SNR in the small-phase analytical model at a fixed 30% clean-loss
> budget. Dashed curves show photon-shot-noise limits and solid curves add
> camera read noise. The ideal PCI and DGI curves are approximately invariant
> with detuning, while read noise causes the quadratic DGI signal to lose this
> invariance. The clean-loss frame budget is an optimistic upper bound and the
> curves do not include frame-dependent condensate depletion, full spatial
> image formation, or technical noise.

### Results discussion

> The ideal cancellation does not remain exact in the deterministic multishot
> simulation. Recoil heating and reabsorption reduce the 30% condensate-depletion
> frame budget by approximately a factor of two relative to clean loss at the
> audited operating points. Recomputing the phase and PCI sensitivity after
> each pulse gives an accumulated RMS SNR 5.6%-10.6% below the identical-frame
> estimate. The scaling argument therefore remains qualitatively informative,
> but quantitative results must use the evolving sequence.

### Limitations section

> The analytical comparison includes photon shot noise and Gaussian camera read
> noise but not correlated probe noise, atom-number fluctuations, DGI stop
> leakage, camera drift, or experimental calibration. Peak-pixel SNR is also not
> equivalent to total image information. These effects may modify absolute
> prefactors, mode ranking, and the degree to which repeated measurements
> accumulate as the square root of frame number.

### Conclusion and future work

> A calibrated optimization should replace the continuous clean-loss frame
> count with an explicit heating-and-reabsorption stopping rule and should
> accumulate frame-dependent information from the full camera likelihood. A
> matched-ROI SNR is the immediate next step; Fisher-information optimization
> becomes appropriate once experimental noise and calibration parameters are
> available.

## Recommended Figure Strategy

Recommendation: **C, retain both figures with a clear separation.**

1. Keep the current idealised figure in the theory or scaling subsection. Its
   purpose is to expose the cancellation and the read-noise breakdown.
2. Add a full evolving multishot comparison in the results section using
   heating plus reabsorption and `sqrt(sum_i SNR_i^2)`.
3. If space is limited, move the ideal scaling figure to the appendix and keep
   the full sequence in the main text.

A two-panel combined figure is possible, but separate figures better prevent
the analytical upper bound from being mistaken for the quantitative result.

## Prioritised Issue Register

### Must correct before thesis submission

- Replace unconditional `independent of detuning` wording with asymptotic,
  model-qualified wording.
- Do not call fractional clean-loss `N_max` a realistic usable-frame count.
- Do not use `SNR_initial*sqrt(N_max)` as the quantitative multishot result.
- State that the 30% heating threshold is condensate depletion, not necessarily
  total atom loss.

### Should clarify

- Label analytical noise as `shot + read noise`, not `realistic`.
- Explain the estimator-dependent PCI/DGI factor of two.
- State that DGI uses the ideal opaque-stop limit.
- Distinguish peak-pixel SNR from integrated information.

### Optional refinement

- Add matched-ROI or matched-filter comparisons.
- Show clean loss and heating as separate bounds.
- Report continuous threshold and integer frame count separately.

### Future model work

- Update reabsorption with evolving density.
- Include OD=4 leakage in DGI noise.
- Add technical-noise covariance and drift.
- Move to parameter-specific Fisher information after calibration.

## Final Statements

1. **Strongest justified statement:** in the far-detuned weak-phase limit, at
   fixed probe photons and under an idealised clean-loss budget, photon-limited
   accumulated SNR is asymptotically detuning-independent; read noise breaks
   this cancellation for a quadratic dark-field observable.
2. **Strongest unjustified statement:** the realistic evolving experiment has
   an exactly detuning-independent accumulated SNR, a universally mode-equal
   shot-noise limit, or a final optimal detuning determined by this figure.
3. **Fate of the current figure:** retain it as a clearly labelled theory/scaling
   figure, paired with a future full multishot result; move it to the appendix
   if only one figure can remain in the main text.
4. **Sentences that must change:** unconditional invariance, `same for every
   detuning`, `realistic usable frames`, exact `N_max proportional to Delta^2`,
   and claims of a fundamental linear-mode advantage.
5. **Next numerical model:** matched-ROI PCI/DGI camera signals evaluated
   frame-by-frame with heating, reabsorption, evolving density, DGI leakage, and
   RMS accumulation.
6. **Suspected code error:** none found. The formulas implement their stated
   models; the problem is model scope and interpretation.
7. **Primary issue category:** modelling scope and estimator choice, followed by
   thesis wording; not a frozen-helper physics bug.
