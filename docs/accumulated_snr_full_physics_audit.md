# Full Physics and Code Audit of Accumulated SNR and N_max

## Verdict

**B. The current figure is physically valid only as an idealised clean-loss
analytical model.**

Its central far-detuned scaling argument is sound, but its `N_max` is not the
later notebook heating-plus-reabsorption frame budget, and
`SNR_shot*sqrt(N_max)` is not the correct quantitative accumulator once the
condensate and per-frame SNR evolve. The figure can remain as a clearly labelled
scaling illustration. It should not be cited as the realistic usable-frame or
accumulated-SNR prediction.

The structured audit outputs are under
`results/accumulated_snr_full_physics_audit/`.

## Parameter Provenance

The audit uses:

- `configs/notebook_v1_defaults.json` for atom, condensate, optical, camera,
  noise, and multishot parameters;
- `configs/dissertation_plots_v1.json` for the plotted operating point;
- `results/notebook_aligned_recovery/parameter_inventory.csv` and
  `unit_inventory.csv` to confirm provenance and units;
- the existing parameter, unit, and linear-approximation audits.

The audited operating point is 3.5 mW, 40 microseconds, imaging axis 0,
`QE=0.4`, read noise `7 e- rms`, `eta_coll=1.3`, and a nominal 30% condensate
depletion threshold. The saturation parameter is about 0.0259. No experimental
calibration is applied.

## Stage 1: Physics-First Audit

### Quantity definitions

`SNR_shot` is the expectation value of the selected per-image signal divided by
its one-frame noise standard deviation. In the current figure it is a
dimensionless peak-pixel analytical quantity, not a cloud-integrated or
resolution-element SNR.

For statistically independent frames with known frame-dependent signal and
variance, optimal RMS accumulation is

```text
SNR_total^2 = sum_i SNR_i^2.
```

The shortcut

```text
SNR_total = SNR_shot * sqrt(N_max)
```

requires all included frames to have the same SNR. It also assumes no
correlated probe noise, atom-number noise, drift, or camera correlations.

The full definitions, units, scope, and time dependence are recorded in
`physics_definition_table.csv`.

### What N_max means in the current figure

The exact implemented expression is

```text
N_max = -log(1-f) / (eta_coll*N_gamma),  f=0.30.
```

This follows the notebook clean-loss model

```text
N0(s) = N0(0)*exp(-eta_coll*N_gamma*s).
```

It is an exact continuous solution for the pulse index `s` under that model.
It is not inherently an integer number of camera frames. It is also not the
later realistic deep-trap model: the notebook explicitly calls clean loss an
optimistic bound and later models recoil energy as heating that depletes the
condensate fraction.

Three different counts must be kept separate:

1. continuous threshold pulse index;
2. integer pulses that do not exceed the budget, normally the floor;
3. notebook sequence array length, which includes frame 0 and records the first
   state at or beyond the threshold.

For example, at 1.5 GHz the clean-loss threshold is 29.58 pulses. The integer
budget not exceeding 30% is 29 pulses, the first threshold-crossing state has
index 30, and a sequence that stores states 0 through 30 has length 31. Calling
all three quantities `N_max` creates an off-by-one ambiguity.

The alternative formula

```text
log(1-f)/log(1-p_loss)
```

would apply if each pulse removed a fixed independent probability
`p_loss=eta*N_gamma`. It is close numerically here but is not the notebook's
chosen exponential model.

### Clean loss versus heating and reabsorption

The notebook heating model uses

```text
dE = N_gamma*(1+reabsorption)*E_rec
T_next = (T^4 + dE/A_E)^(1/4)
N0 = N_total*[1-(T/Tc)^3].
```

The 30% threshold therefore refers to condensate depletion, not necessarily a
30% reduction of total atom number. Heating can destroy condensed population
without ejecting the same fraction of atoms.

Benchmark continuous thresholds are:

| Detuning (GHz) | Clean loss | Heating + reabsorption |
| ---: | ---: | ---: |
| 0.75 | 7.40 | 3.20 |
| 1.50 | 29.58 | 13.76 |
| 3.00 | 118.32 | 56.23 |

Thus the clean-loss figure overstates the realistic frame budget by roughly a
factor of two at these parameters. This is a quantitative model mismatch, not
an algebra or unit error.

### Detuning scaling of N_max

At fixed power and pulse duration,

```text
N_gamma = (Gamma/2)*s/(1+s+delta^2)*tau.
```

Therefore `N_gamma ~ Delta^-2` and clean-loss `N_max ~ Delta^2` only in the
far-detuned limit. The exact denominator contains `1+s`, so the exponent is
asymptotic rather than mathematically exact.

Numerical fits give:

| Range (GHz) | Clean-loss exponent | Heating + reabsorption exponent |
| --- | ---: | ---: |
| 0.75-3.0 | 1.9998 | 2.0628 |
| 0.75-5.0 | 1.9998 | 2.0450 |
| 0.75-1.5 | 1.9996 | 2.1022 |
| 1.5-5.0 | 1.9999 | 2.0202 |

Clean-loss scaling is extremely close to quadratic over this scan. The
heating result grows slightly faster because the initial-density reabsorption
fraction decreases with detuning. In the current notebook sequence,
`N_gamma` is cloud-independent and reabsorption is evaluated from the initial
column densities, so later density evolution changes signal but not the
per-frame deposited-energy coefficient.

### Is N_max mode-independent?

Yes, under the controlled comparison in this figure. PCI and DGI place
different elements after the atom-light interaction. For identical incident
power, pulse duration, detuning, and atomic transition, scattering and heating
occur upstream of the mode-selecting optics. The same destruction budget
therefore applies to both modes.

This statement would not apply to notebook examples that deliberately use
different powers or pulse durations for PCI and DGI.

### Is sqrt(N_max) justified?

Only for the ideal identical-frame argument. In the actual multishot model:

- condensate number decreases;
- Thomas-Fermi radii and column density change;
- scalar phase decreases;
- PCI/DGI signals and SNR change from frame to frame;
- temperature increases in the heating model.

The correct existing repository convention is
`sqrt(sum_i SNR_i^2)`. Independent benchmark sequences show that replacing the
evolving values by the initial SNR overestimates accumulated PCI SNR by about
5.6%-10.6% across the audited points.

For clean loss, explicit PCI accumulated SNR increases from 36.58 at 0.75 GHz
to 38.66 at 3 GHz. For heating plus reabsorption it increases from 24.50 to
26.82. The ideal flatness therefore becomes an approximate trend with about
6%-10% variation, not an exact invariant.

The RMS formula still assumes independent noise. Technical probe fluctuations,
atom-number fluctuations, camera drift, and common-mode correlations are not
included and can prevent indefinite square-root accumulation.

### PCI and DGI signal and noise definitions

The figure's PCI observable is the notebook weak-phase tangent

```text
Delta I_PCI = 2*t_p*phi.
```

The bright carrier sets the photon-shot variance, and read noise is added in
quadrature. The full notebook quantitative PCI path instead uses the exact
transfer function, an axis-dependent blur estimate, and camera photon scale.
The full imaging helper propagates `exp(i*phi_map)` through the pupil.

The figure's DGI observable uses the ideal opaque-stop expression

```text
I_DGI = 4*sin^2(phi/2).
```

In that ideal limit the detected signal photons are also the photon-shot-noise
background. The full notebook DGI path has a finite reference amplitude
`10^(-OD/2)` with `OD=4`, plus pupil propagation and camera binning. At large
detuning the leakage background can matter and should enter the photon-shot
variance.

The figure compares peak-pixel analytical observables, not matched spatial ROIs.
It is therefore suitable for scaling, not for a definitive mode-efficiency
ranking.

The ideal PCI/DGI ratio near 2 is expected under the selected observables. PCI
uses the carrier as a homodyne reference and gives approximately
`2|phi|sqrt(N_ph)`, while ideal DGI gives approximately
`|phi|sqrt(N_ph)`. This is not a coding normalisation bug, but neither is it a
universal experimental factor.

### Noise model findings

- Quantum efficiency is included in the detected photon scale before the
  Poisson variance is formed.
- PCI photon noise uses the bright carrier background.
- Ideal DGI photon noise uses signal photons, which is correct only in the
  opaque-stop/no-leakage limit.
- Read noise is one `7 e- rms` contribution per camera output pixel per frame.
- The figure uses analytical variances, not stochastic Monte Carlo samples.
- Atom-number noise, probe technical noise, common-mode fluctuations, and drift
  are absent.
- Camera gain cancels because SNR is evaluated in electron-count units.

### Re-evaluated invariance

| Case | Audit result |
| --- | --- |
| PCI shot-noise limit | asymptotically constant in the clean-loss weak-phase model |
| PCI shot + read noise | approximately constant here because carrier noise and read-noise terms are detuning-independent |
| DGI shot-noise limit | asymptotically constant only in the opaque-stop, signal-shot-noise regime |
| DGI shot + read noise | decreases; high-detuning finite-range slope approaches `-1` |

The observed DGI high-detuning slope near `-0.872` is consistent with a finite
scan that has not reached a perfectly read-noise-dominated asymptote. Exact
`sin(phi/2)`, the transition between photon and read noise, and the finite
detuning range account for the difference from `-1`. In a full OD=4 model,
leakage-background photon noise would introduce an additional crossover.

## Stage 2: Code-to-Physics Traceability

The complete trace table is in `code_physics_traceability.csv`. Main findings:

- phase and scattering helper formulas exactly match the notebook conventions;
- the figure's `N_max` line correctly implements the stated continuous
  clean-loss expression;
- the multishot helper correctly implements frame-dependent clean-loss and
  heating bookkeeping;
- `accumulate_snr` correctly implements RMS accumulation;
- the ambiguity is semantic and model-level: fractional threshold, integer
  pulses, crossing index, and stored states are not the same quantity;
- the figure uses analytical SNR proxies rather than the full spatial pipeline.

No sign, logarithm, Hz/rad-s, QE, or basic scattering-factor error was found.
No critical physics bug was found in the frozen helper APIs.

## Issue Register

### Quantitative model mismatches

1. The figure uses clean-loss `N_max`; the realistic notebook model is heating
   plus reabsorption and predicts substantially fewer frames.
2. `SNR_initial*sqrt(N_max)` overestimates the explicit evolving-frame RMS
   accumulation by up to about 10.6% in the audited benchmarks.

### Implementation ambiguity

- `N_max` is fractional, while acquired frames are integer. The current
  sequence includes frame 0 and the first threshold-crossing state, so sequence
  length must not be reported directly as allowed pulses.

### Approximations requiring caveats

- quadratic detuning scaling is asymptotic;
- analytical peak-pixel SNR is not full Fourier/camera/ROI SNR;
- DGI stop leakage is omitted;
- independent noise excludes technical correlations and drift.

There are no issues classified as a critical physics error.

## Recommended Dissertation Definition

For a quantitative dissertation result:

1. Use the explicit heating-plus-reabsorption multishot sequence.
2. Define `N_max` as the integer number of acquired pulses whose post-pulse
   condensate depletion remains within the chosen threshold.
3. Keep the continuous threshold value separately as an analytical diagnostic.
4. Accumulate frame-dependent information with
   `sqrt(sum_i SNR_i^2)`, using a matched ROI or estimator for each mode.
5. Report clean loss as an optimistic comparison bound, not the main realistic
   prediction.
6. State whether the stopping criterion is condensate depletion, total atom
   loss, temperature, or another experimentally measurable quantity.

After experimental calibration, reabsorption, camera noise, correlated noise,
and the stopping threshold should be updated from laboratory data.

## Final Answers

- **Is current N_max algebraically correct?** Yes, for the continuous
  exponential clean-loss model.
- **What is currently implemented?** A fractional continuous pulse index for
  30% condensate loss under optimistic clean loss.
- **What should be used for the quantitative dissertation result?** An integer
  stopping budget from the explicit heating-plus-reabsorption condensate
  sequence.
- **Is `SNR_shot*sqrt(N_max)` justified?** Only as an identical independent-frame
  idealisation. Use `sqrt(sum_i SNR_i^2)` for the migrated sequence.
- **Is the PCI/DGI comparison fair?** Fair for the stated analytical peak-pixel
  scaling comparison, not for a final full-pipeline efficiency claim.
- **Does `N_max ~ Delta^2` survive the full model?** Approximately, but not
  exactly. Heating and detuning-dependent reabsorption change the exponent and
  coefficient.
- **Does accumulated-SNR invariance survive the full multishot model?** Only as
  an approximate asymptotic trend; explicit benchmarks vary by several to ten
  percent.
- **Can the current figure be used?** Yes, only as the explicitly labelled
  idealised clean-loss analytical scaling figure.
