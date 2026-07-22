# Reconstruction architecture

- **Status:** active inverse and observable contract
- **Active consumers:** reconstruction code, synthetic benchmarks, Chapters 6
  and 7, and the experimental-analysis plan
- **Update trigger:** a change to the raw measurement model, latent object,
  integration support, reported observable vector, uncertainty contract or
  production reconstruction entry point
- **Retirement rule:** replace only with one indexed successor that preserves
  the current observable and sealed-evidence boundaries

## Purpose

The reconstruction package uses raw camera channels to estimate low-order
condensate observables through the same finite-aperture and detector model used
to generate the synthetic measurements. Its primary scientific output is

```text
q = (A, y_c, z_c, w_major),
```

with uncertainty or a data-consistency interval. It is not a generic image
denoiser, a super-resolution method or an attempt to make a weak image look
sharp.

Resonant absorption imaging (RAI) remains the appropriate terminal reference
when a detailed density distribution is required. Dual-port Faraday imaging
(DPFI) is used for its lower disturbance and repeated access to spatially
resolved physical information. The two methods therefore have different
measurement objectives rather than competing on reconstructed-image clarity.

The baseline inverse does not assume that the cloud is Thomas-Fermi, dipolar,
fragmented or supersolid. Preparation metadata and pulse timing may be supplied
with an observation, but no morphology label is passed to the fit. A flexible
density field is retained internally as a nuisance representation required by
the nonlinear optical forward operator. It is not interpreted as a unique
recovered density image. Physical state models can be compared with the
estimated observables afterward.

## Information supplied to the inverse

The input contract separates preparation metadata from image data:

- a resonant-absorption observation contains atom, reference and dark frames,
  or a pre-normalised image with its calibration convention;
- a non-destructive observation contains an ordered pulse sequence and its raw
  camera channels;
- a dual-port exposure contains simultaneous `H` and `V` arrays;
- a dark-field exposure contains the crossed-analyser count array;
- timing, detuning, fluence and detector metadata remain attached to each
  exposure.

The current synthetic studies use one exposure at a time. A temporal or
destruction model can be added later without changing the single-frame inverse,
and must not force an abrupt physical change to remain smooth in time.

## Package boundaries

```text
src/non_destructive_image/reconstruction/
  observations.py       experimental input and timing contracts
  contracts.py          grid, pupil, ROI, detector and response contracts
  object_models.py      interchangeable column-density representations
  measurements.py       raw dual-port and dark-field forward operators
  noise.py              Poisson-Gaussian camera model
  initialise.py         smooth-reference initial estimates
  density_initialise.py blind shape-flexible initial estimates
  fit.py                 smooth-reference inverse
  density_fit.py         non-negative shape-flexible inverse
  regularisation.py      physically scaled curvature operator
  diagnostics.py         residual and simulation-only truth metrics
  credibility.py         local support and conditional uncertainty tools
  observables.py          physical moments on a fixed object-plane support
  synthetic_morphologies.py independent analytic stress inputs
  resolution.py          physical-camera reduced-grid construction
  studies/               frozen benchmark, provenance and reporting workflows
```

Production fitting is separated from truth-only assessment. Synthetic truth is
used to select and test a reconstruction contract, but is never supplied to an
individual fit.

## Object models

The smooth reference model fits a projected five-parameter Thomas-Fermi-like
profile. It is useful for testing the complete inverse chain and for estimating
low-order size and population quantities when that model is scientifically
appropriate. It cannot establish an unexpected feature that lies outside its
assumed family.

The active benchmark uses `NonnegativeBilinearDensityModel`. Column density is
represented by non-negative nodal coefficients on a coarse object-plane grid.
This representation lets morphology vary while `A`, centroid and major-axis
width are inferred. Its coefficients and interpolated map are nuisance
variables, not individual measured pixels. Knot spacing is part of the inverse
contract: modes finer than the aperture and camera support cannot become
measured information merely by increasing the parameter count.

## Optical and detector contract

Both readouts use the complete finite-aperture circular-polarisation forward
model. The Faraday response remains

\[
\theta_F(y,z)=\kappa_F C_\phi n_{\mathrm{col}}(y,z).
\]

Here `kappa_F=theta_F/phi` uses
`theta_F=(Phi_+-Phi_-)/2`. The approved isolated-transition estimate for the
fully spin-polarised axial `166Er` case is the signed value `-45/91` under the
declared helicity, propagation-direction and H/V-port conventions; an effective
apparatus response still requires calibration. The sealed v4 reconstruction
benchmarks were generated before this correction with `kappa_F=1`. They retain
method provenance but are not current Faraday-amplitude or SNR predictions.

The dual-port likelihood retains the raw `H` and `V` electron counts. The
normalised difference `S` is used for visualisation and initialisation, not as a
noise-free fitted datum. The dark-field likelihood retains the single
crossed-analyser count channel and therefore has the weaker quadratic response.

For each physical camera pixel,

\[
K_i\sim\operatorname{Poisson}(\mu_i),\qquad
Y_i=K_i+\epsilon_i,\qquad
\epsilon_i\sim\mathcal{N}(0,\sigma_r^2).
\]

The active Chapter 4/5 forward scenario uses the Hamamatsu ORCA-Fusion
C14440-20UP at `M=10`: a `6.5 um` sensor pixel represents `0.650 um` in the
object plane, the `100 um` field is sampled on a centred `153 x 153` analysis
crop, `QE=0.65`, effective `NA=0.130`, and the manufacturer-typical Ultra quiet
read noise is `0.7 e- rms` per physical pixel and readout. The sealed v4 inverse
continues to use its historical `NA=0.080`, `kappa_F=1` and `1.4 e- rms`
contract. It must not be described as quantitative evidence for the new active
scenario. A future reconstruction regeneration requires a new result family
and a fresh review of the basis and regularisation at the finer optical
resolution.

The reconstruction calculation uses a `306 x 306` object grid with the same
`100 um` field and exact physical camera pixels. It is a computational
reduction of the canonical `1024 x 1024` propagation, not a different detector.
Across the ten calibration and held-out morphology maps, the maximum signal-map
and peak discrepancies remain below the declared 1% gate.

## Inverse and regularisation

The shape-flexible fit uses frozen-weight iteratively reweighted least squares.
The weights are updated between outer iterations and held fixed inside each
least-squares solve, preventing the optimiser from reducing a residual merely
by increasing its predicted variance.

Non-negativity is an explicit physical constraint. Smoothness is supplied by a
separate thin-plate curvature operator approximating

\[
\int\!\left(n_{,yy}^2+2n_{,yz}^2+n_{,zz}^2\right)\,dA.
\]

Non-uniform divided differences and area weights make the penalty comparable
across basis spacings. Fixed zero ghost knots state the finite-support boundary
assumption explicitly. The data Jacobian is analysed before curvature rows are
added, so numerical stability supplied by the prior is not counted as measured
information.

Dual-port coefficients are initialised from the raw-count Jacobian near zero
density. That construction is not valid for dark-field data because the
quadratic response has zero first-order density information at the origin.
Dark-field initialisation instead uses moments of the square root of the
background-subtracted intensity to construct a starting envelope. The envelope
starts the nonlinear fit but does not constrain its final morphology.

## Calibration and held-out benchmark

The active benchmark contains five calibration and five disjoint held-out
analytic morphology classes. They cover smooth, asymmetric, modulated,
fragmented and notched structures with different positions, phases, widths and
peak numbers across the two splits. These maps are feature stress tests, not
stationary-state predictions for erbium.

All candidates see the same named raw camera draws. Candidate selection first
requires the declared fit-success and data-rank fractions, then minimises the
median pupil-supported relative error on the calibration split. Near-equivalent
candidates are resolved in favour of the lower-dimensional basis. The selected
contract is frozen before held-out observations are reconstructed.

For both readouts, Version 4 selects the 85-coefficient
`resolution_matched_17x5__curvature_30_um2` candidate. Independent range checks
at `100`, `300` and `1000 um2` show increasing calibration error, so `30 um2`
is bracketed against stronger over-smoothing rather than retained only because
it was the endpoint of the first scan.

## Physical observable extraction

The production observable layer operates after a single-frame inverse fit; it
does not change the fit, use a truth image or couple adjacent exposures. Every
map in a comparison is integrated on one immutable object-plane contract that
contains the `y` and `z` cell centres, physical cell areas and a support mask.
The formal outputs are the area-weighted integrated response `A`, the two-
dimensional centroid and the major-axis rms width from the largest eigenvalue
of the complete `2 x 2` covariance tensor. Minor width, aspect ratio and
principal-axis angle may remain diagnostics, but peak density and local
morphology are not primary outputs. Blank or threshold-unsupported inputs
retain their supported integral but return no fabricated centroid or width.

The active sealed held-out study uses `|y| <= 27 um`, `|z| <= 7.5 um` on the
`306 x 306` grid. This is evidence for the old reference-sized benchmark, not a
fixed future domain for wider Oxford or Er conditions. The next contract must
separate the camera fitting ROI, a wider latent-density support with guard
region, and the fixed physical support on which `A`, centroid and width are
defined. A finite observable support remains necessary because an arbitrarily
small remote tail can make a second moment arbitrarily large.

Relative integrated response and reconstructed depletion are defined only for
summaries with exactly identical coordinate grids, cell areas and masks. Their
physical interpretation additionally assumes a stable Faraday response,
detector gain and optical transfer across the sequence; the helper records but
cannot verify those assumptions. A common response scale cancels from the ratio,
but the integrated response is not an absolute atom number without an effective
calibration. A positive fitted map also does not establish cloud presence:
raw-channel improvement and blank evidence remain separate credibility
requirements.

The earlier `credibility.py` feature summary is retained for its frozen
bootstrap study. Its peak and axis-aligned rms features are not the production
parameter contract. The experimental path now requires an explicit
`ObservableIntegrationSupport`, calls `extract_density_observables()` for the
point estimate and passes exactly the same support to every bootstrap draw.
Support coordinates must exactly match the measurement/model object grid, and
the integration mask may be a strict subset of, but cannot extend beyond, the
latent model support. The exact support is retained in the summary and hashed
independently in the experimental result.

## Parameter credibility without a recovered image

Truth error is available only in simulation. The production-facing study
therefore reports four additional controls:

1. **Raw-channel residuals.** The fitted density is propagated back to `H`,
   `V` or dark-field counts and compared at the declared noise scale.
2. **Data/prior mode support.** Generalised local mode fractions and effective
   degrees of freedom separate likelihood support from curvature support.
3. **Conditional detector-noise spread.** Parametric Poisson/read-noise draws
   are refitted with the frozen operator, basis, support and prior.
4. **Prior and blank sensitivity.** Curvature variants and atom-free frames
   expose choice dependence and positive-density false reconstructions.

The conditional bootstrap now retains an aligned joint sample matrix in the
fixed order `(A, y_c, z_c, w_major)`. Unsupported moments remain in their draw
row as `NaN` plus an explicit support mask; the row is not deleted. A formal
marginal interval and joint covariance are reported only when every requested
refit converges and supports the relevant quantities. Otherwise the interval
is marked partial or unsupported, retains requested/successful/supported draw
counts, and has no bounds. This bootstrap still does not include response
calibration, optical-model error, support uncertainty or morphology ambiguity,
so its intervals are descriptive detector-noise spreads rather than calibrated
confidence statements on the physical condensate. Empirical coverage must be
checked on an independent synthetic ensemble. Detector/statistical,
inverse-contract and calibration contributions must remain distinguishable.

Non-negativity converts symmetric noise into a positive density bias. A fit can
therefore converge on a blank frame. Optimiser success and visual smoothness are
not signal-detection criteria. Numerical moment support and raw-channel cloud
evidence are reported separately: the extractor may return a finite centroid
or width, but that number is not physically reportable as a detected cloud
without adequate raw-channel/blank evidence. The current code does not invent
that evidence when no matched blank reference is supplied. If reasonable
latent shapes give materially different values while fitting the counts equally
well, the observable must receive a broad interval, a one-sided bound or an
unsupported status. A regulariser-selected map must not conceal that ambiguity.

## Reproducible entry points

The current commands and required local interpreter are maintained once in
`docs/reproducibility.md`.

The benchmark, credibility and observable-replay studies are long-running.
Each retained result records its config snapshots, source hashes, runtime
versions, Git branch, generation commit and dirty-state flag. Benchmark and
observable figures are generated only after the corresponding sealed numerical
artifact set passes its manifest checks.

## Current boundary

The maintained forward code and active configs implement the approved signed
erbium response `kappa_F=-45/91`. The sealed v4 inverse result family predates
that correction and remains frozen under its recorded `kappa_F=1` operator.
Regeneration under the active response, aperture and read-noise contract would
require a new result family and a fresh inverse-selection review.

The current code does not implement a mixed Poisson-Gaussian exact likelihood,
instrument-background nuisance fits, data-consistency/profile intervals, a full
posterior or a publication-scale coverage ensemble. Its bootstrap remains
conditional on a fixed forward operator, latent basis and support, observable
integration support and regulariser.

Before a formal initial-condition suite run, the user must approve the
condensate subset, camera fitting ROI, latent envelope, guard region and fixed
observable support. The preliminary seven-condition files remain unsealed and
must not be treated as a current result. Point bias and empirical interval
coverage are tested only after that contract is frozen.

The code also does not classify a latent map as a droplet, supersolid, soliton,
vortex or other phase, and it does not claim recovery of sub-resolution local
structure. Such statements require independent physical evidence beyond the
DPFI observables declared here.
