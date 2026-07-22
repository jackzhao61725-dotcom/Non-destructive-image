# Multiframe Heating Model Optimisation Basis

- **Status:** approved implementation basis; not yet implemented
- **Applies to:** Chapter 5 repeated-imaging thermodynamics and the active
  \(^{166}\mathrm{Er}\) DPFI sequence model
- **Established:** 22 July 2026
- **Planned consumers:** one shared multiframe state engine, Figures 5.2 and
  5.4, their regression tests, and the associated dissertation discussion
- **Update trigger:** new Oxford thermodynamic evidence, a verified change to
  the atomic scattering response, or experimental calibration of heating,
  loss, trap depth or inter-frame equilibration
- **Retirement rule:** replace only after one verified successor reproduces
  the direct Oxford initial states, documents its evidence boundary and is used
  by every active Chapter 5 sequence figure

This document defines the smallest defensible replacement for the current
Version 1 multiframe heating model. It is an implementation specification, not
dissertation prose and not an apparatus-level heating calibration.

Until this specification passes Section 12, the stored Figures 5.2 and 5.4 and
their `N_dep`/`N_use` values remain frozen Version 1 screening outputs. The
implemented parameter contract is recorded separately in
`simulation_reference_parameters.md`.

## 1. Decision

The active replacement will be a **quasi-equilibrium, recoil-limited,
fixed-atom-number screening model**. It will estimate the evolution of

\[
(N_0,N_{\rm th},T)
\]

after repeated probe pulses. It will use the measured Oxford thermal
non-saturation behaviour and a bounded first-order energy model. It will not
attempt to reproduce every loss, collision or non-equilibrium process in the
apparatus.

The intended claim is conditional:

> Given the direct Oxford initial state, the measured non-saturation trend,
> complete equilibration between frames, fixed trapped atom number and
> recoil-only energy deposition, the model predicts a sensitivity span for
> condensate depletion across repeated exposures.

The output is an ideal screening estimate. It is not a calibrated prediction of
the number of frames obtainable in an Oxford experiment.

## 2. Why the Version 1 model must be replaced

The Version 1 model infers a total atom number from a specified \(N_0\) and
temperature using

\[
N_0=N_{\rm tot}\left[1-\left(T/T_c\right)^3\right],
\]

then reuses the same ideal saturated-gas relation after every pulse. This is
incompatible with the Oxford state used for the bimodal sample. The currently
selected raw state is

\[
N_0=18292.219843,\qquad
N_{\rm th}=199133.909886,\qquad
T=245.102432\ {\rm nK},
\]

so its directly measured-scale total is \(217426.129729\). Applying the old
ideal relation to this total and temperature predicts approximately \(50413\)
condensed atoms, which is \(2.76\) times the reported \(N_0\).

The mismatch is physical. In an interacting trapped gas the thermal component
continues to grow after condensation. The Oxford analysis measures this
non-saturation and extrapolates to \(N_0\rightarrow0\) to determine the critical
point. Its hold-time trajectory is close to isothermal and includes controlled
evaporative cooling and atom loss. It cannot be reinterpreted as a DPFI
frame-by-frame heating trajectory.

The existing \(N_{\rm dep}\) and \(N_{\rm use}\) values therefore remain frozen
Version 1 screening outputs until the model in this document passes its
acceptance tests.

## 3. Direct initial-state contract

The central trajectory retains the raw 300 ms state already used by the active
Oxford bimodal sample. Two further raw repetitions at the same dipole
orientation, trap, scattering length and hold time define the initial-state
envelope.

| Repetition | \(N_0\) | \(N_{\rm th}\) | \(T\) (nK) |
|---:|---:|---:|---:|
| 1, central | 18292.219843 | 199133.909886 | 245.102432 |
| 2 | 19081.354879 | 198665.044669 | 245.222440 |
| 3 | 18726.235710 | 184464.379726 | 239.046438 |

Each row is one state sample. The three values in a row must remain associated;
they must not be sampled independently. Their spread is an empirical
three-repetition envelope, not a complete measurement uncertainty or a
confidence interval. The accepted 337.5 ms binned state is reserved for a
closure check and will not replace or be mixed with the 300 ms initial
condition.

The current state is for

- \(^{166}\mathrm{Er}\);
- \(a_s=72a_0\);
- dipoles at \(\theta_d=0^\circ\);
- trap frequencies \(2\pi\times(294.0,14.0,233.3)\ {\rm Hz}\), subject to the
  exact measured-scale convention stored with the source data;
- \(\epsilon_{dd}\simeq0.902<1\).

Changing any of these conditions requires a new non-saturation contract.

## 4. Measured-scale non-saturation closure

Define the nominal ideal critical number

\[
\widetilde N_c^\circ(T)
=
\zeta(3)
\left(\frac{k_BT}{\hbar\widetilde{\bar\omega}}\right)^3.
\]

The published 0-degree fit has \(C=1.1961283098\) and
\(S=0.3275456539\). The Oxford combined calibration factor is
\(\beta=0.9202760715\). The implementation will collapse them into two
measured-scale constants:

\[
C_m=\beta C=1.10076826195,
\]

\[
S_m=\beta^{3/5}S=0.31161813155.
\]

The three raw triplets and trap frequencies are controlled by the Trap I,
\(\epsilon_{dd}=0.902\), 0-degree measurement series in the Oxford dataset. The
fit coefficients are controlled by `Fig2c_fits.csv`; the combined scale is
controlled by `calibrations.csv`. When the data are added to the configuration,
their dataset DOI, ZIP member names and source-file SHA-256 values must be
stored with the imported records. The runtime value of
\(\widetilde{\bar\omega}\) is the geometric mean of the three frequencies in
the selected record; it must not be read from the global pure-condensate
surrogate.

This prevents the runtime model from exposing the unresolved atom-number and
trap-frequency scale factors separately. The measured-scale closure is

\[
F_{\rm ns}(T,N_0)
=
\widetilde N_c^\circ(T)
\left[
C_m+S_m
\left(
\frac{N_0}{\widetilde N_c^\circ(T)}
\right)^{2/5}
\right].
\]

A single raw repetition is not required to lie exactly on the global fit. Each
trajectory will therefore use the anchored incremental form

\[
N_{\rm th}(T,N_0)
=
N_{{\rm th},0}
+F_{\rm ns}(T,N_0)
-F_{\rm ns}(T_0,N_{0,0}).
\]

This equation guarantees exact recovery of the selected initial triplet and
uses the Oxford result only for the subsequent local change in thermal
population. It adds no fitted parameter.

The implementation must reject states outside the evidence domain. The first
contract will use

\[
0.0423
\leq
\left(\frac{N_0}{\widetilde N_c^\circ}\right)^{2/5}
\leq
0.4708
\]

and

\[
0.90T_0\leq T\leq1.10T_0.
\]

The \(x\)-coordinate bounds are the minimum and maximum accepted 0-degree
coordinates in the public averaged Figure 2 data. The temperature gate encodes
the approximately 10% rescaling range used in the source analysis as the exact
implementation interval above. A rejected state is unsupported; it must not be
silently extrapolated.

As a same-data implementation check, the maintained 0-degree line gives a
relative RMS residual of approximately 3.1%, a mean absolute residual of 2.2%
and a maximum residual of 7.4% over the accepted averaged points. These values
verify data import and normalisation; they are not independent model
validation.

## 5. Minimal energy closure

The main model conserves trapped atom number:

\[
N_{\rm trap}=N_0+N_{\rm th}={\rm constant}.
\]

It represents heating as redistribution from the condensate into the thermal
component. The thermal energy is

\[
E_{\rm th}=c_E N_{\rm th}k_BT.
\]

Only two physically interpretable endpoints will be run:

\[
c_E=3\frac{\zeta(4)}{\zeta(3)}=2.70118\ldots
\]

for a saturated ideal Bose gas in a harmonic trap, and

\[
c_E=3
\]

for the classical harmonic-gas limit. These endpoints define a two-model
sensitivity span. They have not been shown to be lower and upper bounds for an
interacting, non-saturated dipolar gas. The coefficient \(c_E\) is not a fitted
apparatus parameter.

A preliminary reference-state scale check places the condensate
Thomas--Fermi interaction energy at approximately 0.5% of the total energy in
this bimodal cloud. The implementation must reproduce and record this estimate
from

\[
E_0^{\rm TF}=\frac{5}{7}\mu N_0,
\qquad
\frac{E_0^{\rm TF}}{E_0^{\rm TF}+E_{\rm th}},
\]

using the central triplet, the maintained dipolar Thomas--Fermi chemical
potential and each \(c_E\) endpoint. Its change at 30% condensate depletion is
expected to be a few per cent of the required thermal-energy increase. The
first frame engine will omit it. No condensate interaction term will be added
unless the recorded sensitivity changes the depletion-frame span materially.

The non-saturation closure and the energy expression do not constitute a full
finite-temperature interacting equation of state. The two \(c_E\) endpoints and
the validity gates are the declared control for this approximation.

## 6. Recoil-limited probe disturbance

Let \(n_\gamma^{(2{\rm lvl})}\) be the photons scattered per atom predicted by
the existing unit-strength two-level expression. The scattering response is
method dependent:

\[
n_{\rm sc}=\chi_{\rm sc}n_\gamma^{(2{\rm lvl})}.
\]

For the current axial, fully spin-polarised, linearly polarised Faraday case,

\[
\chi_{\rm sc}
=
\frac{s_-+s_+}{2}
=
\frac{1+1/91}{2}
=
\frac{46}{91}.
\]

This is distinct from the signed rotation coefficient

\[
\kappa_F=\frac{s_+-s_-}{2}=-\frac{45}{91}.
\]

The multiplicative \(46/91\) form is approved only when

- \(n_\gamma^{(2{\rm lvl})}\) is normalised to the total incident linear
  intensity driving a unit-strength cycling transition;
- the incident linear polarisation contains equal \(\sigma^-\) and
  \(\sigma^+\) intensities;
- the two branches use the same dispersive detuning;
- excitation is weak enough that the scattering rate is linear in branch
  strength.

If these conditions fail, the two branch scattering rates must be evaluated
and summed explicitly. The implementation must not apply \(46/91\) on top of a
rate that already contains the branch strengths.

The deposited recoil energy in frame \(q\) is

\[
Q_q
=
N_{\rm trap}
n_{{\rm sc},q}
(2E_{\rm rec}).
\]

The frame update is

\[
E_{q+1}=E_q+Q_q.
\]

PCI or DGI may use a different polarisation and transition branch. Their
\(\chi_{\rm sc}\) values must be derived from their declared optical geometry;
they must not inherit \(46/91\) automatically.

The model assumes that the probe covers every atom included in
\(N_{\rm trap}\). It also assumes that the trap is deep enough for recoil events
not to remove atoms. Until beam coverage and trap depth are supplied, the
fixed-number result remains an ideal screening estimate.

## 7. Per-frame numerical problem

For a selected initial triplet and \(c_E\), initialise

\[
N_{\rm trap}=N_{0,0}+N_{{\rm th},0},
\qquad
E_0=c_E N_{{\rm th},0}k_BT_0.
\]

After adding \(Q_q\), solve the following three equations for
\((N_{0,q+1},N_{{\rm th},q+1},T_{q+1})\):

\[
N_{0,q+1}+N_{{\rm th},q+1}=N_{\rm trap},
\]

\[
N_{{\rm th},q+1}
=F_{\rm ns}^{\rm anchored}(T_{q+1},N_{0,q+1}),
\]

\[
E_{q+1}
=c_E N_{{\rm th},q+1}k_BT_{q+1}.
\]

A one-dimensional bounded root solve is sufficient because number conservation
eliminates one population and the anchored closure eliminates the other. A
general thermodynamic framework, kinetic solver or free-form callback system is
outside the approved implementation.

### 7.1 Interface to the optical model

The thermodynamic state engine and optical forward model remain separate. For
the first replacement of Figures 5.2 and 5.4:

- the optical layer receives the updated \(N_0\);
- it recomputes the condensate Thomas--Fermi chemical potential, radii and
  column-density profile using the existing Oxford core geometry;
- \(N_{\rm th}\) and \(T\) remain reported thermodynamic state variables and do
  not generate a per-frame thermal-halo image;
- the optical signal and SNR therefore describe the condensate-core observable,
  conditional on the broad thermal halo being excluded from the image-plane
  signal model.

This boundary must be stated in the Figure 5.2 and 5.4 metadata. A future
temperature or total-number observable would require an explicit thermal-halo
forward model and is outside this replacement.

## 8. Equilibration condition

The state update is valid only when the cloud re-equilibrates between probe
pulses. A preliminary estimate based on

\[
\tau_{\rm th}=\frac{\alpha}{\bar n\,\bar\sigma\,v_r}
\]

and the cross-dimensional Er model gives 63--65 ms for the current reference
scale. The Ultra quiet catalogue rate gives a minimum full-resolution camera
period of 184.4 ms. These numbers make equilibration plausible for the declared
full-frame scenario, but they do not validate it experimentally. The density
dependence reported for Er and the possibly distinct condensate--thermal
exchange time prevent a general numerical cadence threshold.

The first implementation will therefore support exactly one conditional timing
contract:

- `full_resolution_ultra_quiet`, with a nominal minimum period of 184.4 ms;
- status `conditional_on_rethermalisation`;
- every faster, sub-array or otherwise unverified cadence returns
  `unsupported_cadence`.

The state engine will not compare cadence to the preliminary 63--65 ms estimate
or present the estimate as a calibrated limit. It will not simulate a hot
reservoir or fit a thermalisation time from the Oxford hold-time data.

## 9. Approved outputs

For each of the three initial triplets, each \(c_E\) endpoint and each declared
fluence \(F=90,150,300\ {\rm mW\,\mu s}\), report

- frame index;
- \(T_q\);
- \(N_{0,q}\);
- \(N_{{\rm th},q}\);
- \(N_{0,q}/N_{0,0}\);
- cumulative scattered photons per atom;
- cumulative recoil energy per trapped atom;
- the first excluded frame under the 30% condensate-depletion rule;
- any validity-gate failure.

The primary computation contains six thermodynamic trajectories per fluence:
three initial triplets multiplied by two energy endpoints. Their range is a
**declared sensitivity span**, not a statistical confidence interval or a
guaranteed physical bound.

Define condensate depletion and trapped-atom loss separately:

\[
D_q=1-\frac{N_{0,q}}{N_{0,0}},
\qquad
L_q=1-\frac{N_{{\rm trap},q}}{N_{{\rm trap},0}}.
\]

The fixed-number model has \(L_q=0\) by construction. The retained 30% rule is
\(D_q\geq0.30\). The threshold-crossing pulse is excluded, consistent with the
current strict integer convention. The output name must be
`condensate_depletion_frame`; it must not imply measured atom loss.

`N_use` remains a separate combination of the depletion sensitivity span and the
declared observable-usability criterion. The thermodynamic model does not
redefine reconstruction support or make accumulated SNR a primary observable.

## 10. Deliberately excluded physics

The first implementation will not include

- a finite-temperature dipolar Hartree or HFB--Popov solver;
- condensate and thermal interaction energies in the per-frame update;
- a hot-atom or two-reservoir kinetic model;
- evaporation or finite trap-depth escape;
- background, three-body or light-assisted loss;
- Zeeman-state optical-pumping dynamics;
- a dynamic core-plus-thermal reabsorption calculation;
- an apparatus heating multiplier;
- a fitted Oxford hold-time heating law;
- full thermal-halo image generation after every frame.

These exclusions are model boundaries. They are not assumed to be absent from
an experiment. Light-assisted collisions, stray optical lattices and finite
trap depth require apparatus data. Adding uncalibrated coefficients for them
would increase apparent detail without improving the present prediction.

Two small sensitivity checks are permitted without changing the main model:

1. vary the reabsorption energy correction from 0 to 5%;
2. add the independently calculated condensate interaction-energy change.

If either changes the depletion-frame span by no more than one frame, the
effect will be reported as negligible at the present screening resolution and
will remain outside the frame engine.

## 11. Implementation boundary

The change will be limited to one active thermodynamics implementation and one
test surface.

1. Extend the existing Oxford initial-condition record with the two additional
   300 ms triplets and their source metadata. Do not create a second competing
   sample definition.
2. Implement one small multiframe thermodynamics module containing the
   measured-scale closure, anchored closure, recoil update and bounded state
   solve.
3. Make Figures 5.2 and 5.4 call the same sequence function. Remove their
   duplicated active thermodynamic loops only after both consumers agree with
   the shared-engine states for identical inputs.
4. Keep `scripts/recover_notebook_multishot_stage.py`, notebook-section code and
   notebook-aligned regression outputs unchanged as historical Version 1
   reproduction.
5. Generate one diagnostic table before regenerating any dissertation figure.
   Review all six trajectories at \(F=90,150,300\ {\rm mW\,\mu s}\) before
   accepting replacement values.
6. Update the dissertation figures, captions, parameter register and prose only
   after the diagnostic is accepted.

No abstract plugin architecture, model registry or general-purpose
thermodynamic callback layer is required.

## 12. Verification and acceptance tests

The implementation must pass the following tests.

### State and numerical tests

- A zero-energy pulse leaves the selected triplet unchanged.
- Every one of the three initial triplets is recovered within \(10^{-6}\) atom
  and \(10^{-6}\ {\rm nK}\).
- The relative residual in \(N_0+N_{\rm th}=N_{\rm trap}\) is at most
  \(10^{-10}\) after every accepted frame.
- The relative energy-equation residual is at most \(10^{-10}\), using the
  larger of \(|E_{q+1}|\) and \(k_BT_{q+1}\) as the denominator scale.
- Positive recoil energy produces non-decreasing temperature.
- Populations remain finite and non-negative.
- The solver stops with an explicit reason outside the closure domain.

### Physics and provenance tests

- Setting \(C_m=1\), \(S_m=0\) and removing the anchor reproduces the ideal
  saturated population limit used only as a comparator.
- The measured-scale closure reproduces the accepted Oxford 0-degree
  non-saturation data with relative RMS and mean absolute residuals within
  0.1 percentage point of 3.1% and 2.2%, respectively, and a maximum residual
  no greater than 7.5%.
- The Faraday scattering calculation uses \(\chi_{\rm sc}=46/91\) exactly once.
- The Faraday signal continues to use signed \(\kappa_F=-45/91\); no code path
  substitutes one coefficient for the other.
- Figures 5.2 and 5.4 return identical thermodynamic states for identical pulse
  sequences.
- Stored Version 1 outputs remain labelled and are not silently overwritten.

### Acceptance rule

The minimum model is sufficient when

- all numerical and provenance tests pass;
- the two \(c_E\) endpoints and three initial triplets produce a finite,
  interpretable depletion-frame sensitivity span;
- the permitted reabsorption and condensate-energy checks change each span
  boundary by no more than one frame;
- no reported state leaves the Oxford closure domain before its declared
  stopping point.

If an approved sensitivity changes a span boundary by more than one frame,
report the wider span first. Add a more detailed physical model only when
the dissertation conclusion depends on narrowing it and suitable evidence is
available.

## 13. Dissertation claim language

Allowed descriptions include

- recoil-limited multiframe estimate;
- quasi-equilibrium screening model;
- fixed-trapped-number condensate-depletion sensitivity span;
- conditional prediction under the Oxford non-saturation closure;
- model-supported frame range.

Do not describe these outputs as

- measured heating or measured atom loss;
- an apparatus-calibrated lifetime;
- a complete finite-temperature dipolar-gas simulation;
- validation against the Oxford hold-time trajectory;
- a precise experimental \(N_{\rm use}\).

The central interpretation is that repeated DPFI trades information against a
declared recoil disturbance. The calculation estimates how this disturbance
redistributes a measured bimodal cloud under a local equilibrium closure. It
does not establish how unmodelled apparatus effects alter the sequence.

## 14. Source basis

- M. Krstajić *et al.*, *Interaction shift of the Bose--Einstein condensation
  temperature in a dipolar gas*, Phys. Rev. A **111**, L051303 (2025),
  <https://doi.org/10.1103/PhysRevA.111.L051303>.
- Supplementary derivation and non-saturation equations,
  <https://arxiv.org/pdf/2501.16318>.
- Oxford public data associated with the publication,
  <https://doi.org/10.5287/ora-m8gpvdr2y>.
- A. Patscheider *et al.*, *Determination of the scattering length of erbium
  atoms*, Phys. Rev. A **105**, 063307 (2022),
  <https://doi.org/10.1103/PhysRevA.105.063307>.
- E. Altuntaş and I. B. Spielman, *Weak-measurement-induced heating in
  Bose--Einstein condensates*, Phys. Rev. Research **5**, 023185 (2023),
  <https://doi.org/10.1103/PhysRevResearch.5.023185>.

The Oxford paper and data control the initial state and equilibrium population
closure. Patscheider *et al.* support only the order-of-magnitude
rethermalisation check. Altuntaş and Spielman establish that recoil-only
heating can miss apparatus-specific disturbance; their numerical excess-heating
factors are not transferred to erbium.
