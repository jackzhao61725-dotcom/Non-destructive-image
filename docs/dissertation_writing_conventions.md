# Dissertation Writing Conventions

- **Status:** active prose authority
- **Applies to:** the MSc in Quantum Technologies dissertation on non-destructive imaging of dipolar `166Er` condensates
- **Established:** 22 July 2026
- **Scope:** argument, exposition, equations, figures, source use, terminology and prose review
- **Active consumers:** dissertation drafting, section review and figure QA
- **Update trigger:** an approved change to the dissertation structure, physical
  conventions, reported observables, course requirements or evidence standard
- **Retirement rule:** replace only with one indexed successor that absorbs all
  still-current rules

This file controls how dissertation material is written. Scientific values and
model contracts remain controlled by the latest approved dissertation text,
the current hand-off and the verified repository evidence. A writing rule must
never be used to smooth over a physical inconsistency.

## 1. Project AI-use status and authorship

The author has confirmed that the QTech group and college permit unrestricted
use of AI tools for this research project and support their use as research
infrastructure. This project-specific instruction controls the workflow.

The anti-AI prose gate in Section 13 is a quality standard. Its purpose is to
remove empty rhetoric, vague claims and formulaic transitions. It is not an
attempt to conceal tool use.

The author remains responsible for every equation, source, numerical value and
claim. AI output is not evidence. A source must be opened, read and checked
before it is cited. Tool-assisted calculations and code must remain
reproducible from the recorded inputs and repository state.

The formal course boundary remains the
[Oxford MSc in Quantum Technologies Examination Regulations](https://examregs.admin.ox.ac.uk/Regulation?code=mosbcinquantech).

## 2. Genre, examiner and reader

The dissertation reports an independent four-month MSc research project. The
formal word limit is 20,000 words, and the examiners may require a viva. These
conditions set the required depth.

The assumed reader has a strong undergraduate background in physics,
engineering or a related quantitative subject. The reader may know quantum
mechanics, basic optics, Fourier transforms, probability and experimental
uncertainty. The reader is not assumed to know:

- dipolar Bose-Einstein condensates;
- the `166Er` level structure or Oxford apparatus;
- dispersive cold-atom imaging;
- the repository or its historical notebook conventions;
- inverse problems for finite-aperture camera data.

Explain specialised physics at the point where it changes the argument. Do not
re-teach standard undergraduate material. A familiar result may be stated with
a source and a short physical interpretation.

The quantum-technology connection must be operational. Use measurement
backaction, information per disturbance, repeated observation, optical
hardware, calibration and parameter uncertainty. Generic claims about a
"quantum revolution" add no scientific content.

## 3. Core standard: every sentence must do work

A sentence should perform at least one of these functions:

- define a quantity or convention;
- state a physical or methodological claim;
- supply evidence through an equation, result, figure or source;
- state an assumption or validity condition;
- interpret evidence for the research question;
- make a necessary logical transition.

Delete a sentence if removing it changes no definition, evidence, boundary or
conclusion.

The dissertation's main argumentative chain is:

```text
physical question
  -> required observable
  -> measurement channel
  -> forward-model assumptions
  -> predicted or observed data
  -> inference and uncertainty
  -> physical or experimental consequence
```

Do not reverse this chain by opening with code, an optimiser or a plotted
quantity whose physical role has not been established.

## 4. Scale of argument

### 4.1 Chapter

Each chapter has one research task. Its opening paragraph states:

- the question;
- the input inherited from earlier chapters;
- the output produced here;
- one material scope boundary, if needed.

The final paragraph answers the opening question and establishes the next
dependency. It does not repeat the table of contents.

### 4.2 Section

Each section resolves one part of the chapter task. A section title should name
the physical object, operation or decision. Avoid titles such as "Further
discussion" or "Additional considerations".

### 4.3 Paragraph

One paragraph has one argumentative purpose. A useful pattern is:

1. claim or question;
2. equation, data, figure or primary source;
3. interpretation;
4. validity condition or consequence, when needed.

The pattern is a reasoning check, not a sentence quota. A short paragraph may
need only two sentences. Split a paragraph when its second half answers a new
question.

Use a topic sentence with physical content. Do not open a routine subsection
with "This section discusses...". Chapter openings and essential bridges may
use concise signposting.

## 5. Dissertation structure and chapter contracts

The active structure is:

1. **Introduction:** define the measurement problem, the literature gap, the
   research questions, the contribution and the scope.
2. **Dipolar condensates and target observables:** introduce the `166Er`
   physics needed to motivate the observables; define the reference
   condensate and the limits of its contact-only Thomas-Fermi description.
3. **Non-destructive imaging physics:** derive the scalar and vector responses,
   optical transfer functions, noise scaling and disturbance scaling.
4. **Forward model:** follow the measurement chain from condensate state to
   optical field, finite aperture, raw camera channels, detector noise and
   repeated-exposure evolution.
5. **Acquisition screening:** use explicit per-frame criteria and condensate
   evolution to identify conditions worth passing to the inverse analysis.
6. **Physics-informed inference:** define the raw-channel likelihood, latent
   nuisance representation, physical supports, observables, identifiability
   and uncertainty contract.
7. **Observable estimates and recoverability:** report bias, interval width,
   coverage and unsupported cases for the declared observables.
8. **Initial DPFI implementation plan:** translate the model into probe,
   analyser, detector and commissioning requirements; mark provisional inputs
   and required calibrations.
9. **Conclusions and outlook:** answer the research questions using results
   already established; place TOF, more complex states and sequence-coupled
   inference in the outlook.

The chapter boundary is strict:

- Chapter 3 owns the imaging physics.
- Chapter 4 explains the implemented measurement operator without re-deriving
  Chapter 3.
- Chapter 5 screens acquisition conditions. Its camera-space SNR gate is not a
  final recoverability criterion.
- Chapters 6 and 7 decide which physical observables the raw counts support.
- Chapter 8 is an implementation plan until experimental results exist.

## 6. Physical exposition

### 6.1 Start from the exact case

Use general theory only when it is required to derive or interpret the
`166Er` case. Once the relevant selection rules, geometry and state preparation
are fixed, calculate that case directly. Do not add population-averaged,
arbitrary-`J` or full-tensor formulae solely to display generality.

The reference condensate is the worked example carried through the
dissertation. Introduce its parameters once. Later sections use the name
"reference condensate" and cite the defining section. A new sample appears
only when it tests a distinct physical or numerical boundary.

### 6.2 Define one convention once

Every symbol has one meaning. If two normalisations are common in the
literature, select the one used by the model and state how it relates to the
measured quantity. Do not present both as interchangeable conventions.

For the current imaging chapters:

- `\(\phi(y,z)\)` is the scalar phase object defined in Chapter 3;
- `\(\Phi_\pm(y,z)\)` are the total phases of the two circular components;
- `\(\theta_F=(\Phi_+-\Phi_-)/2\)` is the Faraday rotation;
- `\(\kappa_F=\theta_F/\phi\)` uses the same scalar-phase convention;
- `\(F=P\tau\)` is a power-time coordinate. Call it physical fluence only
  when the beam area is fixed and stated.

Do not revive the obsolete illustrative `\(\kappa_F=1\)` conversion in new
dissertation prose. The approved `166Er` calculation and its stated
approximations control Section 3.3.

### 6.3 Derive to the level required for interpretation

Retain a derivation when it does at least one of the following:

- fixes a sign or normalisation;
- exposes a scaling used later;
- separates two experimentally distinct signals;
- identifies an approximation or a failure mode;
- supplies the measurement operator used by the simulation or inference.

Compress intermediate algebra that has no later role. A citation does not
replace the explanation of a non-standard convention used by this work.

### 6.4 Numerical examples

A numerical substitution follows the equation it tests. State the physical
case, substitute the values, give the result with units and then state the
assumptions that control its interpretation. Do not number a one-off
substitution.

Report only justified precision. Literature constants, apparatus inputs,
fitted parameters and derived values must remain distinguishable.

### 6.5 Assumptions and limits

State an assumption where it first affects an equation or conclusion. A useful
scope statement contains:

1. the conditions under which the result is used;
2. the omitted effect;
3. the consequence of that omission.

Do not attach a generic limitations paragraph to every section. Do not list
possible complications that do not affect the present calculation.

## 7. Equations and notation

Introduce every displayed equation with a grammatical sentence. Treat the
equation as part of the sentence and punctuate it accordingly.

Number an equation only if it is:

- a core definition or model contract;
- used in a later calculation;
- discussed or compared later;
- needed for a cross-reference.

Leave local identities, intermediate algebra and numerical substitutions
unnumbered.

After a core equation, explain its physical content. Useful explanations state
a sign, scaling, limiting case, measurable quantity or model dependency. Do
not translate the equation into a sentence that says exactly the same thing.

Define symbols near first use. Give units for dimensional quantities. If more
than three new symbols appear together, group their definitions by physical
role or split the derivation.

Use these forms consistently:

- isotope: `\(^{166}\mathrm{Er}\)`;
- numbered item: `Equation (3.4)`, `Figure 3.2`, `Table 5.1`, `Section 3.3`;
- value and unit: `401 nm`, `29.5 MHz`, `0.1004 rad`, `5.75^\circ`;
- leading zero: `0.203`, not `.203`;
- uncertainty: estimate and uncertainty use compatible decimal places.

## 8. Evidence and claim strength

Choose verbs from the evidence available.

| Evidence | Appropriate wording | Wording that needs stronger evidence |
|---|---|---|
| definition or declared convention | define, adopt, set | prove, establish physically |
| analytical result under stated assumptions | follows, gives, predicts | universally shows |
| verified code against an analytic or frozen reference | reproduces, agrees with, verifies the implementation | validates the physical model |
| truth-known synthetic test | recovers, estimates, passes in the stated ensemble | measures, experimentally validates |
| parameter scan | screens, indicates, identifies a candidate region | determines the optimum |
| experimental calibration data | calibrates, estimates experimentally | universal coefficient |
| held-out experimental test | validates under the tested conditions | validates generally |
| literature experiment | reported, measured, demonstrated | proves for this apparatus or species |

Use `prove` only for a mathematical result. Use `significant` only with a
defined statistical or physical meaning. Use `robust` only after naming the
perturbations tested and the acceptance criterion.

Keep these distinctions explicit:

- verification tests an implementation against its specification;
- calibration estimates model or instrument parameters;
- validation tests frozen predictions against independent evidence;
- screening removes unsuitable operating regions;
- inference estimates declared observables from measured channels;
- identifiability asks whether the data constrain an observable;
- resolution, accuracy, precision and uncertainty are different quantities.

## 9. Reconstruction and image language

In this dissertation, reconstruction means model-based inference of physical
observables from the raw camera channels. The reported primary vector is

```text
q = (A, y_c, z_c, w_major).
```

Each component carries an uncertainty, a data-consistency interval or an
unsupported flag. The shape-flexible density field is an internal nuisance
representation required by the nonlinear forward operator. Do not present it
as a super-resolved density image or judge success by visual sharpness.

Peak density, peak pixel value, full-map error, minor width, aspect ratio and
principal-axis angle are diagnostics unless a section explicitly promotes one
of them for a stated physical question and demonstrates support from the data.

RAI is the terminal reference for a detailed density measurement. DPFI is
studied for repeated, lower-disturbance estimates of spatial observables. State
these capabilities directly; do not turn them into an image-quality contest.

## 10. Literature and citations

### 10.1 Give every source a role

Use sources for a specific function:

- a review or advanced text establishes accepted background and notation;
- a theory paper supports a derivation or bound;
- an instrument paper supports an implemented measurement and its conditions;
- an apparatus paper or Oxford dissertation supports local parameters and
  procedures;
- an analysis paper supports an estimator or inverse method.

Do not transfer a numerical calibration between species, transitions,
geometries or apparatuses without a derivation that permits the transfer.

### 10.2 Place citations at the supported claim

The citation follows the clause or sentence it supports. A citation cluster at
the end of a long paragraph is acceptable only when every source supports the
same claim. If sources support different parts, split the sentence and place
them separately.

Use an author as the grammatical subject when the research action matters:

> Gajdacz et al. measured the Faraday signal and inferred scattering from the
> observed heating.

Use a parenthetical citation when the physical statement matters more than the
authors:

> Repeated Faraday imaging has followed the evolution of one atomic cloud
> across many exposures (Gajdacz et al., 2013).

### 10.3 Prefer primary evidence

Use the original paper for an exact result, calibration or experimental
demonstration. A review may support broad context. An Oxford DPhil dissertation
may support apparatus-specific information that is absent from a paper, with
the provenance stated.

### 10.4 Citation closure

Before submission:

- every in-text citation has one complete bibliography entry;
- every bibliography entry is cited for a real claim;
- author names, year, title, journal, volume, pages or article number and DOI
  match the primary record;
- diacritics and hyphenation are preserved;
- the cited source has been opened and the supporting passage checked.

The current draft does not yet pass this closure test. Citation cleanup is a
separate factual task; do not hide it inside prose editing.

## 11. Figures, tables and captions

Each figure has one main point. Draft that point before choosing the panels.
Remove a panel that does not support it.

Figure-internal text contains scientific identifiers only. A title or panel
label names the plotted quantity, sample or comparison; it does not narrate a
workflow. Do not use `Stage X`, `Step X`, `notebook-aligned recovery`,
`reference implementation` or explanatory subtitles inside a figure.

Use the dissertation's physical quantity, not the implementation variable:

| Quantity | Figure label | Boundary |
|---|---|---|
| full column-density map | `\(n_{\mathrm{col}}(y,z)\)` or the corresponding imaging plane | reserve `\(\tilde n_x\)`, `\(\tilde n_y\)` and `\(\tilde n_z\)` for peak scalar values |
| scalar phase | `\(\phi\)` in radians | do not call it intensity |
| Faraday rotation | `\(\theta_F\)` in radians | do not call it intensity or a camera signal |
| normalised optical intensity | `\(I/I_0\)` | state the reference used for `\(I_0\)` in metadata or the caption |
| dual-port difference | `\(S=(I_H-I_V)/(I_H+I_V)\)` | `\(S\)` is a dimensionless signal, not intensity or camera counts |
| detected raw data | photon counts or photoelectron counts | state which quantity, gain convention and background correction are used |
| scattering per exposure | `\(N_\gamma\)` | distinguish it from the rate `\(R_{\mathrm{sc}}\)` and state `\(N_\gamma=R_{\mathrm{sc}}\tau\)` when both appear |

The historical synthetic Faraday outputs map notebook port `v` to `H` and
port `u` to `V`. This mapping is provenance for those frozen outputs. It does
not fix the sign or port names of a future experiment; Chapter 8 must determine
those from calibration data.

Use math typesetting for symbols and units in figure text, including
`\(\mu\mathrm{m}\)`, `\(\mathrm{m}^{-2}\)`, `\(\mu\mathrm{s}\)` and
`\(e^-\)`. Round displayed values to the precision needed to read the figure.
Exact values belong in the accompanying CSV, JSON, table or caption.

A caption states:

1. what is shown;
2. the sample and operating conditions;
3. axes, units, normalisation, colour and marker meanings;
4. the estimator and uncertainty or sample-count convention, when relevant;
5. one conclusion supported by the figure;
6. a material boundary needed to read it correctly.

For a normalised or processed image, the figure metadata also records the
plotted quantity, reference image or value, whether the data are absolute,
normalised, binned, noisy or count-valued, where unnormalised data are stored,
and whether the response is calibrated. A caption need not repeat every
metadata field, but it must include every condition needed to evaluate the
claim made in the main text.

Use `Figure n.n.` and `Table n.n.` with a full stop. Do not write "the figure
above" or "the table below". Refer to the numbered object.

Do not use a peak-pixel metric to support a claim about integrated response,
centroid or width. Do not describe finite-ensemble ranges as confidence
intervals. A figure caption cannot strengthen a claim beyond the analysis in
the main text.

Tables are for repeated records with shared fields. Do not use a table to hold
ordinary prose. Report units in the column headings and define all symbols in
the caption or nearby text.

## 12. Voice, tense and English usage

Use British English:

- polarisation;
- normalised;
- centre;
- modelling;
- optimisation;
- characterisation.

Use `dissertation` in official prose. Use `thesis` only in the title or type of
a cited DPhil work when that source uses the term.

Prefer active subjects that identify the agent:

- `The model propagates the field...`
- `The fit estimates the centroid...`
- `Figure 5.2 shows...`
- `This dissertation adopts...`

Avoid `we` in the single-author dissertation. First-person singular may be
used if the course permits it and authorship is material; the current house
style normally uses `this dissertation`, `this work`, `the model` or `the
experiment`.

Use present tense for established theory, definitions, figures and conclusions
that remain true in the document. Use past tense for completed experimental or
computational actions. Use future or conditional tense for Chapter 8 tasks that
have not been performed.

Keep the subject and verb close. A sentence longer than about 30 words is a
review flag. Split it unless its logical structure remains clearer as one
sentence. Replace noun-heavy phrases with direct verbs: `calculate`, not `make
a calculation of`; `assess`, not `make an assessment of`.

## 13. Anti-AI prose gate

The following constructions are review failures by default. Retain one only
when it carries a necessary, explicit scientific relation.

### 13.1 Empty comparison

**Flag:** `X rather than Y`.

Keep the comparison only if both alternatives are real, use the same stated
criterion and lead to different physical or design conclusions. Otherwise,
state the selected object and its scope in separate direct sentences.

Bad:

> The reconstruction aims to recover meaningful physical information rather
> than merely produce a visually appealing image.

Use:

> The reported outputs are `A`, `(y_c,z_c)` and `w_major`, each with an
> uncertainty or support flag. The latent density field is an internal nuisance
> representation.

### 13.2 Decorative expansion

**Flag:** `not only X but also Y`, `not merely`, `not simply`, `more than just`.

Write the two claims separately and give their causal or quantitative
connection.

Bad:

> DPFI not only reduces destruction but also enables repeated imaging.

Use:

> Lower scattering leaves the condensate available for another exposure. The
> usable sequence ends when the framewise observable criterion or disturbance
> limit fails.

### 13.3 Empty importance claims

**Flag:** `important`, `crucial`, `key role`, `cannot be overstated`,
`valuable insight`.

Replace the evaluation with the mechanism or consequence.

Bad:

> Numerical aperture plays a crucial role in reconstruction.

Use:

> The pupil removes spatial frequencies above `NA/lambda`; the lost components
> cannot be restored from the recorded channels without additional prior
> information.

### 13.4 Vague retrospective claims

**Flag:** `This highlights...`, `This underscores...`, `This demonstrates...`,
`Taken together...`, `In this context...`.

Name the evidence and conclusion.

Bad:

> This highlights the robustness of the method.

Use:

> Across the stated perturbation range, the centroid bias remained below
> `0.650 um` in all held-out trials.

### 13.5 Promotional or generic research language

**Flag:** `robust framework`, `comprehensive`, `powerful`, `novel`, `meaningful`,
`offers insight`, `lays the foundation`, `paves the way`, `landscape`, `delve`,
`leverage`, `facilitate`.

Replace these words with the implemented action, tested condition or reported
quantity. Use `novel` only when a literature review establishes the novelty and
the dissertation needs the claim.

### 13.6 Reader coaching and obviousness

**Flag:** `clearly`, `obviously`, `interestingly`, `it should be noted`, `it is
worth noting`, `it is evident`.

Give the equation, number or source that makes the conclusion credible.

### 13.7 Meta-narration

**Flag:** `This section discusses/presents/explores...` and repeated previews of
the next section.

Open with the physical question or result. Keep one chapter-level contract and
only the bridges needed to preserve dependency.

### 13.8 False contrast

**Flag:** `while`, `although`, `whereas`, `by contrast`, `instead` when the two
clauses do not conflict.

Use a contrast marker only when the second clause limits or changes the first.
Two compatible facts should be two direct clauses or sentences.

### 13.9 Inflated purpose phrases

**Flag:** `in order to`, `serves to`, `is designed to allow`, `with the aim of`.

Use `to` plus a direct verb.

### 13.10 Model personification

Do not write that a model `understands`, `discovers` or `reveals` a physical
truth. A model predicts, estimates, constrains, reproduces or fails under stated
conditions.

## 14. Real-comparison test

A scientific comparison is admissible when all four conditions hold:

1. both objects are defined;
2. the metric and normalisation are common;
3. the operating conditions are stated;
4. the result changes an interpretation or decision.

Example:

> At equal incident photon number and the stated reference angle, the two
> readouts give different camera-count distributions. Their observable
> precision is compared with the same likelihood and support.

Do not infer a parameter-estimation advantage from a raw intensity ratio alone.

## 15. Project terminology

Use the following terms consistently:

- **resonant absorption imaging (RAI)**;
- **phase-contrast imaging (PCI)**;
- **dark-ground imaging (DGI)**;
- **dark-field Faraday imaging (DFFI)**;
- **dual-port Faraday imaging (DPFI)**;
- **non-destructive imaging:** operationally leaves the same cloud available
  for another observation; it does not mean zero scattering;
- **reference condensate:** the single defined Oxford baseline used for worked
  calculations;
- **raw channels:** the detected camera counts before latent-object inference;
- **latent field:** the internal shape-flexible object used by the forward
  operator;
- **observable:** a declared physical summary estimated from the measurement;
- **supported:** constrained under the stated data and inverse contract;
- **unsupported:** not assigned a regularisation-selected physical value.

Define every acronym at first use in the abstract and again at first use in the
main text if the separation is large. Do not alternate `non-destructive` and
`nondestructive` in dissertation prose; retain a source's spelling in its
title.

## 16. Chapter-specific prose checks

### Introduction

- Open with the measurement problem, not generic quantum-technology history.
- Identify what existing methods measure and what remains unresolved.
- State contributions at the level actually completed.
- Keep literature claims mapped to primary sources.

### Dipolar-condensate background

- Include physics that motivates `A`, centroid and `w_major` or sets the sample
  range.
- Do not imply that a literature phase, such as a supersolid or roton regime,
  is present in the Oxford sample without direct evidence.
- Separate the contact-only reference generator from dipolar physical claims.

### Imaging theory

- Explain how the atomic response becomes the measured intensity.
- Fix sign, helicity and normalisation conventions before numerical use.
- Carry one reference-condensate calculation through the section.
- State optical and atomic approximations immediately after the result they
  qualify.

### Forward model

- Follow the physical measurement chain.
- Separate physical inputs, numerical choices and fitted/calibrated values.
- Call regression against an analytic result verification.

### Acquisition screening

- State the estimator, spatial support, noise model and stopping rule.
- Treat SNR as a declared camera-space gate.
- Do not accumulate unlike frames without an explicit dynamical parameter
  model.

### Inference and results

- Start from the raw-channel likelihood and declared observables.
- Report estimate, bias, uncertainty or interval, and support status.
- Separate detector uncertainty, inverse ambiguity and calibration uncertainty.
- Report negative and unsupported cases.
- Do not use image appearance as the main evidence.

### Implementation plan

- Distinguish existing hardware, proposed hardware and commissioning
  measurements.
- Use conditional tense for unbuilt components and unmeasured performance.
- Do not convert simulation inputs into apparatus specifications.

### Conclusion

- Answer each research question with an established result.
- State the strongest justified contribution and its boundary.
- Introduce no new calculation, source or numerical result.

## 17. Drafting and review workflow

For each section:

1. write one sentence stating the section question;
2. list the evidence required to answer it;
3. verify each source and numerical input;
4. draft the equations and figures before connective prose;
5. write one paragraph per argumentative purpose;
6. run the claim-strength and anti-AI gates;
7. check notation, units and cross-references;
8. check citation closure;
9. read the section aloud;
10. apply the viva test: for every strong claim, identify the equation, result
    or source that answers "How do you know?"

Phrase scan for plain-text or Markdown drafts:

```powershell
rg -n -i 'rather than|not only|not merely|not simply|more than (just|a)|it is (important|worth noting|clear|evident)|this (highlights|underscores|demonstrates)|plays? (a )?(key|crucial) role|serves as|provides (a|an) (robust|comprehensive|valuable)|in order to|taken together|in this context|landscape|delve|leverage'
```

Every hit requires manual review. Do not perform blind replacement; a small
number of real scientific contrasts may be valid.

## 18. Source-derived writing lessons

The source register below controls how the reviewed literature is used. It is
not a substitute for the bibliography. Exact metadata and every supported
passage must still be checked against the primary record.

### 18.1 Condensate physics and erbium apparatus

| Source | Use and writing pattern | Boundary |
|---|---|---|
| Dalfovo et al. (1999) | Establish Gross--Pitaevskii and Thomas--Fermi approximations. State the controlling regime before the specialised equation, and identify column density as the line-of-sight observable. | The review does not justify a Thomas--Fermi ansatz for arbitrary modulated states; the approximation fails locally near an edge. |
| Ketterle, Durfee and Stamper-Kurn (1999) | Use the sequence physical picture, scale or criterion, then measurable quantity for introductory BEC and dispersive-imaging explanations. | Do not reproduce the lecture-note breadth, conversational style or historical notation. |
| Lahaye et al. (2009) | Support the dipolar Gross--Pitaevskii equation, `\(\epsilon_{dd}\)`, dipolar Thomas--Fermi background and conditional use of time of flight as a magnifier. | Use primary papers for experimental numbers and state the ballistic or interaction assumptions of a time-of-flight mapping. |
| O'Dell, Giovanazzi and Eberlein (2004) | Primary authority for the dipolar Thomas--Fermi and scaling solution under harmonic, spin-polarised, Thomas--Fermi conditions. | The ansatz does not establish the structure of droplet, roton or other modulated states. |
| Aikawa et al. (2012) | Support general erbium properties and a concise result-first pattern whose captions state time of flight, averaging, fit model and fitted quantities. | The experiment used `\(^{168}\mathrm{Er}\)` and is not an Oxford `\(^{166}\mathrm{Er}\)` apparatus authority. |
| Chomaz et al. (2016) | Use fit residuals, sample conditions and model comparison to separate an observed density distribution from a phase interpretation. | Do not import its beyond-mean-field droplet model into the present simplified inverse. |
| Chomaz et al. (2018) | Support the conditional mapping from a roton excitation to time-of-flight momentum peaks, including the simulation check, averaging and fit uncertainty. | A blurred in-situ ripple is not evidence of a roton; the paper's observable is the validated time-of-flight momentum structure. |
| Chomaz et al. (2019) | Separate density modulation from global phase coherence through distinct operational observables and state finite-resolution and ensemble limits. | Multiple peaks alone do not establish supersolidity. |
| Chomaz et al. (2023) | Supply current dipolar-gas context and consistent field notation. | Cite the relevant primary paper for a first observation, apparatus number or method-specific result. |
| Kučera (2024) | Primary local template for Oxford apparatus geometry, 401-nm absorption imaging, time of flight, fitting, calibration and uncertainty reporting. Follow purpose, geometry, model, actual operating limit, data product, fit and uncertainty. | Do not adopt DPhil length, informal wording, manual filtering without a criterion or an uncalibrated `fudge factor`. |
| Krstajić et al. (2023) | Preferred primary source for Oxford `\(^{166}\mathrm{Er}\)` apparatus numbers. Use sequence, independent calibration, selected fit interval, value with uncertainty and model inference. | A drift or cut-off from one data set is not a universal apparatus precision. |
| Krstajić et al. (2025) | Model interleaved measurement series, common-mode cancellation, time-of-flight corrections, binning and orthogonal-distance regression. Use the final title, *Interaction shift of the Bose--Einstein condensation temperature in a dipolar gas*. | Label any continuation beyond the model-valid region as empirical. Do not retain an obsolete pre-publication title. |
| Ulitzsch et al. (2017) | Provides a compact pattern for separating measured, calculated and expected apparatus quantities and validating a condensate with an independent observable. | Its `\(^{168}\mathrm{Er}\)` parameters belong to a different apparatus and cannot set the Oxford reference condensate. |

### 18.2 Atom--light response, Faraday imaging and disturbance

| Source | Use and writing pattern | Boundary |
|---|---|---|
| Andrews et al. (1996) | Support the optical-density scaling and the condensate-as-lens physical picture. | A dispersive signal does not imply zero scattering or zero backaction. |
| Lye et al. (2003) | Compare imaging modes at a common absorption or heating budget. | It does not establish a usable evolving-frame accumulated-SNR criterion. |
| Hope and Close (2004) | State the classical single-pass information bound at fixed spontaneous-scattering cost. | Do not present the bound as independent of the stated probe, detection and prior-information assumptions. |
| McClelland and Hanssen (2006) | Support the open 401-nm erbium level structure, optical-pumping pathways and the absence of a perfect closed cycle. | This paper is not the source of the measured 401-nm natural linewidth. |
| J. J. McClelland (2006) | Cite the early linewidth measurement `\(\Gamma/2\pi=35.6\pm1.2\,\mathrm{MHz}\)` and its separated statistical and systematic uncertainties. | If this work uses `29.5` or `29.7 MHz`, cite the source of that value and explain the selection; do not transfer the number from the McClelland--Hanssen paper. |
| Geremia, Stockton and Mabuchi (2006) | Move from the detector photocurrent through an explicit approximation ladder to scalar, vector and tensor atom--light terms. | Keep angular-momentum algebra not required by the `\(^{166}\mathrm{Er}\)` case out of the main text. |
| Reinaudi et al. (2007) | Separate the atomic optical depth from an apparatus response factor and calibrate the latter through an invariant observable; name polarisation, Zeeman population, viewport and detector effects. | Its rubidium calibration cannot be transferred to erbium. Apparatus-loss factors belong in the implementation plan, not the Section 3.3 atomic derivation or a Faraday-only correction to the comparison simulation. |
| Gajdacz et al. (2013) | Introduce the generic Faraday response, specialise to the prepared state, then separate atomic coupling, measured response and scattering-derived disturbance. | Its response discrepancy is apparatus-specific and its scattering conclusion is conditional on that experiment. Use the numerical loss only in the implementation-plan discussion, not as `\(\kappa_F\)`. |
| Marti and Stamper-Kurn (2016) | Define helicity channels first: `\(\theta_F=(\Phi_+-\Phi_-)/2\)` while `\((\Phi_++\Phi_-)/2\)` is a common phase. State the neglected tensor or linear-birefringence terms and the spin moments retained. | An ideal coupling bound is not an apparatus calibration, and a single Faraday observable does not determine a general `\(F>1\)` spin state. |
| Becher et al. (2018) | Define propagation and polarisation geometry before scalar, vector and tensor polarizabilities; separate statistical and systematic uncertainty and identify the dominant source. | Do not put a full Wigner-symbol expansion in the main text unless the numerical calculation uses it. |
| Kaminski et al. (2012) | Motivate dark-field Faraday imaging from high-optical-depth failure, derive its exact inverse under stated assumptions, then report empirical sensitivity and hard limits separately. | Do not transfer its calibrated sensitivity to this species or apparatus. |
| Wigley et al. (2016) | Connect refractive index, phase, Fresnel propagation and the transport-of-intensity equation; state repetitions, averaging and the experimental detection limit. | The reported number of images and lack of detectable heating are specific to its species, detuning and sensitivity. |

### 18.3 Inverse imaging and atom-cloud analysis

| Source | Use and writing pattern | Boundary |
|---|---|---|
| Turner, Domen and Scholten (2005) | Explain why one propagated intensity does not uniquely determine phase, identify contrast-transfer zeros and show the role and cost of Tikhonov regularisation. | Its reconstruction is not evidence that arbitrary clouds or frequencies removed by the pupil are recoverable. |
| Meppelink et al. (2010) | Use the chain destructive-reference limitation, forward model, apparatus, nested inference models and independent validation for quantitative refractive imaging. | Species-specific matrix elements and calibration factors require separate derivation. |
| Altuntas and Spielman (2021) | Specify a Bayesian forward model, likelihood, prior and finite-support assumption, and treat numerical-aperture or transfer-function zeros as an identifiability boundary. | Self-consistency under a prior is not a universal or super-resolution recovery claim. |
| Hofer et al. (2021) | Support region-of-interest detection, segmentation and fit initialisation for atom-cloud images. | Segmentation does not validate an optical inverse or establish recovery of condensate observables. |

For every model equation derived from these sources, use this order when the
items are relevant: measurable quantity; symbols and geometry; assumptions and
controlling regime; specialised equation; parameter provenance; one
representative numerical evaluation; independent check; failure condition.
Experimental-method prose follows purpose, geometry, model-entering parameters,
calibration, data reduction, uncertainty and limitations. A component inventory
is not an experimental argument.

## 19. Final section gate

A section is ready for integration only if all answers are yes:

- Is its question explicit?
- Does each major claim have an equation, result, figure or checked source?
- Are assumptions and validity conditions local?
- Are the reported quantities defined and given with units?
- Is claim strength matched to evidence?
- Are calibration, verification, validation, screening and inference distinct?
- Are citations and bibliography entries closed in both directions?
- Do all comparisons pass the real-comparison test?
- Have the anti-AI review flags been resolved?
- Does the final paragraph answer the section question without generic praise?

## 20. Active-file retention standard

The active repository is a working surface, not an archive. Git history stores
superseded material. A file stays in the active tree only if it passes every
retention gate below.

### 20.1 Retention gates

A retained file must have:

1. **current authority:** its scope and claims agree with the latest approved
   dissertation and hand-off;
2. **unique necessity:** removing it would delete information required for a
   current scientific, operational or reproducibility task;
3. **an active consumer:** a person, script, test, index or current workflow
   uses it;
4. **no unresolved conflict:** it does not preserve a superseded parameter,
   chapter plan, convention or scientific interpretation;
5. **a reproducible role:** commands, inputs, outputs and provenance are
   sufficiently specific to be used;
6. **a maintenance trigger:** the file states what change requires it to be
   updated;
7. **lower maintenance cost than consolidation:** its unique content cannot be
   carried more clearly by an existing authority file.

Failure of one gate is enough to remove the file from the active tree after an
inbound-reference check. "Might be useful" is not a retention reason.

### 20.2 Files that do not belong in the active tree

Remove or consolidate:

- scratch copies, temporary renders and extracted text;
- names containing `old`, `copy`, `final2`, `backup` or an unexplained date;
- completed transition plans and resolved correction reports;
- superseded chapter outlines and alternative hand-offs;
- duplicate writing, figure or notation conventions;
- legacy audits whose only purpose is historical explanation;
- generated artifacts that can be reproduced and are not sealed evidence;
- documents that retain obsolete physics or parameters;
- one-off meeting notes with no active decision or owner;
- empty directories and index entries that point to removed material.

Do not create an `archive/` directory. Git history is the archive.

### 20.3 Protected active categories

Retention is still conditional, but these categories normally have an active
consumer:

- the single current hand-off;
- the single documentation index;
- this writing convention;
- one current source-of-truth file for each physical or software contract;
- current reproducibility instructions;
- sealed evidence required by a claim in the dissertation;
- licences, citation metadata and release manifests that are currently used.

"Sealed evidence" means an identified result with a stable run identifier,
configuration, code provenance and a current dissertation claim. An old plot
or report is not sealed evidence merely because regeneration would take time.

### 20.4 Admission test for a new file

Before adding a file, answer all of the following:

- Why can this content not update an existing authority file?
- Who or what will read it after the current task?
- What event will update or retire it?
- Which index will link to it?
- Which existing file becomes redundant after this addition?

If the last answer is "none", the default action is to update an existing file.
A cleanup task should reduce the active document count. One new authority file
should normally replace at least two overlapping guidance files.

### 20.5 Safe deletion procedure

Before deleting a tracked file:

1. resolve its absolute path within the repository;
2. search all inbound references, scripts and tests;
3. move any still-current unique fact into the correct authority file;
4. verify that a current result or reproduction command does not depend on it;
5. remove the file from the active tree and update indexes in the same change;
6. run link, test or reproduction checks appropriate to the file;
7. report what was removed and note that Git retains its history.

For untracked scratch material created during a task, verify the exact path and
delete it at task completion. Do not leave a second copy of a user document in
the workspace.
