# Current Codex hand-off

- **Updated:** 22 July 2026
- **Status:** single current hand-off for the dissertation and repository
- **Branch:** `codex/reconstruction-v4-canonical`
- **Clean source/config snapshot:** `498166f5d712c0a1f035fb7dd1cacd62619dc62d`
- **Result/provenance layer:** the current branch tip containing this hand-off
- **Worktree after checkpoint:** clean
- **Publication status:** local only; no push performed
- **Update trigger:** an approved scientific-contract change, a completed
  implementation checkpoint, a branch change or a change to the immediate task

This file records current decisions and the next task. Future superseded
versions will remain in the publication history; pre-publication hand-offs and
transition notes remain only in the local safety backups.

## 1. Operating rules

1. Run `git status --short --branch` before changing anything.
2. Continue on `codex/reconstruction-v4-canonical` unless the user requests a
   different branch.
3. Do not modify, merge or check out `main` without explicit permission.
4. Preserve unrelated user changes. Do not copy files from an older worktree
   merely because they have a later timestamp.
5. Use the active configs, tested code and indexed evidence in this repository
   as numerical authority. The latest dissertation prose remains user-owned and
   is supplied separately when textual integration is needed.
6. Put temporary renders and extracted documents under ignored `.scratch/`
   paths and remove them when the task is complete.
7. Do not create archive directories or new transition reports. Update this
   hand-off; Git retains the previous state.

## 2. Authority order

When sources disagree, use:

1. the user's latest explicit decision;
2. this hand-off;
3. active configs and tested code on this branch;
4. the indexed model and evidence authorities;
5. stored result metadata;
6. historical notebook regression material.

Current authorities are:

- `docs/dissertation_writing_conventions.md` for dissertation prose, figures,
  citations and active-file retention;
- `docs/simulation_reference_parameters.md` for the optical, detector,
  sampling and implemented-sequence contract;
- `docs/multiframe_heating_model_optimisation.md` for the approved replacement
  of the Version 1 Chapter 5 thermodynamics;
- `docs/reconstruction_architecture.md` for the inverse and observable contract;
- `docs/reconstruction_orca_v4_evidence_2026_07_21.md` for sealed historical
  inverse evidence;
- `docs/figure_index.md` for result status and provenance;
- `docs/reproducibility.md` for commands and environment checks;
- `docs/experimental_measurement_plan.md` for the optional laboratory plan.

The original notebook and notebook-aligned exports are historical regression
references. They are not current physical or software authority.

## 3. Dissertation direction

The dissertation connects four tasks:

1. identify physical observables for which repeated non-destructive imaging is
   useful, with particular attention to dipolar `166Er`;
2. propagate a declared condensate state through atom-light response, finite
   aperture, detector noise and repeated-exposure evolution;
3. screen acquisition conditions and estimate supported observables, with
   uncertainty, from synthetic raw camera channels;
4. translate the resulting requirements into an initial dual-port Faraday
   implementation plan.

The dissertation does not depend on experimental Faraday data. Synthetic
observable inference is a result within the declared forward-model contract;
the laboratory chapter remains a design and commissioning plan.

The chapter structure is:

1. Introduction;
2. Dipolar condensates and target observables;
3. Non-destructive imaging physics;
4. Forward model;
5. Acquisition screening;
6. Physics-informed inference;
7. Observable estimates and recoverability;
8. Initial DPFI implementation plan;
9. Conclusions and outlook.

Chapter 5 uses a declared camera-space `SNR_5x5` screen. Chapters 6 and 7 decide
whether the raw channels support the requested physical observables. SNR is not
the final definition of recoverability.

## 4. Active optical and detector contract

The compact optical benchmark uses the Oxford condensate-core surrogate with
`N0=2.5e4`. It is not the thermodynamic initial state for the approved
multiframe replacement.

The active Chapter 4 and single-frame Chapter 5 inputs are:

| Quantity | Active value | Status |
| --- | ---: | --- |
| Species | `166Er` | reference species |
| Reference detuning | `|Delta|/2pi = 1.5 GHz` | scan reference |
| Fluence range | `90-300 mW us` | screening interval |
| Representative fluences | `90, 150, 300 mW us` | A, B and C |
| Effective numerical aperture | `0.130` | optical-design input |
| Camera | ORCA-Fusion C14440-20UP | selected detector |
| Magnification | `10` | design input |
| Object-plane pixel | `0.650 um` | derived from `6.5 um / M` |
| Quantum efficiency | `0.65` | manufacturer-typical input |
| Read noise | `0.7 e- rms` | manufacturer-typical Ultra quiet scenario |
| Faraday conversion | `kappa_F=-45/91` | signed ideal atomic estimate |
| Effective apparatus response | `TBD` | requires calibration |

The maintained `29.5 MHz` linewidth is a reproducible numerical input whose
modern primary-source provenance remains unresolved. It must not be presented
as a closed literature value until that source is identified or the input is
changed through the normal config, figure and test update.

## 5. Repeated-exposure status

The stored Figures 5.2 and 5.4 and their `N_dep`/`N_use` values use the Version 1
ideal saturated-gas heating model. They are **frozen Version 1 screening outputs
pending the approved heating replacement**. They remain reproducible but are not
current quantitative predictions of an Oxford sequence.

The approved replacement is a quasi-equilibrium, recoil-limited,
fixed-trapped-number screening model defined in
`docs/multiframe_heating_model_optimisation.md`. It will:

- initialise from the three direct Oxford 300 ms `(N0, Nth, T)` triplets;
- use the measured-scale anchored thermal non-saturation closure;
- run the two declared energy coefficients `c_E=2.70118...` and `c_E=3` as a
  sensitivity span;
- use signed `kappa_F=-45/91` for signal and `chi_sc=46/91` exactly once for the
  current weak, axial, linearly polarised Faraday scattering response;
- report condensate depletion separately from trapped-atom loss;
- support only the declared full-resolution Ultra quiet cadence, conditional on
  rethermalisation.

No Figure 5.2 or 5.4 regeneration is accepted until one shared state engine
passes the specified tests and its diagnostic trajectories are reviewed. Both
figures are excluded from the default run, and their Version 1 generators refuse
to overwrite the canonical directories.

The retained Figure 5.2 and Figure 5.4 numerical and rendered artifacts were not
regenerated during clean-history consolidation. Their metadata preserves the
original dirty generation base `df00b42ae509921e663cccc0adc0a7a7a240c2a8` as a
pre-publication local checkpoint identifier and records `498166f` only as the
retrospective source/config snapshot. This metadata-only correction is not a
clean-generation claim. Both result records declare
`status=frozen_version_1_screening_output`, `current_oxford_prediction=false`
and `canonical_regeneration_allowed=false`.

## 6. Reconstruction and observable contract

DPFI is not claimed to recover a unique high-resolution density image after the
pupil has removed spatial information. RAI remains the terminal reference for a
detailed density distribution.

The primary inverse output is

```text
q = (A, y_c, z_c, w_major),
```

where `A` is the supported integrated response, `(y_c,z_c)` is the centroid and
`w_major` is the major-axis rms width derived from the complete second-moment
tensor. Each component carries an uncertainty, a data-consistency interval or
an unsupported status. The fitted density field is an internal nuisance object,
not a recovered image interpreted pixel by pixel.

The active code implements the signed Faraday response. The sealed v4 inverse
evidence was generated earlier under `NA=0.080`, `kappa_F=1` and
`sigma_r=1.4 e- rms`. It remains method-development evidence and cannot support
current `166Er` amplitude, SNR or usable-sequence claims.

The planned seven-condition initial-condition suite is a synthetic
method-development replay at `F=90,150,300 mW us`. No accepted result exists.
Its source inverse remains the sealed historical `kappa_F=1` contract, so it
cannot establish current `166Er` performance. Broad Oxford thermal states exceed
the current camera/support envelope, and the Er stress states test resolution
or monotonicity failure rather than ordinary quantitative recovery. Preserve
the preliminary source and tests, but do not present the suite as a dissertation
result until its sample subset and supports are approved and a new result is
sealed.

## 7. Repository and result boundaries

Maintained code and contracts live under:

```text
src/non_destructive_image/                 forward model
src/non_destructive_image/reconstruction/  inverse and observables
configs/                                   active and historical contracts
scripts/                                   generators, orchestration and validation
tests/                                     unit and regression suite
docs/                                      active documentation authorities
```

Current and historical result status is defined only in `docs/figure_index.md`.
The two large NPZ files under `regression/baseline/imaging/` are intentional
regression baselines. A generated file is not retained merely because it is
expensive to regenerate.

The prepared clean history starts at `origin/main`, followed by the consolidated
source/config commit `498166f` and the current result/provenance commit. No
pre-rewrite feature commit is retained in the publication DAG. Historical
identifiers including `5e59406`, `f03802b`, `64a07f7`, `cda577a`, `e7e8b1d` and
`df00b42` are pre-publication local checkpoint IDs and are not expected to
resolve on the eventual public branch. That DAG remains recoverable only through
the local safety refs and verified offline bundle; none may be pushed.

## 8. Immediate task

1. Implement only the minimum heating replacement
   specified in `docs/multiframe_heating_model_optimisation.md`.
2. Produce and review the diagnostic state table before regenerating Figures
   5.2 and 5.4 or changing dissertation numbers.
3. Regenerate the frozen figures only after the shared state engine and its
   diagnostic trajectories have been accepted.

Do not run the unapproved initial-condition suite and do not regenerate the
long sealed inverse studies during documentation or cleanup work.

## 9. Reproduction and verification

Use the repository-local Python 3.12 environment and the exact commands in
`docs/reproducibility.md`. Bare `python`, a globally resolved `pytest` and other
editable clones are not accepted substitutes.

The 22 July 2026 clean-history checkpoint used the repository-local Python 3.12
environment and passed:

- `290` tests in `52.36 s` from a fresh checkout of the staged tree;
- parsing of every retained JSON file;
- local-link validation of the active Markdown files;
- `pip check`, package-origin verification and the default generation dry run;
- validation of the planned seven-condition, `90/150/300 mW us` suite without
  running its fits;
- active Figure 4.2, Figure 5.1 and detuning metadata attribution to `498166f`;
- byte equality of all 15 frozen Figure 5.2/5.4 non-metadata artifacts;
- `git diff --check`, retired-entry scanning, public-history privacy scanning and
  the observable artifact hash/run-ID chain.

The branch is `codex/reconstruction-v4-canonical`, `0` behind and `2` commits
ahead of `origin/main`, with a clean worktree. Two local publication commits were
created. No push, merge or change to `main` occurred; `main` and `origin/main`
remain at `d6dd55855e422779cef3ed30f3aec37df26e6d35`.

## 10. Open inputs

- The dissertation source remains external and read-only until the user
  supplies it for a textual integration pass.
- `NA=0.130`, `QE=0.65` and `0.7 e- rms` remain screening inputs until measured.
- The effective apparatus response and absolute Faraday density scale require
  paired calibration.
- Trap depth, beam coverage and inter-frame equilibration remain conditions on
  the approved fixed-number heating replacement.
- The clean feature branch awaits explicit push approval. The local safety ref
  and offline bundle contain the pre-rewrite DAG and must not be published.
- A repository licence has not yet been selected; public reuse rights must not
  be implied before one is added.
