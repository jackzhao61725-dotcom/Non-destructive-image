# Documentation index

- **Status:** single index of active repository documentation
- **Update trigger:** a document is added, removed, renamed, superseded or
  changes scientific status

This directory is the active documentation surface. Superseded specifications,
audits, conventions and transition reports are retained in Git history, not in
an archive directory.

## Start here

- [Current Codex hand-off](../CODEX_HANDOFF_CURRENT.md) — dissertation
  direction, scientific boundaries, worktree state and immediate task.
- [Dissertation writing conventions](dissertation_writing_conventions.md) —
  single authority for argument, equations, figures, source use, prose review
  and file-retention decisions.
- [Simulation reference parameters](simulation_reference_parameters.md) —
  current condensate, detector, sampling, sequence and Faraday-response
  contract, including frozen-output boundaries.
- [Reproducibility](reproducibility.md) — verified local interpreter, tests,
  regeneration commands and output-retention rules.

## Model and evidence

- [Multiframe heating model optimisation basis](multiframe_heating_model_optimisation.md)
  — approved implementation specification for the direct Oxford initial
  states, non-saturation closure, recoil-limited energy update, validity gates
  and replacement tests for the frozen Version 1 Chapter 5 sequence model.
- [Reconstruction architecture](reconstruction_architecture.md) — raw-channel
  forward operator, latent nuisance field, observable vector, support,
  uncertainty and implementation boundaries.
- [Figure and data index](figure_index.md) — output-to-generator-to-config
  provenance and the current or frozen status of each result family.
- [ORCA reconstruction evidence](reconstruction_orca_v4_evidence_2026_07_21.md)
  — sealed synthetic method-development record. It preserves the historical
  `kappa_F=1` operator and is not a current quantitative `166Er` Faraday
  performance prediction.

## Experimental implementation

- [Experimental measurement plan](experimental_measurement_plan.md) — optional
  HWP–Wollaston dual-port commissioning, calibration, observable-specific
  performance and held-out assessment plan. It is not a dissertation completion
  requirement.

## Admission rule

A new document is admitted only when an existing authority cannot carry the
content, a named consumer will use it, its update or retirement event is known,
and the documentation index will link to it. The complete retention test is in
Section 20 of the writing conventions.
