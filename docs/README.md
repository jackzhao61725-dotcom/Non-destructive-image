# Documentation index

This directory contains the active documentation for the dissertation forward
model and its experimental hand-off. Date-stamped transition notes and abandoned
simulation branches have been removed from the active surface; their history
remains available through Git or the 20 July 2026 cleanup archive.

## Start here

- [Simulation reference parameters](simulation_reference_parameters.md) —
  current provisional optical, detector and sequence contract.
- [Figure index](figure_index.md) — figure-to-script-to-config provenance.
- [Experimental measurement plan](experimental_measurement_plan.md) —
  laboratory tasks, priorities, recorded data and acceptance criteria.
- [Reproducibility](reproducibility.md) — validation and regeneration commands.
- [Repository cleanup report](repository_cleanup_2026_07_20.md) — archived
  material, contract changes and validation evidence for this transition.

## Model

- [Architecture](architecture.md) — maintained layers and data flow.
- [Physics model](physics_model.md) — Version 1 equations and approximations.
- [Faraday model](faraday_model.md) — polarisation readouts and the
  uncalibrated `kappa_F` boundary.
- [Baseline specification](baseline_specification.md) — historical numerical
  baselines retained for regression.
- [Colab validation](colab_validation.md) — portable validation notes.

## Results and audits

- [Results guide](results_readme.md) — active result directories and focused
  regeneration commands.
- [Canonical performance gate](performance_validation_v1_report.md) — current
  10-frame, 228-pixel matched-ROI evidence.
- [Thesis numerical consistency audit](thesis_numerical_consistency_correction_report.md)
  — retained legacy corrections and current contract checks.
- [Linear-approximation audit](linear_approximation_validity_audit.md) —
  validity of weak-phase and small-rotation interpretations.

## Figure conventions

- [Figure language](figure_language_conventions.md)
- [Quantity labels](figure_quantity_label_conventions.md)

## Publication

- [Zenodo release checklist](zenodo_release_checklist.md)

All current Faraday outputs remain structural or screening comparisons while
`kappa_F=1` is uncalibrated. Provisional detector and optical values must be
replaced by apparatus measurements before an experimentally calibrated result
is claimed.
