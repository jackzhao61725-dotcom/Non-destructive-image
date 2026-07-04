# Milestone Report

## Objective

Document the physical model implemented by the notebook, especially the definitions of scalar phase, Faraday rotation, Jones fields, camera images, and the placeholder calibration parameter `kappa_F`.

## Work Completed

- Added `docs/physics_model.md` documenting major notebook quantities and pipelines.
- Added `docs/faraday_model.md` documenting the notebook Faraday model, `kappa_F`, and the model assumptions.
- Documented scalar phase, Faraday rotation, PCI, DGI, and dual-port Faraday pipelines.
- Added an approximation and assumption audit table.
- Explicitly documented that the notebook defines `theta_F` phenomenologically before Jones propagation.

## Files Modified

No existing scientific implementation files were modified.

## New Files

- `docs/physics_model.md`
- `docs/faraday_model.md`
- `docs/milestone_4_report.md`

## Validation Performed

PASS — documentation review against exported notebook sections.

- Checked the definitions against `notebook_sections/02_atomic_model.py`, `notebook_sections/03_light_atom_interaction.py`, `notebook_sections/04_pci.py`, `notebook_sections/05_dgi.py`, `notebook_sections/06_faraday.py`, `notebook_sections/07_camera.py`, and `notebook_sections/09_multishot_simulation.py`.

PASS — no code changes.

- No notebook, helper, script, or test code was modified.

## Behaviour Changes

No observable behaviour changed.

- No equations were modified.
- No algorithms were modified.
- No parameter values were modified.
- No notebook code was modified.
- This milestone is documentation only.

## Remaining Issues

- `kappa_F` remains a placeholder rather than a physical calibration.
- The microscopic relationship between scalar `phi_peak` and circular-component phases remains unresolved.
- PCI, DGI, Faraday, camera, and multi-shot numerical baselines are still required before migration.

## Recommended Next Milestone

Scientific Baseline Generation in a complete scientific Python environment.

## Questions Requiring Review

None.

## Git Summary

- Current branch: `work`
- Latest commit hash at report creation: `ae6425b`
- Latest commit message at report creation: `Add scientific baseline infrastructure`
- PR title: `Add scientific baseline infrastructure`
- PR status: Created for review

## Executive Summary

Milestone 4 made the notebook's physical model explicit without changing any code. The new documentation defines density, column density, scalar phase, phase maps, Faraday rotation, Jones fields, and camera images; traces scalar, PCI, DGI, and Faraday pipelines; and records the current Faraday model as phenomenological in `theta_F_peak`. The main unresolved physics issue remains calibration of `kappa_F` from a proper circular-component or vector-polarizability derivation.

## Review Request

Please review whether the documented quantity definitions match the notebook's intended scientific meaning. Pay particular attention to the distinction between scalar phase, circular-component phase, Faraday rotation angle, and the placeholder calibration `kappa_F`. No migration or code changes should begin until this physical documentation is accepted.

# Physics Review

The notebook implements a scalar dispersive phase model for PCI/DGI and a phenomenological Faraday rotation model for Faraday imaging. The scalar phase is `phi_peak`; the Faraday angle is `theta_F_peak = kappa_F * phi_peak`; Jones propagation then applies circular phases `±theta_F`. The assumptions that were previously implicit are now documented: `kappa_F` is a placeholder, common-mode scalar phase is neglected in Faraday propagation, and Faraday amplitude is not yet derived from Er vector polarizability. Ambiguity remains in the physical calibration connecting scalar `phi_peak` to microscopic circular-component phases.
