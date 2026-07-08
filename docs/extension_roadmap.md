# Extension Roadmap: Experimental Calibration and Beyond-Thomas-Fermi States

## Purpose

This roadmap describes future simulator extensions for experimental
RAI-data-based calibration and condensate states beyond a simple
Thomas-Fermi BEC.

It is a design document only. It does not implement new physics, change helper
APIs, modify notebook sections, or alter regression baselines.

## A. Current Simulator Boundary

### Current Notebook Behaviour

The original notebook remains the authoritative scientific implementation for
the MSc project. Its current scope is continuous non-destructive imaging of an
ultracold `166Er` Bose-Einstein condensate with notebook-defined parameters,
Thomas-Fermi density estimates, light-atom interaction formulas, imaging
models, camera/noise calculations, multi-shot evolution, and analysis plots.

The notebook's current PCI and DGI paths use a scalar dispersive phase model.
The current Faraday path uses a phenomenological rotation model:

```text
theta_F = kappa_F * phi_peak
kappa_F = 1.0
```

`kappa_F = 1.0` is a placeholder calibration parameter in the current notebook
model. It should not be changed during the notebook-equivalent migration.

### Current Migrated Helper Layer

The current helper layer is a conservative support layer underneath the
notebook. It already includes:

- Atomic Model helpers;
- Light-Atom Interaction helpers;
- shared Fourier propagation helpers;
- shared imaging core helpers;
- tested PCI orchestration;
- tested DGI orchestration;
- camera helper functions.

The frozen lower-layer APIs are:

```text
build_thomas_fermi_state(...)
scalar_phase_shift(...)
scattered_photons_per_atom(...)
faraday_rotation_angle(...)
reabsorption_fraction(...)
```

These APIs should remain stable for Version 1 migration work.

### Current Gaps

The migration is not complete. In particular:

- `simulate_faraday_image(...)` is not yet implemented;
- full camera pipeline migration is not complete;
- shot-noise and multi-shot migration are not complete;
- optimisation should wait until notebook-equivalent simulation layers are
  stable and tested.

### Future Extensions

Experimental calibration and beyond-Thomas-Fermi state support are future
extensions. They should be added after Version 1 notebook-equivalent behaviour
is stable. They should be additive layers above or beside the existing helpers,
not replacements for current notebook physics.

## B. Extension 1: RAI-Data-Based Parameter Calibration

RAI data from the laboratory can be used as a destructive or reference imaging
ground truth for calibrating simulator parameters. This should be treated as an
Analysis / Calibration layer, not as a replacement for notebook physics during
the current migration.

### Candidate Calibratable Parameters

RAI data may help calibrate:

- atom number;
- Thomas-Fermi radii or empirical cloud widths;
- column-density scale;
- imaging magnification;
- effective pixel size;
- camera offset;
- camera gain;
- read-noise level;
- probe intensity;
- detuning offset;
- effective saturation intensity;
- optical-depth calibration;
- reabsorption or residual-absorption scaling;
- future Faraday calibration parameters such as `kappa_F`.

`kappa_F` should only be calibrated in a later Faraday calibration milestone.
It should not be changed during current migration.

### Conservative Calibration Workflow

1. Load RAI image data.

2. Preprocess images:
   - subtract dark frames;
   - apply flat-field or reference correction;
   - select regions of interest;
   - handle pixel size and magnification consistently.

3. Extract observables:
   - atom-number proxy;
   - cloud widths;
   - peak optical density;
   - integrated optical density;
   - line profiles;
   - image residuals;
   - noise statistics.

4. Fit simulator parameters to RAI observables.

5. Store fitted parameters in an external calibration file, not as hard-coded
   constants in the simulator.

6. Use calibrated parameters to run non-destructive PCI, DGI, and Faraday
   simulations.

7. Validate against held-out experimental data by comparing predicted images,
   extracted observables, or information-versus-destruction metrics.

### Suitable Fitting Methods

The first implementation should use simple least-squares fitting against
well-defined observables such as widths, peak OD, integrated OD, and line
profiles.

If the camera noise model is known, a likelihood-based fit can be added later.
This would allow offset, gain, read noise, shot noise, and residual image
variance to be included more coherently.

Bayesian inference is a future extension. It could estimate posterior
uncertainties for atom number, cloud width, detuning offset, probe intensity,
and calibration scale factors. Those uncertainties could then propagate into
imaging optimisation.

### Calibration File Policy

Calibration outputs should be stored in versioned external config files such
as:

```text
calibration/rai_calibration_example.json
calibration/lab_run_YYYYMMDD.toml
```

The exact format can be decided in a later milestone. The key policy is that
calibration values should not silently replace notebook constants or frozen
helper defaults. Code should load them explicitly.

### Boundary Recommendation

RAI calibration should live above the simulator core:

```text
Atomic Model
  -> Light-Atom Interaction
  -> Imaging
  -> Camera
  -> Multi-shot
  -> Analysis / Calibration
```

The calibration layer may call simulator helpers, but lower-level helpers
should not depend on RAI file formats or lab-specific image metadata.

## C. Extension 2: Beyond-Thomas-Fermi State Models

Droplets, supersolids, and mixtures should be represented as new state or
profile generators. They should not be implemented by changing
`build_thomas_fermi_state(...)`.

The Thomas-Fermi helper should remain the Version 1 notebook-equivalent BEC
state model. New states should be additive.

### Possible Future Helpers

Future APIs could include:

```text
build_droplet_state(...)
build_droplet_array_state(...)
build_supersolid_state(...)
build_mixture_state(...)
```

A more flexible Version 2 approach would be a state-provider interface whose
output is a density map, column-density map, phase map, or component-resolved
set of maps. Existing imaging code would then consume the maps without knowing
how they were generated.

### Level 1: Phenomenological Profiles

Level 1 models are simple, deterministic profile generators:

- Gaussian droplets;
- arrays of droplets with configurable spacing and amplitude;
- sinusoidally modulated density profiles for supersolid stripes;
- two-component density maps for mixtures;
- optional empirical envelope functions.

This level is likely the most realistic first beyond-TF extension for MSc
scope because it tests imaging consequences without requiring a full many-body
solver.

### Level 2: Semi-Physical Profiles

Level 2 models can include qualitative physics-inspired terms:

- contact-interaction balance;
- dipolar shape anisotropy;
- Lee-Huang-Yang-inspired density stabilisation;
- empirical droplet size and spacing parameters constrained by lab data.

This level is still not a full Gross-Pitaevskii solver. It should be documented
as semi-phenomenological unless validated against separate theory or experiment.

### Level 3: External Numerical Input

Level 3 avoids implementing state physics in this repository. The simulator
loads density maps from:

- experimental reconstruction;
- separate Gross-Pitaevskii simulations;
- lab-provided image analysis;
- archived numerical datasets.

The simulator then handles optical imaging, camera modelling, multi-shot
effects, and information-versus-destruction analysis.

For MSc scope, Level 1 and Level 3 are likely more practical than building a
full droplet or supersolid solver inside this repository.

## D. Mixture Extension

Mixtures require a component-aware architecture. A component may represent a
different species, isotope, spin state, or internal-state population.

Mixture support may require:

- multiple density or column-density profiles;
- different masses and transition parameters;
- different scalar phase shifts;
- different Faraday rotations;
- different scattering rates;
- overlapping spatial profiles;
- component-resolved imaging observables;
- total signal calculated as a scalar sum, coherent field sum, or vector
  polarisation combination depending on the imaging method.

A future data structure could look conceptually like:

```text
MixtureState
  components:
    - label
    - density_map or column_density_map
    - transition_parameters
    - calibration_parameters
```

This should not be implemented until the single-component notebook-equivalent
simulator is stable. When it is implemented, it should be additive and should
not change `build_thomas_fermi_state(...)`.

## E. Relation To Information Versus Destruction

Both extension directions support the central project goal: estimating how much
information is gained per unit destructiveness in continuous non-destructive
imaging.

RAI calibration can anchor the simulator against destructive ground-truth
images. This helps estimate atom number, cloud size, OD scale, camera response,
and residual absorption more realistically.

Beyond-TF states make the simulator relevant to realistic lab targets such as
droplets, supersolids, and mixtures. These systems may have sharper structures,
multiple density peaks, or component-dependent optical responses, so the
information gained per scattered photon may differ from a smooth TF cloud.

The optimisation target should remain explicit:

```text
information gained per scattered photon / destructiveness
```

Future optimisation should compare PCI, DGI, Faraday, and camera/noise choices
using calibrated or experimentally relevant states, but only after the
notebook-equivalent simulator is stable.

## F. Proposed Future Milestone Sequence

### Current Migration Path

1. Complete Faraday orchestration helper against the Faraday baseline.
2. Complete Camera layer migration.
3. Complete shot-noise migration.
4. Complete multi-shot migration.
5. Begin optimisation only after the notebook-equivalent path is validated.

### Future Extension Path

6. Add RAI calibration design details and an example data-loader skeleton.
7. Add a calibration config schema with no hard-coded lab constants.
8. Add deterministic synthetic RAI calibration tests.
9. Add phenomenological droplet and supersolid profile generators.
10. Add external density-map input support.
11. Add mixture-state abstraction.
12. Add calibrated Faraday analysis and optimisation, including a controlled
    treatment of `kappa_F`.

The extension path should not interrupt the current migration path.

## G. Risks And Boundaries

### Risks

- Overfitting a small amount of RAI data.
- Confusing destructive RAI calibration with non-destructive imaging
  simulation.
- Introducing uncontrolled physics before notebook-equivalent behaviour is
  stable.
- Changing `kappa_F` prematurely.
- Implementing a full many-body droplet or supersolid solver beyond MSc scope.
- Adding broad abstractions before the existing notebook behaviour is migrated
  and tested.
- Mixing lab-specific file formats into low-level physics helpers.

### Recommendations

- Keep Version 1 notebook-equivalent.
- Add calibration only after migration layers are stable.
- Keep calibration parameters in external config files.
- Keep new state models additive rather than replacing Thomas-Fermi helpers.
- Keep RAI calibration in Analysis / Calibration, not in Atomic Model or
  Imaging.
- Treat `kappa_F` calibration as a future Faraday-specific calibration stage.
- Prefer Level 1 phenomenological profiles or Level 3 external density maps
  before attempting semi-physical or solver-level droplet/supersolid modelling.

## Summary

The current simulator should remain focused on conservative notebook-equivalent
migration. RAI calibration and beyond-Thomas-Fermi state support are valuable
future extensions, but they should be implemented only after Faraday, camera,
shot-noise, and multi-shot migration are stable.

The safest long-term design is additive:

```text
current TF notebook-equivalent simulator
  + calibration configs and analysis layer
  + optional state/profile providers
  + optional mixture abstraction
  + later calibrated optimisation
```

This keeps the scientific baseline intact while leaving a clear route toward
experimental calibration and realistic lab targets.
