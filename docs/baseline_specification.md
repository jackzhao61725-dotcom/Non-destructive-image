# Scientific Baseline Specification

## Purpose

The scientific baseline will provide numerical reference outputs from the authoritative notebook. Future refactors must compare against these files before notebook-local code is replaced by helper modules.

This specification defines the baseline layout and validation policy only. It does not contain generated baseline data.

## Directory layout

```text
regression/
  baseline/
    README.md
    metadata.json                 # generated later
    notebook_outputs.json          # existing stored-output baseline
    atomic/
      tf_state.npz                 # generated later
      tf_profile.npz               # generated later
      column_density.npz           # generated later
      phase_map.npz                # generated later
      scattering.npz               # generated later
      recoil.npz                   # generated later
    imaging/
      pci_image.npz                # generated later
      dgi_image.npz                # generated later
      faraday_image.npz            # generated later
      ideal_camera_image.npz       # generated later
      shot_noise_image.npz         # generated later
    multishot/
      evolution.npz                # generated later
      shot000.npz                  # generated later
      shot010.npz                  # generated later
      shot020.npz                  # generated later
```

## Required metadata

`metadata.json` should be generated together with the numerical baseline and should contain:

- notebook filename and notebook SHA-256 hash,
- generation timestamp in UTC,
- Python version,
- package versions for `numpy`, `scipy`, `matplotlib`, notebook execution tools, and any other runtime dependencies,
- random seed policy used for stochastic outputs,
- list of generated baseline files and their SHA-256 hashes,
- notes on any intentionally skipped quantities.

## Scientific quantities to save

### Atomic and light-atom quantities

- Thomas-Fermi radii in metres, shape `(3,)`.
- Trap angular frequencies in rad/s, shape `(3,)`.
- Chemical potential in joules, scalar.
- Chemical-potential temperature in kelvin, scalar.
- Peak density in m^-3, scalar.
- Column density along principal axes in m^-2, shape `(3,)`.
- Representative 2D Thomas-Fermi profile arrays on the notebook grid, shape matching the notebook grid.
- Representative optical phase map on the notebook grid, shape matching the notebook grid.
- Scattering quantities, including scattered photons per atom per shot at reference operating points.
- Recoil energy, recoil temperature, and recoil velocity.

### Imaging quantities

- PCI noiseless image arrays.
- DGI noiseless image arrays.
- Faraday dark-field and dual-port image arrays where available.
- Ideal camera-binned images before stochastic noise.
- Camera images with stochastic shot/read noise using the fixed regression seed.

### Multi-shot quantities

- Multi-shot evolution arrays: shot index, condensate number, loss fraction, temperature, phase, and SNR.
- Selected per-shot image states, for example shots `000`, `010`, and `020` when available.

## File format

Use compressed NumPy archives: `.npz` produced by `numpy.savez_compressed`.

Each `.npz` file should include:

- one or more named arrays,
- a `description` string array if useful,
- units encoded in array names or metadata where practical,
- no pickled Python objects.

## Naming conventions

- Use lowercase snake-case file names.
- Use explicit scientific names rather than implementation-local variable names where possible.
- Include units in array names when ambiguity is likely, for example `radii_m`, `column_density_m2`, `chemical_potential_j`.
- Preserve notebook reference operating points in metadata rather than file names unless multiple operating points are saved.

## Random seed policy

The notebook currently uses a fixed random seed for reproducible noise. Scientific regression baselines must preserve this policy.

- Deterministic quantities must be generated without stochastic sampling.
- Stochastic quantities must be generated with a fixed, documented seed.
- Store the random seed and generator type in `metadata.json`.
- Do not compare independent stochastic draws against each other. Compare either the saved deterministic pre-noise image or a saved noisy image generated with the fixed seed.

## Tolerance policy

Default tolerances should be conservative:

- Scalars from algebraic helper extractions: `rtol=1e-12`, `atol=0` where practical.
- FFT/image arrays: `rtol=1e-10`, `atol=1e-12` unless numerical backend differences require documented relaxation.
- Stochastic fixed-seed images: exact equality may be appropriate only if the same NumPy version and generator are used; otherwise compare saved noisy arrays by hash in locked environments and compare deterministic pre-noise arrays numerically across environments.
- Display figures should not be the primary regression target; compare the underlying arrays instead.

## Generation rule

Baseline files must be generated from the original notebook/reference implementation before helper modules are wired into notebook sections. Do not generate baselines from already-refactored helper code.
