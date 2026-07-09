# Camera Noise Milestone Report

## Objective

Add a thin stochastic camera helper around the existing camera noise recipe
with explicit random-number-generator handling.

## Helper Added

`simulate_noisy_camera_image(...)` was added to
`src/non_destructive_image/camera.py` and exported from
`src/non_destructive_image/__init__.py`.

The helper only orchestrates existing camera helpers:

```text
binned_image = bin_to_camera_pixels(image, bin_size)
noisy_counts = add_camera_noise(binned_image, photons_per_pixel, rng, read_noise_electrons)
noisy_image = normalize_camera_counts(noisy_counts, photons_per_pixel)
```

If `input_is_binned=True`, the helper uses the supplied image as an already
binned camera image. If `normalize=False`, it returns noisy electron counts
instead of normalised image units.

## RNG Policy

The caller must pass an explicit `np.random.Generator`.

The helper does not create a hidden global RNG.

The helper does not hard-code a seed.

## Tests Added

`tests/regression/test_camera_noise_pipeline.py` checks:

- same seed gives identical output;
- different seed gives different stochastic output;
- output shape is correct;
- output values are finite;
- helper output agrees with direct composition of existing helpers;
- already-binned input can return noisy counts directly.

## Scope Confirmation

Validation results:

```text
pytest -q: 32 passed
python scripts\validate_notebook_sections.py: passed
```

No notebook sections were changed.

No imaging helpers were changed.

No Atomic or Light-Atom helpers were changed.

No baseline `.npz` files were changed.

No PCI, DGI, or Faraday code was changed.

No multi-shot or optimisation code was changed.

No deliverable zip was regenerated.
