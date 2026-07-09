# Camera Pipeline Milestone Report

## Objective

Migrate only the deterministic camera pipeline around the existing camera
helpers.

## Helper Added

`simulate_camera_image(...)` was added to `src/non_destructive_image/camera.py`
and exported from `src/non_destructive_image/__init__.py`.

The helper performs only deterministic processing:

```text
camera_image = bin_to_camera_pixels(image, bin_size)
```

If `photons_per_pixel` is supplied, it also applies the deterministic count
conversion and normalisation around the existing helper:

```text
deterministic_counts = binned_image * photons_per_pixel
camera_image = normalize_camera_counts(deterministic_counts, photons_per_pixel)
```

No stochastic Poisson noise or Gaussian read noise is added.

## Tests Added

`tests/regression/test_camera_pipeline.py` checks:

- output shape;
- finite output values;
- consistency with `bin_to_camera_pixels(...)`;
- consistency with `normalize_camera_counts(...)`;
- deterministic normalisation behaviour.

## Scope Confirmation

Validation results:

```text
pytest -q: 27 passed
python scripts\validate_notebook_sections.py: passed
```

No notebook sections were changed.

No imaging helpers were changed.

No Atomic or Light-Atom helpers were changed.

No baseline `.npz` files were changed.

No Faraday, PCI, or DGI code was changed.

No shot-noise, multi-shot, or optimisation code was changed.

No deliverable zip was regenerated.
