# Faraday Model Audit

This document describes exactly what the notebook implements for Faraday imaging. It is documentation only.

## Direct answer

The notebook implements option **B**:

> It defines `theta_F` phenomenologically before Jones propagation.

The notebook does not derive `theta_F_peak` from independently computed `phi_+` and `phi_-` values. Instead, it sets `kappa_F = 1.0` as a placeholder and returns `kappa_F * phi_peak(...)`.

## Implementation evidence

The Faraday section defines:

```python
kappa_F = 1.0

def theta_F_peak(Delta_Hz, n_col_peak):
    return kappa_F * phi_peak(Delta_Hz, n_col_peak)
```

The surrounding comment describes `kappa_F` as a placeholder vector-polarizability fraction and `kappa_F=1` as an idealised maximal-coupling bound.

The Jones propagation then receives `theta_F_val` and applies circular phases `+theta_F_val*prof` and `-theta_F_val*prof`:

```python
Pp = 1 + ifft2(fft2(exp(+i theta_F_val*prof)-1)*pupil)
Pm = 1 + ifft2(fft2(exp(-i theta_F_val*prof)-1)*pupil)
Ex, Ey = (Pp + Pm)/2, 1j*(Pp - Pm)/2
```

Thus, in notebook convention, the quantity returned by `theta_F_peak` is already the rotation angle used by Jones propagation.

## Current meaning of `kappa_F`

`kappa_F` currently represents a dimensionless calibration factor mapping the scalar phase returned by `phi_peak` to the Faraday rotation angle used in the notebook:

```text
theta_F_peak = kappa_F * phi_peak
```

## Where `kappa_F` enters

`kappa_F` enters only inside `theta_F_peak`. Downstream Faraday functions receive `theta_F` or `theta_F_val`; they do not multiply by `kappa_F` again.

## Assumptions required by current implementation

- `phi_peak` is the scalar dispersive phase from the scalar light-atom interaction section.
- `theta_F_peak` is the physical Faraday rotation angle used by the Jones model.
- The small common-mode scalar phase is neglected in the Faraday pure-rotation model.
- `kappa_F` absorbs the missing atomic-structure calibration.

## What `kappa_F` does not currently represent

- It is not derived in the notebook from Clebsch-Gordan coefficients.
- It is not derived in the notebook from vector polarizability.
- It is not applied after Jones propagation.
- It is not a separately computed circular-component phase.
- It is not a measured experimental calibration.

## Approximation and assumption audit

| Approximation / assumption | Reason given or implied in notebook | Where used | Possible impact |
|---|---|---|---|
| Thomas-Fermi condensate | Notebook models condensate using TF chemical potential, peak density, radii, and column density. | Atomic model section. | Sets all density and phase scales. |
| Peak column-density model | Peak column densities are computed as `(4/3)*n_peak*R`. | Atomic model and phase calculations. | Imaging signals use peak values and TF spatial profile. |
| Scalar dispersive phase | Absorption is treated as negligible for far detuned probe in scalar phase section. | `phi_peak`, PCI, DGI. | Ignores absorption in core phase model except residual OD diagnostics. |
| Phenomenological `kappa_F` | Atomic-structure calibration is pending. | `theta_F_peak`. | Faraday amplitude is a placeholder, not a calibrated prediction. |
| Pure Faraday rotation | Common-mode scalar phase is neglected in Faraday Jones model. | `sim_faraday_fields`. | Faraday model isolates polarization rotation only. |
| Circular components have phases `±theta_F` | Implements rotation angle directly. | `sim_faraday_fields`. | Factor-of-two is embedded in the definition of `theta_F`. |
| Finite numerical aperture | Imaging arm has finite NA pupil. | `sim_image`, `sim_faraday_fields`. | Blurs and dilutes contrast. |
| Fourier-plane carrier manipulation | PCI and DGI act on carrier differently. | `sim_image`. | Distinguishes PCI phase plate from DGI dark stop. |
| Camera binning | High-resolution grid is averaged into camera pixels. | `to_camera`, `camera_from_image`. | Reduces spatial resolution and changes noise level per pixel. |
| Poisson shot noise plus Gaussian read noise | Camera model includes photon counting and read noise. | `to_camera`, `camera_from_image`. | Stochastic images depend on fixed RNG seed. |
| Heating and clean-loss bracketing | Multi-shot section brackets destruction models. | `run_sequence`, `Nmax_heating`, `Nmax_cleanloss`. | Affects frame count, condensate depletion, and SNR evolution. |

## Self review and remaining ambiguity

- The notebook is internally explicit about how `theta_F` is used in Jones propagation, but the physical derivation of `kappa_F` is not present.
- The relationship between scalar `phi_peak` and microscopic circular phases remains the main unresolved definition.
- Future documentation should include a derived or referenced Er vector-polarizability calibration before replacing `kappa_F`.
