# Physics Model Implemented by the Notebook

This document records the physical meanings used by the current notebook. It is descriptive only; it does not change equations, algorithms, or parameters.

## Quantity definitions

| Quantity | Physical meaning | Mathematical definition in notebook | Computed in | First used in |
|---|---|---|---|---|
| Density | Thomas-Fermi condensate density represented through peak density and radii. | `n_peak = mu*m/(4*np.pi*hbar**2*a_s)` with `R = sqrt(2*mu/(m*omega**2))`. | `notebook_sections/02_atomic_model.py` lines 20-28. | Used to form peak column density `n_col`. |
| Column density | Peak column density along each principal axis. | `n_col = (4/3)*n_peak*R`. | `notebook_sections/02_atomic_model.py` lines 25-28. | Passed to `phi_peak(Delta_Hz, n_col[axis])`. |
| `phi_peak` | Scalar dispersive column phase for PCI/DGI-style imaging. | `sigma0 * n_col_peak * d / (2*(1 + d**2))`, with `d = delta_of(Delta_Hz)`. | `notebook_sections/03_light_atom_interaction.py` lines 25-31. | Used in phase tables and image simulations. |
| Phase map | Spatially varying scalar phase over the image plane. | `phi_peak_val * prof`, where `prof = _tf_profile(...)`. | `notebook_sections/04_pci.py` lines 100-106. | Used inside `exp(1j*phi_peak_val*prof)-1`. |
| `theta_F_peak` | Peak Faraday rotation angle used by the notebook's Faraday model. | `kappa_F * phi_peak(Delta_Hz, n_col_peak)`. | `notebook_sections/06_faraday.py` lines 33-47. | Passed to `sim_faraday_fields`. |
| Faraday rotation map | Spatially varying polarization rotation over the image plane. | `theta_F_val * prof`. | Implied by `sim_faraday_fields` using `exp(±1j*theta_F_val*prof)`. | Used to form circular-component fields. |
| Jones field | NA-limited recombined linear electric field after circular components acquire opposite phases. | `Pp = 1 + ifft2(fft2(exp(+i theta_F prof)-1)*pupil)`, `Pm = 1 + ifft2(fft2(exp(-i theta_F prof)-1)*pupil)`, then `Ex=(Pp+Pm)/2`, `Ey=1j*(Pp-Pm)/2`. | `notebook_sections/06_faraday.py` lines 106-125. | Used by `faraday_maps`. |
| Camera image | Binned and noisy image in normalized intensity units. | Bin by 15 grid cells, draw Poisson counts plus Gaussian read noise, divide by detected photons per pixel. | `notebook_sections/05_dgi.py` lines 19-38 and `notebook_sections/07_camera.py` lines 69-76. | Used in PCI/DGI/Faraday displayed camera frames. |

## Scalar phase pipeline

```text
Thomas-Fermi density
↓
column density n_col
↓
phi_peak(Delta_Hz, n_col[axis])
↓
phase map = phi_peak * _tf_profile(...)
↓
object field exp(i phase_map)
↓
scattered field exp(i phase_map) - 1
↓
finite-NA Fourier propagation
↓
mode-specific image formation
```

## Faraday rotation pipeline

```text
Thomas-Fermi density
↓
column density n_col
↓
phi_peak(Delta_Hz, n_col[axis])
↓
theta_F_peak = kappa_F * phi_peak
↓
rotation map = theta_F_peak * _tf_profile(...)
↓
circular fields exp(+i theta_F map), exp(-i theta_F map)
↓
finite-NA propagation of each circular component
↓
linear fields Ex, Ey
↓
Faraday analyzer model
```

## PCI image pipeline

```text
column density
↓
phi_peak
↓
phase map
↓
scattered field exp(i phase_map) - 1
↓
finite-NA pupil propagation
↓
PCI carrier t_p exp(i theta) + propagated scattered field
↓
image intensity |E|^2
↓
camera binning and noise
```

## DGI image pipeline

```text
column density
↓
phi_peak
↓
phase map
↓
scattered field exp(i phase_map) - 1
↓
finite-NA pupil propagation
↓
DGI residual carrier 10^(-OD/2) + propagated scattered field
↓
image intensity |E|^2
↓
camera binning and noise
```

## Dual-port Faraday image pipeline

```text
column density
↓
theta_F_peak
↓
rotation map
↓
propagate circular fields with phases ±theta_F
↓
recombine to Ex, Ey
↓
I_u = |Ex + Ey|^2 / 2
I_v = |Ex - Ey|^2 / 2
↓
S = (I_v - I_u) / (I_v + I_u)
↓
camera binning and noise where simulated
```

## Self review and remaining ambiguity

- `phi_peak` is documented and implemented as a scalar phase, but the Faraday section also discusses circular-component phase differences. Future derivations must not silently reinterpret `phi_peak` as a circular-component phase without documenting the convention change.
- `theta_F_peak` is implemented phenomenologically as `kappa_F * phi_peak`; the notebook does not derive `kappa_F` from Er vector polarizability.
- The camera image definition is split across multiple sections and duplicated for state-dependent images.
