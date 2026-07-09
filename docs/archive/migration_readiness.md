# Migration Readiness Review

This review answers whether each current helper area has enough baseline information to replace notebook-local code safely.

| Area | Status | Explanation | Missing baseline |
|---|---|---|---|
| Atomic model | BLOCKED | Helper formulas exist, but migration should wait for fresh notebook-generated scalar and array baselines. | `atomic/tf_state.npz`, `atomic/column_density.npz`, `atomic/recoil.npz` |
| Thomas-Fermi profile | BLOCKED | Helper exists; needs profile arrays from the original notebook grid. | `atomic/tf_profile.npz` |
| Fourier propagation | BLOCKED | Helper exists; needs reference propagated fields/images from original notebook execution. | `imaging/pci_image.npz`, `imaging/dgi_image.npz` |
| Camera | BLOCKED | Helper exists; needs deterministic pre-noise and fixed-seed noisy camera arrays. | `imaging/ideal_camera_image.npz`, `imaging/shot_noise_image.npz` |
| Light-atom interaction | BLOCKED | Helper exists; needs reference phase, OD, scattering, and recoil quantities saved from notebook execution. | `atomic/phase_map.npz`, `atomic/scattering.npz`, `atomic/recoil.npz` |
| PCI | BLOCKED | No PCI-specific helper migration should occur before image baselines exist. | `imaging/pci_image.npz` |
| DGI | BLOCKED | No DGI-specific helper migration should occur before image baselines exist. | `imaging/dgi_image.npz` |
| Faraday | BLOCKED | No Faraday helper migration exists yet and no Faraday array baseline exists. | `imaging/faraday_image.npz` |
| Multi-shot simulation | BLOCKED | Multi-shot code remains notebook-local and needs sequence baselines. | `multishot/evolution.npz`, selected shot files |

## Summary

No notebook section is ready for complete migration yet. The helper package is useful preparation, but scientific baseline arrays must be generated first.
