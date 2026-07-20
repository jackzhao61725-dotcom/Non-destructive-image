# Simulation reference parameters

Last updated: 20 July 2026
Status: active dissertation screening contract

## Decision

The dissertation simulation will use one reproducible reference description of the condensate, optical system and detector. Parameters that have not yet been measured on the installed apparatus are assigned experimentally plausible values; they are not presented as measurements or calibrated apparatus specifications.

Chapter 5 will vary only the principal experimental design variables:

- absolute probe detuning, `|Delta|`;
- probe fluence, `F = P tau`;
- imaging readout mode.

The remaining condensate, optical and detector quantities are held fixed in the main screening calculation. They will be replaced by measured values after commissioning, after which the same scripts will regenerate the figures. Broad scans of every uncertain input are deliberately excluded because they would obscure the central signal-versus-disturbance comparison.

## Reference configuration

### Condensate and atomic inputs

| Quantity | Reference value | Status |
| --- | ---: | --- |
| Species | `166Er` | fixed |
| Transition wavelength | `401 nm` | literature/atomic input |
| Natural linewidth | `Gamma/2pi = 29.5 MHz` | literature input used by the maintained model |
| Initial condensate population | `2.5e4` | reference condensate input |
| Initial temperature | `200 nK` | reference condensate input |
| Trap frequencies `(x,y,z)` | `(293, 14, 233) Hz` | reference condensate input |
| Scattering length | `72 a0` | reference condensate input |

These quantities define the reference cloud used for the screening figures. They are not scanned in the main Chapter 5 design study.

### Optical and detector inputs

| Quantity | Reference value | Status |
| --- | ---: | --- |
| Numerical aperture | `NA = 0.080` | retained unchanged as a simulation reference; not an installed-arm measurement |
| Magnification | `M = 4` | provisional screening value; not an installed-arm measurement |
| DCC3260M physical pixel pitch | `5.86 um` | hardware value |
| Object-plane pixel pitch | `1.465 um` | `5.86 um / M` |
| Simulation camera binning | `15 x 15` high-resolution cells | represents `1.4648 um` in the object plane |
| Binned frame size | `68 x 68` pixels | after truncating the `1024 x 1024` grid to `1020 x 1020` |
| Quantum-efficiency factor | `QE = 0.60` | plausible screening value; not a DCC3260M measurement |
| Camera read noise | `sigma_r = 3 e- rms` per pixel and readout | plausible screening value; not a DCC3260M measurement |
| PCI phase-plate amplitude transmission | `t_p = 0.95` | retained Version 1 reference |
| PCI retardance | `eta = pi/2` | retained Version 1 reference |
| DGI stop optical depth | `OD_s = 4` | retained Version 1 reference |
| Faraday response coefficient | `kappa_F = 1` | illustrative structural setting only; not a calibrated `166Er` response |

The previous `QE = 0.40` value is no longer used because it is an unnecessarily low detector reference. The previously used `7 e- rms` is a DCC3260M catalogue upper bound rather than a nominal or measured read noise and is therefore not used as the main simulation value. The revised `QE = 0.60` and `sigma_r = 3 e- rms` values are intentionally simple reference assumptions. Their purpose is to provide realistic camera images and SNR estimates before detector characterisation, not to predict the installed camera exactly.

The official CS235MU response data give a typical quantum efficiency close to 0.6 at 401 nm and motivate the provisional QE scale. The DCC3260M datasheet specifies read noise only as less than `7 e- rms`; it does not establish that `7 e- rms` is the operating value of the installed camera.

Manufacturer sources:

- DCC3260M archived datasheet: <https://www.thorlabs.com/catalogpages/obsolete/2019/DCC3260M.pdf>
- CS235MU product page: <https://www.thorlabs.com/item/CS235MU>
- CS235MU raw QE data: <https://media.thorlabs.com/globalassets/family-pages/sharedassets/c/cs/cs235mu_quantum_efficiency.xlsx?v=1117122645>

### Reference operating point and sequence model

| Quantity | Reference value | Status |
| --- | ---: | --- |
| Absolute detuning | `|Delta|/2pi = 1.5 GHz` | marked reference point within the Chapter 5 scan |
| Probe power | `P = 1.0 mW` | reference division of the fluence |
| Exposure time | `tau = 90 us` | reference division of the fluence |
| Fluence | `F = 90 mW us` | reference operating fluence |
| Imaging direction | `x` | common comparison geometry |
| Recoil-energy convention | `2 E_rec` per absorption-spontaneous-emission cycle | retained Version 1 sequence convention |
| Reabsorption model | initial-density three-axis estimate, fixed within a sequence | retained Version 1 approximation |
| Condensate-depletion threshold | `30%` | analysis definition of an accepted frame |
| Stopping rule | strict integer; threshold-crossing pulse excluded | retained |

The split between `P` and `tau` is retained so that the reference exposure can be related to hardware. Within the far-detuned, low-saturation screening model, the principal probe-strength coordinate is the fluence `F = P tau`.

The current full heating-plus-reabsorption model gives ten strict accepted frames at this reference operating point. This frame count is unaffected by the QE and read-noise update because those detector parameters change the measured SNR, not the photon scattering or recoil heating. It remains conditional on the Version 1 disturbance model.

## Consequences for figures and numbers

The current maintained outputs use `M = 4`, `QE = 0.60` and
`sigma_r = 3 e- rms`. This contract controls:

- the detector-noise curves in Figure 3.2;
- all detected-electron count scales in Section 4.4;
- all Chapter 5 single-frame SNR maps and quoted SNR values;
- all stochastic camera-image examples;
- all accumulated SNR values that include camera noise.

The following quantities are unchanged by this detector update:

- the scalar phase and Faraday rotation maps;
- the finite-aperture fields and the noiseless Figure 4.2 images, because `NA = 0.080` is not changed in this update;
- the propagated field `u` and the normalised PCI, DGI and Faraday image-plane values reported in Section 4.3;
- photon scattering, recoil heating and the strict ten-frame sequence length;
- the ideal transfer-function and fluence-scaling arguments of Chapter 3;
- the statement that absolute Faraday values remain uncalibrated while `kappa_F` is unknown.

Outputs generated with `M = 2`, `QE = 0.5851647` or
`sigma_r = 7 e- rms` are superseded and must not be mixed with the active
dissertation figures.

In the dissertation, `NA = 0.080`, `M = 4`, `QE = 0.60` and
`sigma_r = 3 e- rms` are therefore simulation reference values, not measured
properties of the installed apparatus.

## Commissioning update route

The post-commissioning update should replace, at minimum:

1. the effective pupil or measured point-spread function;
2. the object-plane magnification and pixel sampling;
3. the effective photon-to-electron conversion, including optical transmission and camera response;
4. the measured read noise and gain;
5. the effective Faraday coefficient `kappa_F`;
6. the net condensate disturbance observed across repeated exposures.

After these inputs are frozen, the same configuration-driven scripts should regenerate Figures 3.2, 4.2 and the Chapter 5 figures. Agreement with data used to determine the inputs is calibration; assessment against separate operating conditions is required before calling the updated model experimentally validated.
