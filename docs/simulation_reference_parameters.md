# Simulation reference parameters

Last updated: 22 July 2026
Status: active optical/detector contract and implemented Version 1 sequence record

Update trigger: a change to the reference condensate, atomic response,
detector, sampling, aperture, operating point or disturbance model. Dependent
configs, figures and regression tests must be updated in the same numerical
change.

## Decision

The dissertation simulation will use one reproducible reference description of the condensate, optical system and detector. Parameters that have not yet been measured on the installed apparatus are assigned experimentally plausible values; they are not presented as measurements or calibrated apparatus specifications.

Chapter 5 will vary only the principal experimental design variables:

- absolute probe detuning, `|Delta|`;
- probe fluence, `F = P tau`;
- imaging readout mode.

The remaining condensate, optical and detector quantities are held fixed in the main screening calculation. They will be replaced by measured values after commissioning, after which the same scripts will regenerate the figures. Broad scans of every uncertain input are deliberately excluded because they would obscure the central signal-versus-disturbance comparison.

The active detector, sampling and Faraday-response contract is implemented in
`configs/dissertation_v3_orca_fusion.json` and the dependent Figure 4.2-5.4
configs. It uses signed `kappa_F=-45/91`, effective `NA=0.130` and the
manufacturer-typical ORCA-Fusion Ultra quiet value `sigma_r=0.7 e- rms`. The
sealed reconstruction outputs retain `NA=0.080`, `kappa_F=1` and
`sigma_r=1.4 e- rms` as a historical operator and are not current performance
predictions.

The `N0=2.5e4`, `T=200 nK` state below is the compact Oxford condensate-core
surrogate used by the optical screening calculation. It is not the initial
thermodynamic state of the approved multiframe replacement. That replacement
uses direct Oxford 300 ms `(N0,Nth,T)` triplets and is specified in
`multiframe_heating_model_optimisation.md`; it has not yet replaced the stored
Version 1 Figures 5.2 and 5.4.

## Reference configuration

### Condensate and atomic inputs

| Quantity | Reference value | Status |
| --- | ---: | --- |
| Species | `166Er` | fixed |
| Transition wavelength | `401 nm` | literature/atomic input |
| Natural linewidth | `Gamma/2pi = 29.5 MHz` | maintained numerical input; its modern primary source must be identified before dissertation use |
| Initial condensate population | `2.5e4` | reference condensate input |
| Initial temperature | `200 nK` | reference condensate input |
| Trap frequencies `(x,y,z)` | `(293, 14, 233) Hz` | reference condensate input |
| Scattering length | `72 a0` | reference condensate input |

These quantities define the reference cloud used for the screening figures. They are not scanned in the main Chapter 5 design study.

### Optical and detector inputs

| Quantity | Reference value | Status |
| --- | ---: | --- |
| Effective numerical aperture | `NA = 0.130` | active optical-design input; not an installed-arm measurement |
| Camera | Hamamatsu ORCA-Fusion C14440-20UP | selected detector |
| Magnification | `M = 10` | design input; the code does not infer it from the legacy relay focal lengths |
| Physical camera pixel pitch | `6.5 um` | manufacturer value |
| Object-plane pixel pitch | `0.650 um` | `6.5 um / M` |
| Simulation camera sampling | centred physical-pixel area integration | avoids treating the non-integer high-resolution-to-camera ratio as block binning |
| Analysis crop | `153 x 153` pixels per port | centred coverage of the numerical field; not a directly programmed hardware ROI |
| Quantum-efficiency factor | `QE = 0.65` | manufacturer-typical value at 400 nm, used at 401 nm |
| Camera readout | Ultra quiet scan, 16-bit output | active low-noise operating scenario |
| Camera read noise | `sigma_r = 0.7 e- rms` per physical pixel and readout | manufacturer-typical Ultra quiet value; requires installed-camera dark-frame verification |
| Camera integration | reference `360.6 us`; minimum input `280 us`, `80 us` timing increment | `360.6 us` is the shortest calculated Ultra quiet setting containing the `300 us` reference optical pulse |
| Full-resolution readout | `184.4 ms` (`5.42 fps`) | Ultra quiet catalogue timing; sub-array rate depends on vertical extent |
| Hardware sub-array step | `4` pixels/lines | Normal Area mode through DCAM-API |
| First-frame count scale | `735.2695759 e- / I0 pixel` at `F = 300 mW us` | derived from the optical/detector model, not a catalogue camera specification |
| PCI phase-plate amplitude transmission | `t_p = 0.95` | retained Version 1 reference |
| PCI retardance | `eta = pi/2` | retained Version 1 reference |
| DGI stop optical depth | `OD_s = 4` | retained Version 1 reference |
| Atomic Faraday conversion | `kappa_F = -45/91 = -0.4945...` | signed ideal `166Er` estimate for the fully spin-polarised axial case and isolated 401-nm transition |
| Effective apparatus response | `TBD` | requires paired RAI/Faraday calibration; includes measured optical and multilevel departures without redefining `kappa_F` |

The detector values are declared screening inputs. The `0.7 e- rms` value is a
manufacturer-typical Ultra quiet scenario, not a guaranteed or measured value
for the installed camera. The `0.650 um` object-plane pixel is integrated
directly over the numerical field and is not approximated by a rounded number
of simulation cells. The `153 x 153` product is an analysis crop. Commissioning
must acquire a legal four-pixel-step ROI containing both analyser ports and
must verify either a global-reset trigger or a common-exposure interval in
which the short optical pulse reaches every relevant row.

The preceding DCC3260M values (`M=4`, `QE=0.60`, `sigma_r=3 e- rms`) remain in
the historical v2 config. The sealed ORCA-Fusion reconstruction studies use the
separate historical contract `NA=0.080`, `kappa_F=1` and
`sigma_r=1.4 e- rms`. Neither historical contract is mixed with the active
screening figures.

### Faraday-response convention and frozen-output boundary

The dissertation defines the total phases of the two circular components as
`Phi_+` and `Phi_-`, and the polarisation rotation as

```text
theta_F = (Phi_+ - Phi_-) / 2,
kappa_F = theta_F / phi.
```

For the fully spin-polarised `166Er` reference state probed along its
quantisation axis, the isolated-transition calculation gives
`kappa_F = -45/91`, so the reference scalar phase `phi_0 = 0.203 rad`
corresponds to `theta_F,0 = -0.1004 rad`; its magnitude is `5.75 degrees`.
This is an atomic-response estimate, not a measurement of the installed
apparatus. It assumes a weak
probe, the selected isolated 401-nm transition, the stated helicity and detuning
convention, negligible optical pumping within one exposure and negligible
circular dichroism. Additional excited levels, a non-axial geometry or an
imperfectly prepared Zeeman state require a multilevel calculation or an
effective experimental calibration.

The historical notebook and sealed reconstruction operator instead set
`kappa_F=1` before Jones propagation and
passed the resulting `theta_F` directly to circular phases `+theta_F` and
`-theta_F`; no later stage supplied the missing factor of one half. Stored
reconstruction studies generated under that convention remain frozen
method-development evidence. Their absolute Faraday amplitudes, SNRs and
usable-frame predictions are not current `166Er` performance estimates. The
active Figure 4.2-5.4 family has been regenerated under the signed coefficient.

The maintained linewidth input also needs bibliographic closure. McClelland and
Hanssen (2006) supports the open-transition and optical-pumping physics, not the
linewidth value. The early direct measurement by J. J. McClelland (2006) gives
`Gamma/2pi = 35.6 +/- 1.2 MHz`; continued use of `29.5 MHz` therefore requires
the exact modern primary source and an explicit reason for selecting it.

Manufacturer sources:

- ORCA-Fusion product page: <https://camera.hamamatsu.com/us/en/product/camera/C14440-20UP.html>
- ORCA-Fusion technical datasheet: <https://camera.hamamatsu.com/content/dam/hamamatsu-photonics/sites/documents/99_SALES_LIBRARY/sys/SCAS0136E_C14440-20UP.pdf>
- ORCA-Fusion instruction manual: <https://www.hamamatsu.com/content/dam/hamamatsu-photonics/sites/static/sys/en/manual/C14440-20UP_IM_En.pdf>

### Reference operating point and implemented Version 1 sequence

| Quantity | Reference value | Status |
| --- | ---: | --- |
| Absolute detuning | `|Delta|/2pi = 1.5 GHz` | marked reference point within the Chapter 5 scan |
| Fluence scan | `90-300 mW us` | active Chapter 5 interval |
| Representative fluences | `90, 150, 300 mW us` | A, B and C conditions in Figures 5.1 and 5.4 |
| Probe power | `P = 1.0 mW` | reference division of the fluence |
| Optical probe duration | `tau = 300 us` | one division of the first-frame reference; not the general camera integration time |
| First-frame reference | `F = 300 mW us` | high-SNR reference for Figure 5.1 and the screen |
| Imaging direction | `x` | common comparison geometry |
| Recoil-energy convention | `2 E_rec` per absorption-spontaneous-emission cycle | frozen Version 1 sequence convention |
| Reabsorption model | initial-density three-axis estimate, fixed within a sequence | frozen Version 1 approximation |
| Condensate-depletion threshold | `30%` | analysis definition of an accepted frame |
| Stopping rule | strict integer; threshold-crossing pulse excluded | retained |

This sequence implementation infers a total atom number from the compact
surrogate and applies an ideal saturated-gas population relation after each
pulse. Its generated `N_dep` and `N_use` values are frozen Version 1 screening
outputs pending the approved non-saturation replacement. They remain available
for regression and provenance; they are not current quantitative predictions of
an Oxford bimodal sequence.

The split between `P` and `tau` is retained so that the optical pulse can be
related to hardware. Within the far-detuned, low-saturation screening model,
the principal probe-strength coordinate is `F = P tau`. The optical pulse may
be shorter than the camera integration only when it lies inside a verified
common-exposure interval for both analyser ports. The `300 us` reference pulse
therefore uses an approximately `360.6 us` camera integration under the Ultra quiet exposure
quantisation.

At `F = 300 mW us`, the active first-frame analytic central-`5x5` SNR values
are `7.8000` for dark-field Faraday and `15.3446` for dual-port Faraday. At
`F = 90, 150, 300 mW us`, the frozen Version 1 strict `(N_dep, N_use)` pairs are
`(10, 3)`, `(6, 6)` and `(3, 3)`. Detector parameters change SNR and therefore
`N_use`; they do not change photon scattering or the depletion-only `N_dep` for
a fixed optical pulse within that implementation. All sequence counts remain
conditional on the Version 1 disturbance model until replacement.

## Consequences for figures and numbers

The detector-dependent Figure 4.2-5.4 family uses `M=10`, `QE=0.65`, signed
`kappa_F=-45/91`, effective `NA=0.130` and the Ultra quiet
`sigma_r=0.7 e- rms` scenario. This contract controls:

- the regenerated Faraday panels in Figure 4.2;
- the Figure 5.1 scan over `F=90-300 mW us`, with A, B and C at
  `90, 150, 300 mW us` and the first-frame reference at `F=300 mW us`;
- the frozen Version 1 Figure 5.2 usable-frame screen and fixed-detuning
  operating band pending heating replacement;
- the frozen Version 1 Figure 5.4 noisy sequences at the same representative
  fluences pending heating replacement;
- all detected-electron count scales, stochastic frames and analytic SNR values
  reported for this result family, subject to the sequence-status distinction
  above.

Changing `NA` changes the propagated scalar and Faraday fields. Changing the
signed atomic coefficient changes the Faraday field amplitudes and the DPFI
sign. Camera QE and read noise change detected counts and SNR. Within the frozen
Version 1 implementation, none of these changes alters photon scattering,
recoil heating or `N_dep` at fixed detuning and optical fluence. Absolute
Faraday density inference still requires an effective apparatus-response
calibration.

Outputs generated with the DCC3260M/M4 detector contract, the earlier
`M=2`, `QE=0.5851647`, `sigma_r=7 e- rms` trial contract, or the sealed
ORCA-Fusion reconstruction contract `NA=0.080`, `kappa_F=1`,
`sigma_r=1.4 e- rms` are historical. They must not be mixed with the current
signed-response Figure 4.2-5.4 family.

In the dissertation, `NA=0.130`, `M=10`, `QE=0.65`, Ultra quiet 16-bit output
and `sigma_r=0.7 e- rms` are explicit simulation inputs. Only the atomic
coefficient is analytically fixed; the effective NA, QE, read noise, timing and
apparatus response remain commissioning quantities.

## Replaceable experimental inputs

The configuration keeps the following quantities replaceable without changing
the reconstruction or figure-generation code:

1. the effective pupil or measured point-spread function;
2. the object-plane magnification and pixel sampling;
3. the effective photon-to-electron conversion, including optical transmission and camera response;
4. the measured read noise and gain;
5. multilevel corrections to the atomic Faraday response and the independently
   calibrated effective apparatus response;
6. the net condensate disturbance observed across repeated exposures.

The same configuration-driven scripts regenerate the dependent figures when
any of these inputs changes. Figure 4.2 is unaffected by detector-noise changes
because it contains noiseless normalised image-plane observables, but it must be
regenerated when `NA` or the atomic Faraday coefficient changes.
