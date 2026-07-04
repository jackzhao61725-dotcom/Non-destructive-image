# 00 Imports
#
# Exported from: 1 calculations revised 2  multishot  6  extended.ipynb
# This file is a mechanical notebook-section export for refactoring.
# Keep physics, numerical algorithms, parameter values, and execution order equivalent to the notebook.


# %% [cell 0: markdown]
# # Phase-contrast and dark-ground imaging on the Oxford Er apparatus — design calculations (revised)
#
# **Revision note.** This notebook supersedes the original design calculations. The
# derivations and the reference operating point of the original are preserved and
# reproduced exactly; the revisions are:
#
# 1. **Scattering uses the on-axis (peak Gaussian) intensity**, not the area
#    average — atoms sit at the beam centre, where a Gaussian probe is $2\times$
#    the area-averaged intensity. This doubles photons scattered per atom and
#    halves every destruction-budget shot count. The original convention is kept
#    behind a flag (`use_peak_intensity`).
# 2. **The destruction limit is an atom-loss model, not a heating model.** For
#    $^{166}$Er at 401 nm, $E_\mathrm{rec}/k_B=359$ nK is $7.5\times$ the chemical
#    potential $\mu/k_B=48$ nK: one scattered photon ejects an atom rather than
#    warming the cloud. A family of physically motivated thresholds is compared (§9).
# 3. **A realistic photon-detection model** — finite quantum efficiency, camera
#    read noise, and contrast dilution by the finite-NA PSF — replaces the
#    idealised $\mathrm{QE}=1$, noiseless, infinite-resolution SNR (§8).
# 4. **A new image-formation section** (§7) propagates the actual Thomas-Fermi
#    cloud through the 4f system, the Fourier-plane optic, and the camera by
#    explicit Fourier optics.
# 5. **A parameter-sensitivity analysis** (§8.5) quantifies how the headline SNR and
#    shot-budget numbers move over the plausible range of the uncertain inputs (camera
#    QE and read noise, collisional multiplier $\eta$), with a numerical demonstration
#    of the accumulated-SNR invariance and motional / depth-of-focus sanity checks.
# 6. **Normalisation conventions are made explicit** (§5.3, §7.2): the linear PCI signal
#    is $t_p^2+2t_p\varphi$ relative to the incident intensity, or $1+(2/t_p)\varphi$
#    relative to the plate-in reference; the transfer-curve figure now plots the genuine tangent.
#
# **Destruction model — refined in §12–§18.** Item 2 and §9 treat destruction as clean atom *loss*, argued from $E_\mathrm{rec}\gg\mu$. Sections 15 onward correct the criterion: what decides whether a recoiled atom leaves is the recoil energy versus the **ODT trap depth** ($\sim\mu$K), not versus $\mu$. Since $E_\mathrm{rec}=359$ nK is well *below* the trap depth, the atom stays trapped and **heats** the cloud (melting a condensate that is only $\sim$22 % condensed) instead of escaping. The realistic model from §12 on is therefore *heating*, with clean loss as an optimistic upper bound; the two bracket the truth, pinned by one measured loss-vs-shots point.
#
# Apparatus baseline: Kucera (DPhil thesis, 2024), **K24** — §3.1.3 imaging arm;
# §3.2.2 401 nm laser and transition; §6.3.1 BEC operating point.
#
# ### Standard references
# * W. Ketterle, D. S. Durfee, D. M. Stamper-Kurn, Varenna lecture notes (1999), arXiv:cond-mat/9904034 — dispersive imaging, §3.
# * R. Meppelink et al., Phys. Rev. A **81**, 053632 (2010) — quantitative PCI methodology (primary).
# * M. R. Andrews et al., Science **273**, 84 (1996) — first non-destructive BEC observation (DGI).
# * M. Gajdacz et al., Rev. Sci. Instrum. **84**, 083105 (2013) — minimally-destructive imaging; destructiveness as fractional loss per image.
#
#
# **Second extension (this revision).** §17 adds Faraday (polarization-rotation) imaging as a
# third dispersive modality alongside PCI/DGI, comparing the dark-field (Gajdacz et al. 2013) and
# dual-port (Kaminski et al. 2012) detection schemes. §18 walks the complete simulation pipeline
# stage by stage -- from the Thomas-Fermi condensate through to the final camera frame -- for all
# four imaging modes side by side, so every headline number in this notebook is traceable back to
# where it came from. A full reference list closes the notebook.
#
# **Third extension (this revision).** Part IV (§19) replaces the earlier validation-style anatomy
# with a guided, step-by-step walk of a single probe shot through the machine: every figure is the
# *actual simulated field or image at that exact line of the algorithm*, in the order the code
# computes it, and each step carries its theory, its physical meaning, and the code that executes
# it. It ends with the multi-shot run rendered as a filmstrip of the dying condensate.

# %% [cell 1: markdown]
# ## 1. Imports and fundamental constants

# %% [cell 2: code]
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# %matplotlib inline
plt.rcParams.update({
    'figure.figsize': (8, 5), 'font.size': 11, 'axes.labelsize': 12,
    'axes.titlesize': 12, 'legend.fontsize': 10, 'figure.dpi': 110,
})
rng = np.random.default_rng(7)   # reproducible noise for the image simulation

# Fundamental constants (SI)
hbar = 1.054571817e-34       # J s
h    = 2*np.pi*hbar
c    = 2.99792458e8          # m/s
kB   = 1.380649e-23          # J/K
amu  = 1.66053907e-27        # kg
a0   = 5.29177211e-11        # Bohr radius, m
