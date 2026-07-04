# 01 Parameters
#
# Exported from: 1 calculations revised 2  multishot  6  extended.ipynb
# This file is a mechanical notebook-section export for refactoring.
# Keep physics, numerical algorithms, parameter values, and execution order equivalent to the notebook.


# %% [cell 3: markdown]
# ## 2. Apparatus and condensate parameters
#
# All inputs are gathered here. New entries relative to the original are flagged
# **[rev]**: the detection model and the loss-model controls.

# %% [cell 4: code]
# ---- 166Er and the 401 nm transition (K24 3.2.2) ----
m      = 166 * amu                # atomic mass
lam    = 401.0e-9                 # transition wavelength, m
k      = 2*np.pi/lam              # wavevector
Gamma  = 2*np.pi * 29.5e6         # natural linewidth, rad/s

# ---- Existing 401 nm imaging arm (K24 3.1.3) ----
P_probe_mW = 10.0                 # available probe power at cloud, mW
D_probe    = 24.0e-3              # collimated beam 1/e^2 diameter, m
f1, f2     = 150e-3, 300e-3       # 4f telescope focal lengths, m
pix_cam    = 2.93e-6              # effective camera pixel after 2x mag, m
tau        = 40e-6                 # reference imaging pulse duration, s (sections 6-15 baseline; section 15 folds in tau explicitly)

# ---- BEC operating point (K24 6.3.1) ----
N0     = 2.5e4                    # condensate atom number
a_s    = 72 * a0                  # s-wave scattering length
trap_Hz = np.array([293.0, 14.0, 233.0])   # trap frequencies (x, y, z), Hz
T_cloud = 200e-9                  # thermal cloud temperature, K

# ---- Phase-plate parameters for PCI ----
t_p   = 0.95                      # supplier-spec amplitude transmittance
theta = np.pi/2                   # plate retardation

# ---- [rev] Detection model (verify vs DCC3260M datasheet at lab visit) ----
QE_cam   = 0.40                   # quantum efficiency at 401 nm (CMOS ~0.3-0.5)
read_e   = 7.0                    # read noise, e- rms / pixel
QE_ideal = 1.0                    # reproduces the original notebook

# ---- [rev] Scattering / loss-model controls ----
use_peak_intensity = True         # on-axis Gaussian intensity (atoms at beam centre)
eta_coll = 1.3                    # effective atoms lost per scattered photon (collisional secondaries)

# %% [cell 5: markdown]
# ## 3. Atomic and optical constants
# $\sigma_0=3\lambda^2/(2\pi)$;  $I_\mathrm{sat}=\pi h c\Gamma/(3\lambda^3)$;
# $E_\mathrm{rec}=(\hbar k)^2/2m$, $T_\mathrm{rec}=E_\mathrm{rec}/k_B$.

# %% [cell 6: code]
sigma0 = 3*lam**2 / (2*np.pi)
Isat   = np.pi*h*c*Gamma / (3*lam**3)
E_rec  = (hbar*k)**2 / (2*m)
T_rec  = E_rec/kB
v_rec  = hbar*k/m
E_phot = h*c/lam
print(f"sigma_0  = {sigma0*1e4:.3e} cm^2")
print(f"I_sat    = {Isat*1e-1:.2f} mW/cm^2   (K24 quotes ~60)")
print(f"E_rec/kB = {T_rec*1e9:.1f} nK")
print(f"v_rec    = {v_rec*1e3:.3f} mm/s")
