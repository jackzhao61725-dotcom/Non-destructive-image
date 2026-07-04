# 03 Light-Atom Interaction
#
# Exported from: 1 calculations revised 2  multishot  6  extended.ipynb
# This file is a mechanical notebook-section export for refactoring.
# Keep physics, numerical algorithms, parameter values, and execution order equivalent to the notebook.


# %% [cell 9: markdown]
# ## 5. Dispersive phase shift
#
# ### 5.1 Two-level polarisability
# With dimensionless detuning $\delta\equiv2\Delta/\Gamma$,
# $\mathrm{Re}\,\alpha\propto-\delta/(1+\delta^2)$, $\mathrm{Im}\,\alpha\propto1/(1+\delta^2)$.
# For $|\delta|\gg1$ absorption ($\propto1/\delta^2$) is suppressed faster than refraction
# ($\propto1/\delta$): a non-resonant probe acquires phase without much attenuation.
#
# ### 5.2 Column phase
# $\varphi=\dfrac{\sigma_0\tilde n}{2}\dfrac{\delta}{1+\delta^2}\;\xrightarrow{|\delta|\gg1}\;\dfrac{\sigma_0\tilde n}{2\delta}.$
#
# ### 5.3 Linear vs Meppelink regime
# $I/I_0\approx t_p^2+2t_p\varphi$ (incident-$I_0$ convention; see §7.2) needs $\varphi\lesssim0.5$;
# otherwise use the full periodic $I/I_0=|t_p e^{i\theta}+e^{i\varphi}-1|^2$ (Meppelink et al. 2010),
# fitting $(t_p,\theta)$ from a detuning scan, and beyond $\varphi=\pi$ the inversion is multivalued.

# %% [cell 10: code]
from non_destructive_image import (
    dimensionless_detuning,
    faraday_rotation_angle,
    intensity_at_atoms as _intensity_at_atoms_helper,
    reabsorption_fraction,
    residual_optical_depth,
    scalar_phase_shift,
    scattered_photons_per_atom,
)

def delta_of(Delta_Hz):
    return dimensionless_detuning(Delta_Hz, Gamma)

def phi_peak(Delta_Hz, n_col_peak):
    return scalar_phase_shift(Delta_Hz, n_col_peak, sigma0, Gamma)

def od_resonant_equiv(Delta_Hz, n_col_peak):
    # residual on-resonance-scaled optical depth at this detuning: OD = sigma0 n / (1+delta^2)
    return residual_optical_depth(Delta_Hz, n_col_peak, sigma0, Gamma)

print(f"At Delta=1.5 GHz: delta = {delta_of(1.5e9):.1f}")
print(f"  phi_x (across cigar) = {phi_peak(1.5e9, n_col[0]):.3f} rad")
print(f"  phi_y (along cigar)  = {phi_peak(1.5e9, n_col[1]):.3f} rad")
print(f"  residual OD_x        = {od_resonant_equiv(1.5e9, n_col[0]):.4f}  (absorption negligible)")

# %% [cell 11: markdown]
# ### Table 1 — Peak phase shift at representative detunings

# %% [cell 12: code]
detunings_GHz = [0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0]
fig, ax = plt.subplots(figsize=(9.5, 3.3)); ax.axis('off')
col_labels = ['D (GHz)', 'delta', 'phi_x\n(across, 1.19um)',
              'phi_y\n(along, 24.8um)', 'phi_z\n(across, 1.49um)', 'Regime']
rows = []
for Dg in detunings_GHz:
    d = delta_of(Dg*1e9)
    px, py, pz = (phi_peak(Dg*1e9, n_col[i]) for i in range(3))
    pmax = max(px, pz)
    reg = 'linear PCI' if pmax < 0.5 else ('Meppelink' if pmax < np.pi else 'wrapped')
    rows.append([f'{Dg:.2f}', f'{d:.0f}', f'{px:.3f}', f'{py:.2f}', f'{pz:.3f}', reg])
tb = ax.table(cellText=rows, colLabels=col_labels, cellLoc='center', loc='center',
              colWidths=[.10,.12,.20,.20,.20,.18])
tb.auto_set_font_size(False); tb.set_fontsize(9.5)
for (i,j),cl in tb.get_celld().items():
    cl.set_height(0.16 if i==0 else 0.11)
    if i==0: cl.set_facecolor('#1a3a5c'); cl.set_text_props(color='white', weight='bold')
ax.set_title('Table 1.  Peak phase shift through the in-trap BEC (N0=2.5e4, a_s=72a0)',
             fontsize=10.5, pad=10, loc='left')
plt.tight_layout(); plt.show()

# %% [cell 13: markdown]
# ## 6. Photon scattering and per-shot loss
#
# ### 6.1 Scattering rate and on-axis intensity
# $R_\mathrm{scatt}=\dfrac{\Gamma}{2}\dfrac{s}{1+s+\delta^2}$, $s=I/I_\mathrm{sat}$.
# **[rev]** The relevant intensity is the value *at the atoms*. For a Gaussian probe of
# $1/e^2$ diameter $D$ carrying power $P$, the on-axis intensity is
# $I_\mathrm{peak}=\dfrac{2P}{\pi (D/2)^2}$ — twice the area-averaged
# $\bar I = P/[\pi(D/2)^2]$ used in the original. The cloud ($R\le25\,\mu$m) sits well
# within the beam centre, so $I_\mathrm{peak}$ is correct for scattering; the flag
# `use_peak_intensity` recovers the original $\bar I$ for traceability.
#
# ### 6.2 Photons per atom per pulse
# $N_\gamma=R_\mathrm{scatt}\,\tau$. In the far-detuned limit $N_\gamma\approx(\Gamma/2)(s/\delta^2)\tau$.
#
# ### 6.3 Why heating is the wrong picture, and the loss model
# $E_\mathrm{rec}/k_B=359$ nK while $\mu/k_B=48$ nK, so a single recoil ejects an atom
# from the condensate. Rather than accumulating temperature, the condensate number
# decays as $N_0(s)=N_0\exp(-\eta\,N_\gamma\,s)$ over $s$ shots, where $\eta\ge1$ counts
# **collisional secondaries**: a recoiling atom traversing the cloud can knock out
# further atoms. The mean free path $\ell=1/(n_\mathrm{peak}\sigma_\mathrm{el})$,
# $\sigma_\mathrm{el}=8\pi a_s^2$, sets whether secondaries matter along each axis.
#
# *Refined in §12:* this clean-loss picture is the **optimistic bound**. Because $E_\mathrm{rec}=359$ nK is below the ODT trap depth ($\sim\mu$K), a recoiled atom usually stays trapped and *heats* the cloud rather than escaping; §12 adopts that heating model as the realistic one, with clean loss as the upper bound.

# %% [cell 14: code]
def intensity_at_atoms(P_mW):
    return _intensity_at_atoms_helper(P_mW, D_probe, use_peak_intensity)

def N_scatt(Delta_Hz, P_mW, tau_s=None):
    if tau_s is None: tau_s = tau
    return scattered_photons_per_atom(
        Delta_Hz,
        P_mW,
        tau_s,
        Isat,
        Gamma,
        D_probe,
        use_peak_intensity,
    )

def Nmax_loss(Delta_Hz, P_mW, frac, tau_s=None, eta=None):
    # shots until a fraction `frac` of the condensate is lost:  N0(s)=N0 exp(-eta N_gamma s)
    if eta is None: eta = eta_coll
    return -np.log(1-frac) / (eta*N_scatt(Delta_Hz, P_mW, tau_s))

def Nmax_heat(Delta_Hz, P_mW, eps_K, tau_s=None):
    # original-style cumulative-heating threshold (eta=1): shots until accumulated energy = eps
    return (eps_K/T_rec) / N_scatt(Delta_Hz, P_mW, tau_s)

# collisional mean free path vs cloud size
sigma_el = 8*np.pi*a_s**2
mfp = 1.0/(n_peak*sigma_el)
print(f"elastic cross-section sigma_el = {sigma_el*1e4:.2e} cm^2")
print(f"mean free path at n_peak       = {mfp*1e6:.1f} um")
print(f"  vs  R = ({R[0]*1e6:.2f}, {R[1]*1e6:.2f}, {R[2]*1e6:.2f}) um")
print(f"  -> across-cigar recoils ({R[0]*1e6:.1f}, {R[2]*1e6:.1f} um) mostly escape;")
print(f"     along-cigar recoils ({R[1]*1e6:.1f} um) collide -> eta ~ 1.3-1.5")
print(f"\nN_gamma at (1.5 GHz, 2 mW), peak intensity = {N_scatt(1.5e9, 2.0):.2e}")
print(f"  (original area-averaged convention would give {N_scatt(1.5e9,2.0)/2:.2e})")

# --- motional sanity check: is the cloud frozen during the probe pulse? ---
disp_recoil = v_rec * tau                     # displacement of a single-recoil atom during tau
T_trap_min = 1/trap_Hz.max()                  # shortest trap period
print(f"\n[check] recoil displacement during {tau*1e6:.0f} us pulse = {disp_recoil*1e9:.0f} nm "
      f"(<< pixel {pix_cam/(f2/f1)*1e6:.2f} um, resolution {0.61*lam/((D_probe/2)/f1)*1e6:.2f} um)")
print(f"[check] shortest trap period {T_trap_min*1e3:.2f} ms >> pulse {tau*1e6:.0f} us "
      f"-> in-trap dynamics frozen during the shot; no motional blur")

# %% [cell 32: code]
# ===== CELL D: pulse duration & realistic destruction models =====
from scipy.special import zeta as _zeta
zeta3, zeta4 = float(_zeta(3)), float(_zeta(4))

def reabs_frac(Delta_Hz):
    """Angle-averaged reabsorption probability of a spontaneously emitted (Rayleigh)
    photon: residual OD along each principal axis, averaged over emission solid angle.
    Grows as ~1/delta^2, so it bites hardest at low detuning / high OD."""
    return reabsorption_fraction(Delta_Hz, n_col, sigma0, Gamma)

def Nmax_cleanloss(Delta_Hz, P_mW, frac, tau_s=None, eta=None):
    """Optimistic bound: every recoiled atom is promptly LOST (eta secondaries each)."""
    if eta is None: eta = eta_coll
    return -np.log(1-frac)/(eta*N_scatt(Delta_Hz, P_mW, tau_s))

def Nmax_heating(Delta_Hz, P_mW, frac, tau_s=None, reabs=None):
    """Realistic for a deep trap (E_rec = 359 nK << ODT depth ~uK): the recoiled atom
    does NOT leave - it stays, thermalises, and HEATS the cloud, melting the (~22%)
    condensate. Closed form from the trapped-Bose-gas energy E(T)=A T^4 below Tc,
    A = 3(zeta4/zeta3) kB/Tc^3; condensate number N0(s) = Ntot[1-(T/Tc)^3]."""
    if reabs is None: reabs = reabs_frac(Delta_Hz)
    T0, Tc = T_cloud, Tc_sc
    fc0 = 1-(T0/Tc)**3
    Ts  = Tc*(1-(1-frac)*fc0)**(1/3)
    dE  = N_scatt(Delta_Hz, P_mW, tau_s)*(1+reabs)*E_rec
    A   = 3*(zeta4/zeta3)*kB/Tc**3
    return A*(Ts**4 - T0**4)/dE

# ---- pulse-duration caps ----
T_trap_min = 1.0/trap_Hz.max()
NA_e   = (D_probe/2)/f1; res_e = 0.61*lam/NA_e
res_hi = 0.61*lam/0.40
tau_trap    = 0.1*T_trap_min
tau_blur_e  = res_e/v_rec
tau_blur_hi = res_hi/v_rec
tau_max_now = min(tau_trap, tau_blur_e)
tau_max_hi  = min(tau_trap, tau_blur_hi)

print("Realistic destruction at the reference point (1.5 GHz, 2 mW, 15 us), 30% condensate loss:")
print(f"   clean-loss (eta=1.3, optimistic)   N_max = {Nmax_cleanloss(1.5e9,2.0,0.30):5.0f}")
print(f"   heating-melt (recoil stays)        N_max = {Nmax_heating(1.5e9,2.0,0.30,reabs=0.0):5.0f}")
print(f"   heating-melt + reabsorption        N_max = {Nmax_heating(1.5e9,2.0,0.30):5.0f}   <- realistic")
print(f"\nreabsorption fraction: {reabs_frac(1.0e9):.3f} @1GHz, {reabs_frac(1.5e9):.3f} @1.5GHz, "
      f"{reabs_frac(2.0e9):.3f} @2GHz  (grows ~1/delta^2)")
print("\nPulse-duration caps:")
print(f"   shortest trap period {T_trap_min*1e3:.2f} ms -> tau < 0.1 T = {tau_trap*1e6:.0f} us (no in-trap motional blur)")
print(f"   recoil blur < resolution: tau < {tau_blur_e*1e6:.0f} us (NA={NA_e:.3f}, dx={res_e*1e6:.2f} um) "
      f"| {tau_blur_hi*1e6:.0f} us (NA=0.40, dx={res_hi*1e6:.2f} um)")
print(f"   => tau_max = {tau_max_now*1e6:.0f} us (existing arm),  {tau_max_hi*1e6:.0f} us (future objective)")
print("\nKEY INVARIANT  N_max x tau is fixed by (Delta, P)  [destruction budget = total integration time]:")
print(f"   at (1.5 GHz, 2 mW): N_max x tau = {Nmax_heating(1.5e9,2.0,0.30)*15:.0f} us, however it is sliced into shots.")

# %% [cell 49: code]
# ============================================================================
# 17.1 FARADAY ROTATION ANGLE -- same dispersive lineshape as the scalar phase (Sec. 5),
#      scaled by the vector-coupling factor kappa_F
# ============================================================================
#  >>>>>  kappa_F is a PLACEHOLDER pending Er atomic-structure input -- see markdown above  <<<<<
kappa_F = 1.0    # [TBD/rev] vector-polarizability fraction for the stretched-state 401 nm
                 # cycling transition; kappa_F=1 is the idealised maximal-coupling bound.
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def theta_F_peak(Delta_Hz, n_col_peak):
    """Peak Faraday-rotation angle (rad): identical dispersive lineshape to phi_peak (Sec. 5),
    scaled by kappa_F. We isolate the pure-rotation observable (phi_+ = +theta_F, phi_- = -theta_F),
    i.e. a fully spin-polarized column with the small common-mode scalar phase neglected."""
    return faraday_rotation_angle(Delta_Hz, n_col_peak, sigma0, Gamma, kappa_F)

# ---- Table 2: rotation angle vs phase shift at the same detunings as Table 1 ----
fig, ax = plt.subplots(figsize=(9.2, 3.0)); ax.axis('off')
col_labels = ['D (GHz)', 'delta', 'phi_x (PCI)', 'theta_F,x (Faraday)', 'ratio theta_F/phi']
rows = []
for Dg in detunings_GHz:
    d = delta_of(Dg*1e9)
    px = phi_peak(Dg*1e9, n_col[0]); tx = theta_F_peak(Dg*1e9, n_col[0])
    rows.append([f'{Dg:.2f}', f'{d:.0f}', f'{px:.3f}', f'{tx:.3f}', f'{tx/px:.2f}'])
tb = ax.table(cellText=rows, colLabels=col_labels, cellLoc='center', loc='center',
              colWidths=[.14,.16,.22,.26,.20])
tb.auto_set_font_size(False); tb.set_fontsize(9.5)
for (i,j),cl in tb.get_celld().items():
    cl.set_height(0.17 if i==0 else 0.13)
    if i==0: cl.set_facecolor('#1a3a5c'); cl.set_text_props(color='white', weight='bold')
ax.set_title('Table 2.  Faraday rotation vs PCI phase shift, same detunings as Table 1  (kappa_F=1)',
             fontsize=10.5, pad=10, loc='left')
plt.tight_layout(); plt.show()
print(f"theta_F and phi share the same delta/(1+delta^2) lineshape -> the ratio column is exactly "
      f"kappa_F={kappa_F} at every detuning (a direct check on the implementation above).")
