# 09 Multi-shot Simulation
#
# Exported from: 1 calculations revised 2  multishot  6  extended.ipynb
# This file is a mechanical notebook-section export for refactoring.
# Keep physics, numerical algorithms, parameter values, and execution order equivalent to the notebook.


# %% [cell 39: markdown]
# ---
# # Part II — The multi-shot run
#
# Sections 1–12 sized a **single** frame and the closed-form budget $N_\mathrm{max}$. But a
# non-destructive run is a *sequence* of frames on the *same* atoms, and that is where the interesting
# physics lives. Part II turns the budget into an explicit frame-by-frame simulation and uses it to
# answer the questions an experimentalist actually asks: how fast does the condensate fade, what does
# the camera see, which detuning to choose, how many atoms are needed, and how to resolve the cloud.
#
# ## 13. The multi-shot sequence engine
#
# Every probe pulse scatters $N_\gamma=R_\mathrm{scatt}\,\tau$ photons per atom, so the condensate is
# progressively depleted: the Thomas–Fermi cloud shrinks ($R_i\propto N_0^{1/5}$) and the dispersive
# phase fades ($\varphi\propto$ column density $\propto N_0^{3/5}$). **Each successive frame is dimmer
# and smaller than the last**, and the run ends when the condensate has lost the tolerated fraction.
#
# Two models bracket the truth (§12):
#
# * **heating (realistic).** $E_\mathrm{rec}=359$ nK is below the ODT trap depth, so a recoiled atom
#   stays trapped, thermalises, and *heats* a cloud that is only $\sim$22 % condensed; the condensate
#   *melts* via $N_0/N_\mathrm{tot}=1-(T/T_c)^3$ as $T$ climbs through
#   $E(T)=3(\zeta_4/\zeta_3)k_BT(T/T_c)^3$.
# * **clean loss (optimistic).** Each recoil ejects its atom: $N_0(s)=N_0\,e^{-\eta N_\gamma s}$.

# %% [cell 40: code]
# ============================================================================
# 13. MULTI-SHOT SEQUENCE ENGINE
# ============================================================================
#
#  >>>>>>>>>>>>>>>  USER CONTROLS -- EDIT THESE AND RE-RUN 13-17  <<<<<<<<<<<<<<
SEQ_Delta_GHz = 1.5     # detuning                (GHz)   try 1.0, 1.5, 2.5, 3.0
SEQ_P_mW      = 3.5     # probe power / intensity  (mW)   try 2.0, 3.5, 5.0
SEQ_tau_us    = 40.0    # pulse duration           (us)   try 15, 40, 80
SEQ_axis      = 0       # 0 = x across-cigar (resolved w/ TOF),
                        # 1 = y along-cigar (atom-number probe, phase-wrapped @1.5GHz),
                        # 2 = z across-cigar
SEQ_loss_frac = 0.30    # stop run at this CONDENSATE-loss fraction (protocol = 0.30)
SEQ_reabs     = None    # None -> realistic reabs_frac(Delta);  set 0.0 for recoil-only heating
SEQ_max_shots = 400     # safety cap on number of frames
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
SEQ_Delta, SEQ_tau = SEQ_Delta_GHz * 1e9, SEQ_tau_us * 1e-6


def tf_state(N0_now):
    """Thomas-Fermi mu, n_peak, R(x,y,z), peak column density(x,y,z) for a given
    condensate number -- identical formulas to section 4, just re-evaluated as
    atoms are lost."""
    mu_   = 0.5 * (15 * N0_now * a_s / a_ho)**(2/5) * hbar * omega_bar
    npk   = mu_ * m / (4 * np.pi * hbar**2 * a_s)
    R_    = np.sqrt(2 * mu_ / (m * omega**2))
    ncol_ = (4/3) * npk * R_
    return mu_, npk, R_, ncol_


def SNR_pixel_phi(phi, P_mW, axis, tau_s, blur):
    """Realistic per-pixel PCI SNR for an explicit phase value (full transfer
    model, peak intensity, plate-t_p^2 background, fixed axis NA blur).
    Returns NaN if phi >= pi (phase-wrapped: no single-valued linear SNR)."""
    if not (0.0 < phi < np.pi):
        return np.nan
    Iimg = bg_plate + blur * (I_full(phi) - bg_plate)
    Nph  = N_phot_pix(P_mW, tau_s)
    return abs(Iimg - bg_plate) * Nph / np.sqrt(abs(Iimg) * Nph + read_e**2)


def run_sequence(Delta_Hz, P_mW, tau_s, axis=0, loss_frac=0.30,
                 model='heating', reabs=None, max_shots=400):
    """Step through a continuous run frame-by-frame. Frame s sees the cloud
    AFTER s previous probe pulses of destruction (frame 0 = fresh cloud).

    Returns a dict of per-frame arrays: shot index, condensate number N0,
    condensate-loss fraction, temperature T (K, heating only), peak phase, and
    realistic per-pixel SNR. The frame counts reproduce Nmax_heating /
    Nmax_cleanloss (defined in section 12) by construction.
    """
    Ng   = N_scatt(Delta_Hz, P_mW, tau_s)         # photons/atom/shot (cloud-independent)
    blur = blur_axis[axis]                         # ~const: R changes <10% out to 30% loss
    N0_0 = N0                                       # fresh condensate number
    # --- heating bookkeeping: SAME coefficients as Nmax_heating() in section 12 ---
    Tc = Tc_sc
    if reabs is None:
        reabs = reabs_frac(Delta_Hz)
    A_E = 3 * (zeta4 / zeta3) * kB / Tc**3          # per-atom energy coeff: E(T) = A_E * T^4
    dE  = Ng * (1 + reabs) * E_rec                  # per-atom energy deposited per shot
    T   = T_cloud

    out = {key: [] for key in ('shot', 'N0', 'frac', 'T', 'phi', 'snr')}
    s = 0
    while s <= max_shots:
        if model == 'heating':
            N0_now = N_tot_sc * (1 - (T / Tc)**3)
            T_now  = T
        else:                                       # clean-loss (optimistic)
            N0_now = N0_0 * np.exp(-eta_coll * Ng * s)
            T_now  = np.nan
        if N0_now <= 0:
            break
        frac = 1 - N0_now / N0_0
        _, _, _, ncol = tf_state(N0_now)
        phi = phi_peak(Delta_Hz, ncol[axis])
        out['shot'].append(s);   out['N0'].append(N0_now);  out['frac'].append(frac)
        out['T'].append(T_now);  out['phi'].append(phi)
        out['snr'].append(SNR_pixel_phi(phi, P_mW, axis, tau_s, blur))
        if frac >= loss_frac:
            break
        T = (T**4 + dE / A_E)**0.25                 # advance heating by one shot
        s += 1
    return {key: np.asarray(val) for key, val in out.items()}


def accumulate(snr):
    """RMS-accumulated SNR vs frame: sqrt(cumsum(SNR_s^2)). NaN frames (phase-
    wrapped) contribute nothing."""
    return np.sqrt(np.nancumsum(np.where(np.isfinite(snr), np.asarray(snr, float)**2, 0.0)))


# ---- run the user-selected operating point under BOTH models ----
seq_h = run_sequence(SEQ_Delta, SEQ_P_mW, SEQ_tau, SEQ_axis, SEQ_loss_frac, 'heating', SEQ_reabs, SEQ_max_shots)
seq_l = run_sequence(SEQ_Delta, SEQ_P_mW, SEQ_tau, SEQ_axis, SEQ_loss_frac, 'loss',    SEQ_reabs, SEQ_max_shots)

sat  = intensity_at_atoms(SEQ_P_mW) / Isat
Ng0  = N_scatt(SEQ_Delta, SEQ_P_mW, SEQ_tau)
fh   = Nmax_heating(SEQ_Delta, SEQ_P_mW, SEQ_loss_frac, SEQ_tau)
fl   = Nmax_cleanloss(SEQ_Delta, SEQ_P_mW, SEQ_loss_frac, SEQ_tau)
reab = SEQ_reabs if SEQ_reabs is not None else reabs_frac(SEQ_Delta)

print("="*70)
print(f"CONTINUOUS RUN @  Delta={SEQ_Delta_GHz} GHz   P={SEQ_P_mW} mW   tau={SEQ_tau_us:.0f} us"
      f"   axis={'xyz'[SEQ_axis]}")
print("-"*70)
print(f"  intensity at atoms (peak)  : I/Isat = {sat:.3f}")
print(f"  photons / atom / shot      : N_gamma = {Ng0:.2e}   (reabsorption factor {reab:.3f})")
print(f"  initial peak phase         : phi0 = {seq_h['phi'][0]:.3f} rad   ({regime_of(seq_h['phi'][0])})")
print(f"  stop run at                : {SEQ_loss_frac:.0%} condensate loss")
print("-"*70)
print(f"  usable frames to {SEQ_loss_frac:.0%} loss   :  HEATING (realistic) ~{fh:.0f}"
      f"   |   LOSS (optimistic) ~{fl:.0f}")
print(f"  total non-destructive time :  {fh*SEQ_tau_us:.0f} us (heating)   |   {fl*SEQ_tau_us:.0f} us (loss)")
print(f"  Nmax x tau (set by Delta,P):  {fh*SEQ_tau_us:.0f} us   (invariant -- tau only re-slices it)")
print("="*70)

# %% [cell 41: markdown]
# ## 14. Evolution of the run
#
# Four views of the same sequence; the realistic **heating** curve (red) and the optimistic
# **clean-loss** curve (blue) bracket every quantity.
#
# * **(a) Condensate depletion.** Surviving condensate vs frame; the right axis (heating) shows $T$
#   climbing toward $T_c$ — the mechanism, since the condensate fraction $1-(T/T_c)^3$ collapses as
#   $T\!\to\!T_c$, melting the cloud in roughly *half* the frames of clean loss.
# * **(b) Phase fades.** $\varphi\propto N_0^{3/5}$ shrinks frame-to-frame; the orange line is the
#   linear-PCI ceiling $\varphi=0.5$.
# * **(c) Per-frame SNR** declines as the cloud dims and shrinks.
# * **(d) Accumulated SNR** $=\sqrt{\sum_s\mathrm{SNR}_s^2}$ **saturates**: once the condensate is gone
#   there is no more signal to add. The plateau value — not the per-shot SNR — is the figure of merit
#   for a fixed destruction budget, and adding frames beyond $N_\mathrm{max}$ buys nothing.

# %% [cell 42: code]
# ---- 14. Evolution of the run: condensate, temperature, phase, SNR ----------
fig, axs = plt.subplots(2, 2, figsize=(12.5, 8.4), constrained_layout=True)

# (a) condensate surviving (%) vs frame, with T (nK) on a twin axis (heating) -----
axA = axs[0, 0]
survH = 100 * (1 - seq_h['frac'])
survL = 100 * (1 - seq_l['frac'])
axA.plot(seq_h['shot'], survH, 'o-', color='#c4161c', lw=2.2, ms=4, label='heating (realistic)')
axA.plot(seq_l['shot'], survL, 's--', color='#1f5fa8', lw=2.0, ms=3, label='loss (optimistic)')
axA.axhline(100 * (1 - SEQ_loss_frac), color='k', ls=':', lw=1.1)
axA.annotate(f'{SEQ_loss_frac:.0%} loss budget', (0.5, 100 * (1 - SEQ_loss_frac) + 1.2), fontsize=8.5)
axA.set_xlabel('frame number  s')
axA.set_ylabel('condensate surviving  (% of initial N0)')
axA.set_title('(a) Condensate depletion over the run')
axA.grid(alpha=0.25)
axA.legend(fontsize=8.6, loc='lower left')
axT = axA.twinx()
axT.plot(seq_h['shot'], seq_h['T'] * 1e9, ':', color='#c4161c', lw=1.4, alpha=0.7)
axT.axhline(Tc_sc * 1e9, color='gray', ls='--', lw=0.9)
axT.annotate('T_c', (seq_h['shot'][-1] * 0.88, Tc_sc * 1e9 - 7), color='gray', fontsize=8.5)
axT.set_ylabel('T (nK) -- heating only', color='#c4161c')
axT.tick_params(axis='y', colors='#c4161c')

# (b) peak phase vs frame -----------------------------------------------------
axB = axs[0, 1]
axB.plot(seq_h['shot'], seq_h['phi'], 'o-', color='#c4161c', lw=2.2, ms=4, label='heating')
axB.plot(seq_l['shot'], seq_l['phi'], 's--', color='#1f5fa8', lw=2.0, ms=3, label='loss')
_phi_top = max(0.55, np.nanmax(np.concatenate([seq_h['phi'], seq_l['phi']])) * 1.15)
axB.set_ylim(0.0, _phi_top)
axB.axhline(0.5, color='#e08020', ls=':', lw=1.2)
axB.text(seq_l['shot'][-1] * 0.30, 0.5, 'linear-PCI ceiling  phi = 0.5',
         color='#b06010', fontsize=8.3, va='bottom', ha='left')
axB.set_xlabel('frame number  s')
axB.set_ylabel('peak phase  phi  (rad)')
axB.set_title('(b) Dispersive signal fades as the cloud shrinks')
axB.grid(alpha=0.25)
axB.legend(fontsize=8.6)

# (c) per-shot SNR vs frame ---------------------------------------------------
axC = axs[1, 0]
axC.plot(seq_h['shot'], seq_h['snr'], 'o-', color='#c4161c', lw=2.2, ms=4, label='heating')
axC.plot(seq_l['shot'], seq_l['snr'], 's--', color='#1f5fa8', lw=2.0, ms=3, label='loss')
axC.set_xlabel('frame number  s')
axC.set_ylabel('per-shot PCI SNR / pixel')
axC.set_title('(c) Per-frame SNR declines shot-to-shot')
axC.grid(alpha=0.25)
axC.legend(fontsize=8.6)

# (d) accumulated SNR vs frame ------------------------------------------------
axD = axs[1, 1]
accH = accumulate(seq_h['snr'])
accL = accumulate(seq_l['snr'])
axD.plot(seq_h['shot'], accH, 'o-', color='#c4161c', lw=2.2, ms=4,
         label=f'heating  ->  {accH[-1]:.0f} total')
axD.plot(seq_l['shot'], accL, 's--', color='#1f5fa8', lw=2.0, ms=3,
         label=f'loss  ->  {accL[-1]:.0f} total')
axD.set_xlabel('frame number  s')
axD.set_ylabel('accumulated SNR = sqrt( sum SNR_s^2 )')
axD.set_title('(d) Accumulated SNR saturates as atoms run out')
axD.grid(alpha=0.25)
axD.legend(fontsize=8.6, loc='lower right')

fig.suptitle(f"Continuous run @ Delta={SEQ_Delta_GHz} GHz, P={SEQ_P_mW} mW, tau={SEQ_tau_us:.0f} us, "
             f"axis {'xyz'[SEQ_axis]}  --  heating melts the condensate in ~half the frames of clean loss",
             fontsize=11)
plt.savefig('fig_sequence_evolution.png', dpi=140, bbox_inches='tight', facecolor='white')
plt.show()

# %% [cell 43: markdown]
# ## 15. The cloud fading on camera
#
# The curves are abstract; here is what the *camera* sees. Using the full Fourier-optics chain of §7
# (finite-NA pupil, phase plate on the carrier, binning to camera pixels, $\mathrm{QE}=0.40$, $7\,e^-$
# read noise), we render the PCI frame at the **start, middle and end** of the heating run, each with
# its *current* depleted condensate — so both the phase contrast and the Thomas–Fermi radii shrink as
# the run proceeds. The dispersive PCI peak visibly flattens and narrows: atom loss made concrete.

# %% [cell 44: code]
# ---- 15. The cloud fading on camera: PCI frames across the run (heating model)
def sim_image_state(axis, phi_val, R_state, mode='PCI', OD=4.0):
    """Same Fourier-optics propagation as section 7's sim_image(), but with the
    Thomas-Fermi radii supplied explicitly so the cloud can SHRINK as N0 falls."""
    plane = [i for i in range(3) if i != axis]
    prof  = np.maximum(0, 1 - GA**2 / R_state[plane[0]]**2 - GB**2 / R_state[plane[1]]**2)**1.5
    Esc   = np.fft.ifft2(np.fft.fft2(np.exp(1j * phi_val * prof) - 1) * pupil)
    if   mode == 'PCI':  E = t_p * np.exp(1j * theta) + Esc
    elif mode == 'DGI':  E = 10**(-OD/2) + Esc
    else:                E = 1 + Esc
    return np.abs(E)**2


def camera_from_image(Iratio, P_mW, tau_s):
    """Bin to camera pixels (15 grid cells = one object pixel) and add Poisson
    photon noise + Gaussian read noise -- identical recipe to section 7.4."""
    nb = (Ngrid // 15) * 15
    binned = Iratio[:nb, :nb].reshape(nb // 15, 15, nb // 15, 15).mean(axis=(1, 3))
    Nd = N_phot_pix(P_mW, tau_s)
    counts = rng.poisson(np.clip(binned, 0, None) * Nd) + rng.normal(0, read_e, binned.shape)
    return counts / Nd, binned


# pick three frames along the heating run: start, middle, last
Slast = len(seq_h['shot']) - 1
frame_idx = sorted(set([0, max(1, Slast // 2), Slast]))
ext = [-FOV/2 * 1e6, FOV/2 * 1e6, -FOV/2 * 1e6, FOV/2 * 1e6]
nbc = Ngrid // 15
ycam = (np.arange(nbc) - nbc // 2 + 0.5) * pix_obj * 1e6
mid = nbc // 2

ncols = len(frame_idx)
fig, axs = plt.subplots(2, ncols, figsize=(4.6 * ncols, 6.8), constrained_layout=True)
if ncols == 1:
    axs = axs.reshape(2, 1)

for col, fi in enumerate(frame_idx):
    _, _, Rst, _ = tf_state(seq_h['N0'][fi])
    phi_fi = seq_h['phi'][fi]
    Iimg = sim_image_state(SEQ_axis, phi_fi, Rst, 'PCI')
    cam, binned = camera_from_image(Iimg, SEQ_P_mW, SEQ_tau)

    im = axs[0, col].imshow(cam, extent=ext, origin='lower', cmap='inferno', vmin=0.85, vmax=1.15)
    axs[0, col].set_xlim(-45, 45)
    axs[0, col].set_ylim(-12, 12)
    axs[0, col].set_title(f"frame {fi}:  {100*(1-seq_h['frac'][fi]):.0f}% condensate left\n"
                          f"phi = {phi_fi:.3f} rad,  N0 = {seq_h['N0'][fi]:.0f}")
    axs[0, col].set_xlabel('y (um)')
    axs[0, col].set_ylabel('z (um)')
    plt.colorbar(im, ax=axs[0, col], fraction=0.032, pad=0.02, label='I/I0')

    axs[1, col].plot(ycam, binned[mid], color='#1f5fa8', lw=2, label='noiseless')
    axs[1, col].plot(ycam, cam[mid], '.', color='#c4161c', ms=4, alpha=0.7, label='single shot + noise')
    axs[1, col].axhline(bg_plate, color='gray', ls=':', lw=1, label='plate background t_p^2')
    axs[1, col].set_xlim(-45, 45)
    axs[1, col].set_ylim(0.82, 1.18)
    axs[1, col].set_xlabel('y (um)')
    axs[1, col].set_ylabel('I/I0')
    axs[1, col].grid(alpha=0.25)
    if col == 0:
        axs[1, col].legend(fontsize=8)

fig.suptitle("PCI camera frames at the start, middle and end of the run -- the dispersive PCI peak "
             "shrinks and flattens as the condensate is depleted", fontsize=11)
plt.savefig('fig_sequence_frames.png', dpi=140, bbox_inches='tight', facecolor='white')
plt.show()

# %% [cell 45: markdown]
# ## 16. Detuning sweep: the SNR–destruction trade-off
#
# This is the central design trade-off of dispersive imaging, and it reproduces the structure of K24
# Fig. 4. **Figure 1** establishes *which* imaging regime applies and *where* the usable window sits;
# **Figure 2** runs the full detuning sweep as a 100-image-per-detuning experiment.
#
# **The sweep (Fig. 2).** This study deliberately uses a short 15 µs pulse — shorter than the 40 µs operating baseline of §11–§12 — so the low-detuning destruction plays out *within* the 100-frame window and the gradient is visible. For each detuning we run a fixed 100-image sequence (heating model) and record
# the per-image SNR; stacking the rows gives the SNR-vs-(detuning, image) map, and four representative
# detunings are rendered as strips of consecutive in-trap PCI frames. The result is the Fig.-4 picture:
# a **bright but instantly-destroyed tongue at low detuning** (large signal, tiny $N_\mathrm{max}$),
# grading into **dim signal that survives all 100 frames at high detuning** (small signal, large
# $N_\mathrm{max}$). The accumulated SNR is nearly the same along the sweep — §12's invariance made
# visible. *Practical reading:* pick low detuning for a few high-contrast frames of fast dynamics, high
# detuning for a long, gently-sampled run. (SNR is single-shot per-pixel, not the region-integrated,
# multi-set value of the source figure.)

# %% [cell 46: code]
# ============================================================================
# 16. DETUNING SWEEP: THE SNR-DESTRUCTION TRADE-OFF  (cf. K24 Fig. 4)
#     + the absorption<->dispersion crossover that bounds the usable window
# ============================================================================

#  >>>>>  controls for this study (independent of the section-13 controls)  <<<<<
MAP_P_mW   = 3.5      # probe power (mW) -- low/short so the low-detuning end dies fast
MAP_tau_us = 40.0     # pulse duration (us)
MAP_axis   = 0        # across cigar (resolves the transverse cigar shape)
N_IMAGES   = 60      # images per detuning (as in the paper)
strip_dets = [0.5, 1.0, 1.5, 2.0]   # detunings shown as image strips (GHz)
n_strip    = 15       # consecutive frames per strip
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
MAP_tau = MAP_tau_us * 1e-6
blurM   = blur_axis[MAP_axis]
nc0     = n_col[MAP_axis]            # fresh-cloud peak column density (condensate)

# ---- unified single-shot imaging SNR, valid from resonance to far detuning ----
def snr_abs_pix(D, ncol, P, tau_s):
    """Resonant-absorption-imaging SNR: signal is the optical density, with shot+read
    noise on both the with-atom frame (transmits e^-OD) and the light frame."""
    OD = sigma0 * ncol / (1 + delta_of(D)**2)
    N  = N_phot_pix(P, tau_s); Na = N * np.exp(-OD)
    var = (Na + read_e**2) / Na**2 + (N + read_e**2) / N**2
    return OD / np.sqrt(var), OD

def snr_pci_pix(D, ncol, P, tau_s, blur):
    """Dispersive PCI SNR with absorption folded into the complex transmission
    t = exp(i*phi - OD/2) -- so it stays correct as the probe approaches resonance."""
    d = delta_of(D); phi = sigma0 * ncol / 2 * d / (1 + d**2); OD = sigma0 * ncol / (1 + d**2)
    I = abs(t_p * np.exp(1j*theta) + (np.exp(1j*phi - OD/2) - 1))**2
    I = bg_plate + blur * (I - bg_plate); N = N_phot_pix(P, tau_s)
    return abs(I - bg_plate) * N / np.sqrt(abs(I) * N + read_e**2), phi

def usable_snr(D, ncol, P, tau_s, blur):
    """Best achievable single-shot SNR at this detuning: PCI where it is single-valued
    (phi < pi), otherwise resonant absorption."""
    sa, OD  = snr_abs_pix(D, ncol, P, tau_s)
    sp, phi = snr_pci_pix(D, ncol, P, tau_s, blur)
    return max(sp if phi < np.pi else 0.0, sa), phi, OD, sa, sp, (phi < np.pi)

# ---- frame-0 signal vs detuning, decomposed (the crossover) ----
D_ax = np.logspace(np.log10(0.04), np.log10(4.0), 220)
rows = [usable_snr(D*1e9, nc0, MAP_P_mW, MAP_tau, blurM) for D in D_ax]
best = np.array([r[0] for r in rows]); phiA = np.array([r[1] for r in rows])
ODA  = np.array([r[2] for r in rows]); saA  = np.array([r[3] for r in rows])
spA  = np.array([r[4] for r in rows]); valA = np.array([r[5] for r in rows], bool)
D_wrap = D_ax[np.argmin(np.abs(phiA - np.pi))]
D_od1  = D_ax[np.argmin(np.abs(ODA - 1.0))]
D_peak = D_ax[np.argmax(best)]

# ---- depleting per-image SNR over a fixed run (fills the map) ----
def usable_curve(D, n_images):
    """Per-image usable SNR over a fixed N-image run as the condensate depletes
    (heating model); decays to ~0 once the BEC is gone."""
    Ng = N_scatt(D, MAP_P_mW, MAP_tau); Tc = Tc_sc; A_E = 3*(zeta4/zeta3)*kB/Tc**3
    dE = Ng*(1+reabs_frac(D))*E_rec; T = T_cloud; snr = np.zeros(n_images)
    for s in range(n_images):
        fc = max(N_tot_sc*(1-(T/Tc)**3), 0.0) / N0
        snr[s] = usable_snr(D, fc*nc0, MAP_P_mW, MAP_tau, blurM)[0]
        T = (T**4 + dE/A_E)**0.25
    return snr

def make_strip(Delta_GHz, n_frames):
    """Concatenate n_frames consecutive single-shot PCI camera frames (vertical
    cigars) into one strip; the condensate depletes frame-to-frame (uses section 15)."""
    D = Delta_GHz*1e9; Ng, Tc = N_scatt(D, MAP_P_mW, MAP_tau), Tc_sc
    A_E = 3*(zeta4/zeta3)*kB/Tc**3; dE = Ng*(1+reabs_frac(D))*E_rec; T = T_cloud
    yhw = int(round(30e-6/pix_obj)); zhw = int(round(6e-6/pix_obj)); c = (Ngrid//15)//2
    sep = np.full((2*yhw, 2), np.nan); cols = []
    for s in range(n_frames):
        N0n = max(N_tot_sc*(1-(T/Tc)**3), 0.0)
        if N0n > 0:
            _, _, Rst, ncol = tf_state(N0n)
            Iimg = sim_image_state(MAP_axis, phi_peak(D, ncol[MAP_axis]), Rst, 'PCI')
        else:
            Iimg = np.full((Ngrid, Ngrid), bg_plate)
        cam, _ = camera_from_image(Iimg, MAP_P_mW, MAP_tau)
        cols.append(cam.T[c-yhw:c+yhw, c-zhw:c+zhw] - bg_plate)
        if s < n_frames-1: cols.append(sep)
        T = (T**4 + dE/A_E)**0.25
    return np.hstack(cols)

# ---- build the map (to resonance) + the four strips ----
det_lo, det_hi = 0.1, 2.0
det_axis = np.logspace(np.log10(det_lo), np.log10(det_hi), 70)
snr_map = np.array([usable_curve(D*1e9, N_IMAGES) for D in det_axis])
strips  = {D: make_strip(D, n_strip) for D in strip_dets}

# ---- FIGURE 2: the Fig.4 reproduction (map + image strips) ----
fig = plt.figure(figsize=(14, 6.6), constrained_layout=True)
gs = fig.add_gridspec(len(strip_dets), 2, width_ratios=[1.05, 1.0])
axMap = fig.add_subplot(gs[:, 0])
im = axMap.imshow(snr_map, aspect='auto', extent=[0, N_IMAGES, np.log10(det_hi), np.log10(det_lo)],
                  cmap='jet', vmin=0, vmax=np.percentile(snr_map, 99.5))
yt = [0.1, 0.2, 0.5, 1.0, 2.0]; axMap.set_yticks(np.log10(yt)); axMap.set_yticklabels([f'{v:g}' for v in yt])
for D in strip_dets: axMap.axhline(np.log10(D), color='0.8', lw=2.0, alpha=0.8)
axMap.set_xlabel('image number'); axMap.set_ylabel('detuning (GHz)')
cb = fig.colorbar(im, ax=axMap, location='top', fraction=0.05, pad=0.03, aspect=38)
cb.set_label('usable SNR / pixel (single shot)')
smax = max(np.nanpercentile(s, 99.5) for s in strips.values())
for i, D in enumerate(strip_dets):
    axS = fig.add_subplot(gs[i, 1])
    axS.imshow(strips[D], aspect='auto', cmap='jet', vmin=-0.02, vmax=smax, interpolation='nearest', origin='lower')
    axS.set_xticks([]); axS.set_yticks([])
    axS.text(0.012, 0.90, f'{int(round(D*1000))} MHz', transform=axS.transAxes, color='white',
             fontsize=10, va='top', ha='left', fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.18', fc='black', ec='none', alpha=0.45))
    if i == 0: axS.set_title(f'{n_strip} consecutive single-shot PCI frames  (time ->)', fontsize=10)
fig.suptitle('Detuning sweep of an in-trap PCI run (cf. K24 Fig. 4): a bright but instantly-destroyed tongue '
             'at low detuning,\ngrading into dim signal that survives all 100 frames at high detuning', fontsize=11)
plt.savefig('fig_detuning_map.png', dpi=140, facecolor='white'); plt.show()

# %% [cell 92: markdown]
# ### 19.6 Step 14 — the price of looking: the run as a filmstrip
#
# **Theory.** Each frame scatters $N_\gamma$ photons per atom, depositing a fixed energy that walks
# the temperature up ($T^4_{s+1} = T^4_s + \Delta E/A_E$, §12–13) and eats the condensate,
# $N_0(T) = N_\mathrm{tot}[1-(T/T_c)^3]$. The shrinking cloud imprints a shrinking $\varphi$ — the
# signal the *next* frame has to work with.
#
# **What is happening.** Below is the §13 heating run (Δ = 1.5 GHz, 3.5 mW, 40 µs) rendered as what
# it actually is: a sequence of camera frames of one condensate being watched to its 30% loss
# budget. Each panel is a full noisy PCI frame simulated at that shot's phase, with the shot's
# $N_0$, $T$ and $\varphi$ as captions. Frame by frame the signal drains — a fifth of $\varphi$
# over the run, read off the captions — because non-destructive imaging does not mean free imaging:
# it means *budgeted* imaging, and this is the budget being spent. (The Thomas–Fermi radii change by <10% over the run, so the fresh-cloud profile is reused
# with each shot's own $\varphi_s$ — the same approximation as the engine itself, §13.)
#
# **What the code does.** `run_sequence(...)` supplied the per-shot $\varphi_s$ (stored in `seq_h`);
# each panel is one `sim_image` + Poisson/read-noise draw at the run's own power and pulse length.

# %% [cell 93: code]
# ---- STEP 14: watching the condensate die, frame by frame ------------------------
st_show = [0, 5, 10, 14]                                   # shots along the heating run
st_Nd_run = N_phot_pix(SEQ_P_mW, SEQ_tau)                  # the run's own power AND pulse length

fig, axs = plt.subplots(1, len(st_show), figsize=(15.2, 3.6))
for a, s in zip(axs, st_show):
    phi_s = seq_h['phi'][s]
    I_s, _ = sim_image(SEQ_axis, phi_s, 'PCI')
    b = I_s[:_nb2,:_nb2].reshape(_nb2//15,15,_nb2//15,15).mean(axis=(1,3))
    frame = (rng.poisson(np.clip(b,0,None)*st_Nd_run) + rng.normal(0, read_e, b.shape)) / st_Nd_run
    im = a.imshow(frame, extent=ext, origin='lower', cmap='inferno', vmin=0.90, vmax=1.18)
    a.set_xlim(-45, 45); a.set_ylim(-12, 12); a.set_xlabel('y (um)')
    a.set_title(f'shot {s}\n$N_0$ = {seq_h["N0"][s]/1e3:.1f}k   T = {seq_h["T"][s]*1e9:.0f} nK'
                f'   $\\varphi$ = {phi_s:.2f}', fontsize=9.5)
axs[0].set_ylabel('z (um)')
plt.colorbar(im, ax=axs[-1], fraction=0.04, label='$I/I_0$')
fig.suptitle(f'Step 14 — the same condensate across the run  '
             f'(PCI, $\\Delta$ = {SEQ_Delta_GHz} GHz, {SEQ_P_mW} mW, {SEQ_tau_us:.0f} us; '
             f'stops at {SEQ_loss_frac:.0%} loss)', fontsize=12)
plt.tight_layout(); plt.show()
print(f"signal spent per frame: phi falls {seq_h['phi'][0]:.3f} -> {seq_h['phi'][-1]:.3f} rad over "
      f"{len(seq_h['shot'])} frames while T climbs {seq_h['T'][0]*1e9:.0f} -> {seq_h['T'][-1]*1e9:.0f} nK. "
      f"Every mode of Steps 6-13 spends this same budget; only what it buys per frame differs.")
