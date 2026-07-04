# 10 Analysis
#
# Exported from: 1 calculations revised 2  multishot  6  extended.ipynb
# This file is a mechanical notebook-section export for refactoring.
# Keep physics, numerical algorithms, parameter values, and execution order equivalent to the notebook.


# %% [cell 27: markdown]
# ## 10. Spatial resolution
# $\mathrm{NA}=(D/2)/f_1=0.08$, $\Delta x=0.61\lambda/\mathrm{NA}=3.06\,\mu$m. The across-cigar
# radii ($R_x=1.19$, $R_z=1.49\,\mu$m) are unresolved; $R_y=24.8\,\mu$m is resolved. Brief
# time-of-flight $R_i(t)=R_i(0)\sqrt{1+(\omega_i t)^2}$ restores transverse resolution.

# %% [cell 28: code]
NA_e = (D_probe/2)/f1
res_e = 0.61*lam/NA_e
res_obj = 0.61*lam/0.40
t_arr = np.linspace(0, 2e-3, 2000)
Rx_t = R[0]*np.sqrt(1+(omega[0]*t_arr)**2)
t_cross = t_arr[np.argmax(Rx_t > res_e)]
print(f"existing arm NA={NA_e:.3f}: dx = {res_e*1e6:.2f} um")
print(f"trichroic obj NA=0.40:     dx = {res_obj*1e6:.2f} um")
print(f"Rx(t) crosses dx at t = {t_cross*1e3:.2f} ms")

# depth of focus vs cloud extent along the imaging axis
DOF = lam/NA_e**2
print(f"\n[check] depth of focus  ~ lam/NA^2 = {DOF*1e6:.0f} um")
print(f"        cloud half-depth along x (across-cigar imaging) = {R[0]*1e6:.2f} um  -> deep within DOF")
print(f"        cloud half-depth along y (along-cigar imaging)  = {R[1]*1e6:.2f} um  -> within DOF ({R[1]*1e6:.0f}<{DOF*1e6:.0f}), no defocus blur")

# %% [cell 29: markdown]
# ## 11. Operating-point summary
#
# The report cell (folded in from the former §13) tabulates every derived quantity at a chosen $(\Delta,P,\text{axis})$ in both the idealised and realistic detection models.
#
# $\Delta=1.5$ GHz, $P=2$ mW, $\tau=40\,\mu$s (the balanced operating point of §12.3), across-cigar ($x$), revised conventions
# (peak intensity, $\eta=1.3$, QE = 0.40, $7\,e^-$ read, full transfer model, axis-aware
# NA blur). Numbers below are emitted by the report cell.
#
# | quantity | original | revised ($\tau=40\,\mu$s) |
# |---|---|---|
# | $N_\gamma$ /atom/shot | $1.0\times10^{-3}$ | $5.3\times10^{-3}$ (peak intensity) |
# | single-shot SNR across cigar @ 2 mW | 6.5 (QE=1) | ~5.8/pixel, ~11.7/res-elem |
# | same @ 3.5 mW (cascade) | — | ~7.6/pixel, ~15.5/res-elem |
# | **DGI peak signal @ 5 mW** | — | **~24 e⁻ — quantitative (> 3× read)** |
# | $N_\mathrm{max}$ (30% loss, clean-loss / heating) | — | ~52 / ~24 |
# | $N_\mathrm{max}$ ($\varepsilon=\mu$, 13% loss, $\eta$=1.3) | 134 | ~21 |
# | across-cigar NA blur $b_x$ | (not modelled) | 0.61 |
# | recommended detuning method | Route A | **E** injection-lock / **B** (calibration, not SNR) |
#
# **Along-cigar caveat.** Imaging *along* the cigar at 1.5 GHz is **phase-wrapped**
# ($\varphi_y=4.25$ rad), so a linear/Meppelink SNR is not meaningful there — the report
# returns "n/a" and flags the detuning needed to reach the linear regime ($\gtrsim13$ GHz,
# where $N_\gamma$ also drops $\sim80\times$). Both *transverse* radii then sit below the
# resolution limit ($b_y=0.24$), so along-cigar imaging resolves no structure; it is an
# integrated **atom-number** probe, not a structure probe. The across-cigar axis with
# brief time-of-flight remains the route to resolved in-trap profiles.

# %% [cell 30: code]
# ---- operating-point report (folded from former section 13) ----
def operating_point_report(Delta_GHz, P_mW, axis=0, tau_s=None):
    Dh = Delta_GHz*1e9; ax_n = 'xyz'[axis]
    phi = phi_peak(Dh, n_col[axis])
    reg = 'linear PCI' if phi<0.5 else ('Meppelink' if phi<np.pi else 'phase-wrapped')
    ng = N_scatt(Dh, P_mW, tau_s)
    print('='*58)
    print(f"Operating point:  Delta={Delta_GHz:.2f} GHz (delta={delta_of(Dh):.0f})  "
          f"P={P_mW:.2f} mW  tau={(tau_s or tau)*1e6:.0f} us  axis={ax_n} (R={R[axis]*1e6:.2f} um)")
    print('-'*58)
    print(f"peak phase shift          phi   = {phi:.3f} rad  ({reg}, blur_axis={blur_axis[axis]:.2f})")
    print(f"photons/atom/shot         N_g   = {ng:.2e}   (eta={eta_coll}, peak intensity={use_peak_intensity})")
    print(f"SNR/shot ideal (QE=1, linear)   = {SNR_shot_ideal(Dh,P_mW,axis,tau_s):.2f}")
    if phi < np.pi:
        print(f"SNR/shot realistic (per pixel)  = {SNR_pixel(Dh,P_mW,axis,tau_s):.2f}")
        print(f"SNR/shot realistic (res-elem)   = {SNR_reselem_sim(Dh,P_mW,axis,tau_s):.2f}")
    else:
        print(f"SNR/shot realistic              = n/a  (phi>pi: phase-wrapped, inversion ambiguous;")
        dreq = sigma0*n_col[axis]/(2*0.5); Dreq = dreq*Gamma/(2*2*np.pi)/1e9
        print(f"                                       linear PCI needs Delta>~{Dreq:.1f} GHz, or full Meppelink / TOF)")
    print(f"N_max  eps=mu (=13% loss)       = {Nmax_loss(Dh,P_mW,T_mu/T_rec,tau_s):.0f}")
    print(f"N_max  30% loss (recommended)   = {Nmax_loss(Dh,P_mW,0.30,tau_s):.0f}")
    print(f"N_max  50% loss                 = {Nmax_loss(Dh,P_mW,0.50,tau_s):.0f}")
    print('='*58)


print("Reference summary (revised), across cigar (x):")
operating_point_report(1.5, 2.0, axis=0)
print("\nAlong cigar (y) at 1.5 GHz -- phase-wrapped, see guard:")
operating_point_report(1.5, 2.0, axis=1)
print("\nAlong cigar (y) detuned into the linear regime (Delta = 13 GHz):")
operating_point_report(13.0, 2.0, axis=1)

# %% [cell 31: markdown]
# ## 12. Pulse duration and the $(\Delta, P, \tau)$ operating point [new]
#
# Sections 6–11 fixed the reference $\tau=40\,\mu$s. Folding $\tau$ in changes two things at once, in opposite
# directions, and resolves both open issues from the review:
#
# * **Destruction.** $N_\gamma = R_\mathrm{scatt}\,\tau \propto \tau$, so longer pulses scatter more
#   and $N_\mathrm{max}\propto 1/\tau$ — the budget gets *smaller*, the realistic direction. We also
#   drop the over-optimistic clean-loss picture: with $E_\mathrm{rec}=359$ nK $\ll$ ODT depth
#   ($\sim\mu$K) a recoiled atom does **not** leave — it stays, thermalises, and **heats** a cloud
#   whose condensate is only $\sim$22 %, melting it in roughly half the clean-loss shots. Right-substance
#   references: Aikawa et al. PRL **108**, 210401 (2012) (Er, 401 nm imaging line); the Stuttgart Dy
#   group PRX **9**, 011051 (2019) (single-shot in-situ phase-contrast on a strongly dipolar condensate
#   at the analogous 421 nm line, cross-checked against resonant absorption to $<$10 %); Gajdacz et al.
#   RSI **84**, 083105 (2013) (destructiveness as measured fractional loss). A reabsorption term
#   $\propto1/\delta^2$ makes low detuning disproportionately worse.
# * **Detection.** Detected photons $\propto P\tau$, so longer pulses raise per-shot SNR and lift the
#   dark-ground signal above the $7\,e^-$ read-noise floor — what makes the DGI simulation quantitative.
#
# The pivot: **$N_\mathrm{max}\times\tau$ is fixed by $(\Delta,P)$** (the total non-destructive
# integration time), and the photon-limited accumulated SNR is invariant. So $\tau$ only chooses how to
# *slice* a fixed budget — short $\tau$ = many dim frames (dynamics), long $\tau$ = few bright frames
# (snapshots / quantitative DGI). $\tau$ is capped by in-trap motion ($\tau<0.1\,T_\mathrm{trap}=341\,\mu$s)
# and recoil blur ($v_\mathrm{rec}\tau<$ resolution).

# %% [cell 37: markdown]
# ### 12.3 Recommendation and revised headline numbers
#
# **Most reasonable combination (default Route A / AOM cascade): $\Delta\approx1.5$ GHz,
# $P\approx3.5$ mW, $\tau\approx40\,\mu$s** — about **14 non-destructive frames**, each with PCI
# SNR/resolution-element $\approx16$ and DGI $\approx17\,e^-$ (above the read-noise floor, so PCI *and*
# DGI are quantitative). This sits at the read-noise knee: long enough for usable DGI, short enough to
# keep $\sim$14 timepoints.
#
# * **Dynamics (more timepoints):** $\tau\approx15$–$20\,\mu$s → $\sim$37 frames, PCI-only (DGI floored).
# * **High-SNR snapshots / quantitative DGI:** $\tau\approx80\,\mu$s → $\sim$5 bright frames.
# * **More frames at the cost of phase size:** $\Delta\approx2.5$ GHz, $P\approx5$ mW, $\tau\approx40\,\mu$s
#   → $\sim$27 frames at $\varphi_x=0.12$ (needs the injection-locked or beat-locked source).
#
# **Manual updates.** The PDF manual's "$N_\mathrm{max}\approx140$ shots" figure was the optimistic
# clean-loss model at $\tau=15\,\mu$s; at the $\tau\approx40\,\mu$s operating point now used throughout
# (§11) the budget **lands in the tens of frames** (clean-loss $\sim$52, realistic heating $\sim$14–24),
# so update §1.1, §4, §12 of the manual. Two lab inputs turn the band into a
# number: the **ODT trap depth** (sets whether recoils leave or heat) and one measured **loss-vs-shots**
# point at a known $(\Delta,P,\tau)$.

# %% [cell 38: code]

# %% [cell 56: markdown]
# ---
# ## 18. From condensate to camera: the full simulation pipeline, stage by stage [new]
#
# Sections 1-17 built every *piece* of this simulation, but introduced them in the order needed for
# the derivations rather than the order light actually experiences them. This section reassembles
# the pieces into one explicit pipeline and walks a single concrete example through every stage, so
# that a "result" anywhere else in this notebook is visibly traceable back to the physical object it
# started from.
#
# | Stage | What happens | Derived in | Implemented as |
# |---|---|---|---|
# | 1. Construct the condensate | Thomas-Fermi ground state $n(\mathbf r)$ from $(N_0,a_s,\omega_i)$ | §4 | `tf_state`, `_tf_profile` |
# | 2. Imprint the optical signal | probe accumulates $\varphi(\mathbf r_\perp)$ (PCI/DGI) or $\theta_F(\mathbf r_\perp)$ (Faraday) from the column density | §5, §17.1 | `phi_peak`, `theta_F_peak` |
# | 3. Propagate through the imaging system | finite-NA Fourier optics band-limits the field | §7.3 | `sim_image`, `sim_faraday_fields` |
# | 4. Generate mode-specific contrast | Fourier-plane optic (PCI/DGI) or polarization analyzer (Faraday) turns the field into an intensity | §7.1, §17.2 | `sim_image(...,mode=...)`, `faraday_maps` |
# | 5. Detect | bin to camera pixels, apply QE, add shot + read noise | §7.4 | `to_camera` |
#
# The four subsections below walk stages 1-5 in order for one reference point, ending with all four
# imaging modes side by side.

# %% [cell 57: code]
# ============================================================================
# 18.1  STAGE 1 -- constructing the condensate (Sec. 4's Thomas-Fermi solution, revisited)
# ============================================================================
print("Stage 1 inputs -> the fresh Thomas-Fermi condensate (Sec. 4):")
print(f"  N0 = {N0:.2e} atoms,  a_s = {a_s/a0:.0f} a0,  "
      f"trap = ({trap_Hz[0]:.0f}, {trap_Hz[1]:.0f}, {trap_Hz[2]:.0f}) Hz")
print(f"  => mu/kB = {T_mu*1e9:.1f} nK,  n_peak = {n_peak*1e-6:.2e} cm^-3,  "
      f"R = ({R[0]*1e6:.2f}, {R[1]*1e6:.2f}, {R[2]*1e6:.2f}) um")

axis_demo = 0                                       # across-cigar imaging, as used from Sec. 7 on
plane = [i for i in range(3) if i != axis_demo]
xax = np.linspace(-1.5*R.max(), 1.5*R.max(), 400)

fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.8, 4.6))
for i, (lab, col) in enumerate(zip(['x (across)', 'y (along)', 'z (across)'], ['C0','C1','C2'])):
    n1d = n_peak*np.maximum(0, 1 - (xax/R[i])**2)   # 3D TF density cut through the centre, along axis i
    axA.plot(xax*1e6, n1d*1e-6, col, lw=2, label=f'{lab}, R={R[i]*1e6:.2f} um')
axA.set_xlabel('position (um)'); axA.set_ylabel('n(r)  (cm^-3)')
axA.set_title('(a) 3D density cuts through the trap centre'); axA.legend(fontsize=8.5); axA.grid(alpha=0.25)

col_dens_map = n_col[axis_demo] * _tf_profile(R[plane[0]], R[plane[1]])
im = axB.imshow(col_dens_map*1e-4, extent=ext, origin='lower', cmap='viridis')
axB.set_xlim(-45, 45); axB.set_ylim(-12, 12)
axB.set_xlabel(f'{"xyz"[plane[0]]} (um)'); axB.set_ylabel(f'{"xyz"[plane[1]]} (um)')
axB.set_title('(b) Column density along '+"xyz"[axis_demo]+'\n'
              '(what every imaging mode actually integrates over)')
plt.colorbar(im, ax=axB, fraction=0.032, label='n_col (cm^-2)')
fig.suptitle('Stage 1: from trap parameters to a Thomas-Fermi condensate', y=1.04, fontsize=11.5)
plt.tight_layout(); plt.show()

# %% [cell 58: markdown]
# ### 18.2 Stage 2 -- the optical signal imprinted on the probe
# Both signals below are generated from the *identical* Stage-1 column-density map -- they are
# literally the same spatial profile as panel (b) above, just multiplied by a different peak value
# ($\varphi_\mathrm{peak}$ from §5.1's scalar dispersion, or $\theta_{F,\mathrm{peak}}$ from
# §17.1's vector dispersion). This is the moment the atoms actually "write" information onto the
# probe beam; everything from here on is optics and detection engineering, not atomic physics.

# %% [cell 59: code]
# ============================================================================
# 18.2  STAGE 2 -- the optical signal the probe picks up from that condensate
# ============================================================================
D_demo, P_demo = 1.5e9, 5.0     # single reference point carried through the rest of Sec. 18
phi_demo   = phi_peak(D_demo, n_col[axis_demo])
theta_demo = theta_F_peak(D_demo, n_col[axis_demo])
phase_map = phi_demo   * _tf_profile(R[plane[0]], R[plane[1]])
rot_map   = theta_demo * _tf_profile(R[plane[0]], R[plane[1]])

fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.8, 4.6))
im0 = axA.imshow(phase_map, extent=ext, origin='lower', cmap='magma')
axA.set_title(f'(a) Scalar phase  phi(r) -- PCI/DGI read this\npeak = {phi_demo:.3f} rad')
plt.colorbar(im0, ax=axA, fraction=0.032, label='phi (rad)')
im1 = axB.imshow(rot_map, extent=ext, origin='lower', cmap='PuOr')
axB.set_title(f'(b) Rotation  theta_F(r) -- Faraday reads this\n'
              f'peak = {theta_demo:.3f} rad  (kappa_F={kappa_F})')
plt.colorbar(im1, ax=axB, fraction=0.032, label='theta_F (rad)')
for a in (axA, axB):
    a.set_xlim(-45,45); a.set_ylim(-12,12); a.set_xlabel('y (um)'); a.set_ylabel('z (um)')
fig.suptitle(f'Stage 2: same cloud, two different optical signals imprinted at '
             f'Delta={D_demo/1e9:.1f} GHz', y=1.04, fontsize=11.5)
plt.tight_layout(); plt.show()
print("Both maps share the identical spatial shape (Stage 1's column-density profile) -- only the "
      "physical origin and the size of the peak value differ (Sec. 5 vs Sec. 17.1).")

# %% [cell 60: markdown]
# ### 18.3 Stages 3-4 -- propagation and mode-specific contrast
# The phase/rotation maps of Stage 2 are not yet visible as an intensity -- a pure phase or rotation
# produces no contrast on its own. Stage 3 propagates the field through the finite-NA 4f imaging
# system, band-limiting high spatial frequencies (§7.3); Stage 4 is where the four modes
# diverge, converting that band-limited field into an intensity by four different routes (a
# phase-shifted reference plate, an opaque beam block, a crossed polarizer, or a polarizing beam
# splitter). All four panels below use one common reference probe power and detuning purely so the
# comparison is visual and fair; §17.3 already compared SNR mode-by-mode at each mode's own
# typical operating power.

# %% [cell 61: code]
# ============================================================================
# 18.3  STAGES 3-4 -- NA-limited propagation, then mode-specific contrast generation
# ============================================================================
I_pci_nl, _ = sim_image(axis_demo, phi_demo, 'PCI')
I_dgi_nl, _ = sim_image(axis_demo, phi_demo, 'DGI')
fmap_demo = faraday_maps(D_demo, axis_demo)
I_dark_nl = fmap_demo['I_dark']
S_nl = (fmap_demo['I_v'] - fmap_demo['I_u']) / (fmap_demo['I_v'] + fmap_demo['I_u'])

panels = [(I_pci_nl,  'PCI',                'inferno', (0.85, 1.15), 'I/I0'),
          (I_dgi_nl,  'DGI',                'inferno', (-0.01, 0.05), 'I/I0'),
          (I_dark_nl, 'Dark-field Faraday', 'inferno', (-0.005, 0.05), 'I_dark/I0'),
          (S_nl,      'Dual-port Faraday',  'RdBu_r',  (-0.5, 0.5), 'S')]
fig, axs = plt.subplots(1, 4, figsize=(16.5, 4.0))
for a, (img, title, cmap, vv, lab) in zip(axs, panels):
    im = a.imshow(img, extent=ext, origin='lower', cmap=cmap, vmin=vv[0], vmax=vv[1])
    a.set_xlim(-45,45); a.set_ylim(-12,12); a.set_xlabel('y (um)')
    a.set_title(title, fontsize=10.5)
    plt.colorbar(im, ax=a, fraction=0.045, label=lab)
axs[0].set_ylabel('z (um)')
fig.suptitle(f'Stage 3-4: noiseless NA-limited image-plane intensity, all four modes\n'
             f'(common reference point: Delta={D_demo/1e9:.1f} GHz, P={P_demo:.0f} mW, same cloud)',
             fontsize=11.5)
plt.tight_layout(); plt.show()
print("PCI and dual-port Faraday are LINEAR in their signal (phi, theta_F) and so visibly modulate "
      "around a nonzero background; DGI and dark-field Faraday are QUADRATIC and so appear as a "
      "faint bright sliver on an otherwise dark field -- the same linear/quadratic split as Sec. 17.3's "
      "SNR table, made visible.")

# %% [cell 62: markdown]
# ### 18.4 Stage 5 -- detection and the final result
# The last stage is the camera: bin the continuous field onto discrete pixels, convert to detected
# photons via the quantum efficiency, and add Poisson shot noise and Gaussian read noise (§7.4).
# This is the only stage that is genuinely random -- everything upstream is a deterministic function
# of the condensate and the chosen probe parameters -- which is why every SNR number in this notebook
# is a *ratio* of a deterministic signal to a well-defined noise floor, not a guess.

# %% [cell 63: code]
# ============================================================================
# 18.4  STAGE 5 -- detection: binning to camera pixels, QE, shot noise, read noise
# ============================================================================
cam_pci_nl, _  = to_camera(I_pci_nl, P_demo)
cam_dgi_nl, _  = to_camera(I_dgi_nl, P_demo)
cam_dark_nl, _ = to_camera(I_dark_nl, P_demo)
cam_u_nl, _ = to_camera(fmap_demo['I_u'], P_demo)
cam_v_nl, _ = to_camera(fmap_demo['I_v'], P_demo)
cam_S_nl = (cam_v_nl - cam_u_nl) / (cam_v_nl + cam_u_nl)

panels2 = [(cam_pci_nl,  'PCI',                'inferno', (0.85, 1.15)),
           (cam_dgi_nl,  'DGI',                'inferno', (-0.01, 0.05)),
           (cam_dark_nl, 'Dark-field Faraday', 'inferno', (-0.005, 0.05)),
           (cam_S_nl,    'Dual-port Faraday',  'RdBu_r',  (-0.5, 0.5))]
fig, axs = plt.subplots(1, 4, figsize=(16.5, 4.0))
for a, (img, title, cmap, vv) in zip(axs, panels2):
    im = a.imshow(img, extent=ext, origin='lower', cmap=cmap, vmin=vv[0], vmax=vv[1])
    a.set_xlim(-45,45); a.set_ylim(-12,12); a.set_xlabel('y (um)')
    a.set_title(title, fontsize=10.5)
    plt.colorbar(im, ax=a, fraction=0.045)
axs[0].set_ylabel('z (um)')
fig.suptitle(f'Stage 5: the actual camera frame each mode delivers, single shot\n'
             f'(same reference point as Stage 3-4: Delta={D_demo/1e9:.1f} GHz, P={P_demo:.0f} mW)',
             fontsize=11.5)
plt.tight_layout(); plt.savefig('fig_pipeline_final_frames.png', dpi=140, bbox_inches='tight', facecolor='white'); plt.show()

# ---- summary comparison, pulling every stage together ----
Nmax_demo   = Nmax_heating(D_demo, P_demo, 0.30)
s_pci_demo  = SNR_pixel(D_demo, P_demo, axis_demo)
Ndgi_demo   = I_dgi_nl.max()*N_phot_pix(P_demo)
s_dgi_demo  = Ndgi_demo/np.sqrt(Ndgi_demo + read_e**2)
s_dark_demo = SNR_faraday_sim(D_demo, P_demo, axis_demo, scheme='dark')
s_dual_demo = SNR_faraday_sim(D_demo, P_demo, axis_demo, scheme='dual')

fig, ax = plt.subplots(figsize=(10.8, 2.6)); ax.axis('off')
rows = [
    ['PCI',                'scalar phase phi',  'phase-plate carrier',     'linear',    f'{s_pci_demo:.1f}'],
    ['DGI',                'scalar phase phi',  'beam-block (dark gnd)',   'quadratic', f'{s_dgi_demo:.1f}'],
    ['Dark-field Faraday', 'rotation theta_F',  'crossed polarizer',       'quadratic', f'{s_dark_demo:.1f}'],
    ['Dual-port Faraday',  'rotation theta_F',  'PBS @ 45deg, 2 ports',    'linear',    f'{s_dual_demo:.1f}'],
]
tb = ax.table(cellText=rows,
              colLabels=['Mode','Observable','Contrast mechanism','Order','SNR/pixel (this point)'],
              cellLoc='center', loc='center', colWidths=[.20,.18,.28,.14,.20])
tb.auto_set_font_size(False); tb.set_fontsize(9.3)
for (i,j),cl in tb.get_celld().items():
    cl.set_height(0.22 if i==0 else 0.18)
    if i==0: cl.set_facecolor('#1a3a5c'); cl.set_text_props(color='white', weight='bold')
ax.set_title(f'Summary: all four modes share N_max={Nmax_demo:.0f} frames (30% loss budget) at '
             f'Delta={D_demo/1e9:.1f} GHz, P={P_demo:.0f} mW -- only the SNR each extracts differs',
             fontsize=10, pad=14, loc='left')
plt.tight_layout(); plt.show()

# %% [cell 64: markdown]
# Reading the five stages back to front: a headline SNR number anywhere in this notebook is the
# size of a deterministic optical signal (Stage 2, set by how many atoms are in the beam and how far
# detuned the probe is) after it has survived being band-limited by a finite lens (Stage 3),
# converted to an intensity by a specific piece of contrast-generating hardware (Stage 4), and
# finally weighed against a noise floor set by the number of photons collected and the camera's read
# noise (Stage 5). Every simulated image in §7, §15 and §17 is one particular instance
# of exactly this pipeline; §9/§12's destruction budget is simply the constraint on how many
# times Stage 2 can be repeated before the condensate itself is gone.

# %% [cell 65: markdown]
# ---
# # Part IV — One shot through the machine
#
# ## 19. The simulation, step by step [new]
#
# Everything below follows **one probe shot** through the imaging system, in the exact order the
# simulation computes it. Every figure shows the *actual state of the simulated field at that step* —
# no summary plots, no sweeps: what you see in Step $n$ is literally the array the code holds after
# executing the line quoted in Step $n$. Each step carries three short threads:
#
# > **Theory** — the equation governing this step.
# > **What is happening** — the physical event, in words.
# > **What the code does** — the line(s) of this notebook that execute it.
#
# The shot is the notebook's reference point throughout: fresh condensate, probe along $x$
# (across-cigar), $\Delta = 1.5$ GHz, $P = 5$ mW — so $\varphi_\mathrm{peak} = \theta_{F,\mathrm{peak}} = 0.203$ rad.
#
# **The plot of the story.** Steps 1–5 are shared by every mode and end in an apparent dead end: even
# through a perfect lens system, the atoms are *invisible*. Steps 6–13 are the four ways out — PCI,
# DGI, dark-field Faraday, dual-port Faraday — each one picking up the same field where Step 5 left
# it. Step 14 closes with the price of looking: the multi-shot run as a filmstrip.

# %% [cell 66: markdown]
# ### 19.1 Steps 1–5: the shared road to invisibility
#
# #### Step 1 — the atoms hand us a phase map, and nothing else
#
# **Theory.** Far from resonance the cloud's susceptibility is almost purely real (§5–6): a photon
# crossing it is *delayed*, not absorbed. The accumulated phase is proportional to the column
# density, $\varphi(\mathbf r) = \varphi_\mathrm{peak}\cdot\big(1 - y^2/R_y^2 - z^2/R_z^2\big)^{3/2}$,
# the integrated Thomas–Fermi profile of §4, with $\varphi_\mathrm{peak} = 0.203$ rad at this
# detuning.
#
# **What is happening.** The condensate acts as a weak, cloud-shaped pane of glass 50 µm wide and
# 3 µm tall. This phase map is the *only* imprint the atoms make on the light — every method in this
# notebook is a different way of cashing exactly this one map into photons on a camera.
#
# **What the code does.** `prof = _tf_profile(R[1], R[2])` builds the profile on the 1024² grid;
# `phi_peak(1.5e9, n_col[0])` sets its height.

# %% [cell 67: code]
# ---- STEP 1: the phase map -------------------------------------------------
st_Delta, st_P = 1.5e9, 5.0
st_phi0 = phi_peak(st_Delta, n_col[0])
st_prof = _tf_profile(R[1], R[2])                        # probe along x -> (y,z) plane
st_phi_map = st_phi0 * st_prof

fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.4, 3.6), gridspec_kw=dict(width_ratios=[1.5, 1]))
im = a1.imshow(st_phi_map, extent=ext, origin='lower', cmap='inferno')
a1.set_xlim(-45, 45); a1.set_ylim(-12, 12)
a1.set_xlabel('y (um)'); a1.set_ylabel('z (um)')
a1.set_title('the phase map  $\\varphi(\\mathbf{r})$  the atoms imprint')
plt.colorbar(im, ax=a1, fraction=0.03, label='$\\varphi$ (rad)')
a2.plot(gax*1e6, st_phi_map[Ngrid//2, :], 'C1', lw=2, label='cut along y')
a2.plot(gax*1e6, st_phi_map[:, Ngrid//2], 'C0', lw=2, label='cut along z')
a2.axhline(st_phi0, color='gray', ls=':', lw=1)
a2.annotate(f'$\\varphi_\\mathrm{{peak}}$ = {st_phi0:.3f} rad', (12, st_phi0*0.93), fontsize=9, color='gray')
a2.set_xlim(-45, 45); a2.set_xlabel('position (um)'); a2.set_ylabel('$\\varphi$ (rad)')
a2.set_title('centre cuts'); a2.legend(fontsize=8.5); a2.grid(alpha=0.25)
fig.suptitle('Step 1 — the only thing the atoms give us', fontsize=12)
plt.tight_layout(); plt.show()

# %% [cell 68: markdown]
# #### Step 2 — the field after the atoms: the camera would see *nothing*
#
# **Theory.** The probe leaves the cloud as $E(\mathbf r) = e^{i\varphi(\mathbf r)}$ — a **pure phase
# object**. Its intensity is $|e^{i\varphi}|^2 = 1$ *identically*: the information survives only in
# $\arg E$, which no photodetector responds to.
#
# **What is happening.** Below, the same field twice. Left: its intensity — put a camera here (or at
# any image plane of a perfect system) and the frame is blank. Right: its phase — the condensate,
# perfectly formed, hiding in the one quantity light detectors cannot read. The entire imaging
# problem is getting the right-hand panel to leak into the left-hand one.
#
# **What the code does.** `np.exp(1j*phi_peak_val*prof)` — the first line inside `sim_image`.

# %% [cell 69: code]
# ---- STEP 2: intensity blank, phase structured -------------------------------
st_E_atoms = np.exp(1j * st_phi_map)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.4, 3.4))
im1 = a1.imshow(np.abs(st_E_atoms)**2, extent=ext, origin='lower', cmap='inferno', vmin=0.8, vmax=1.2)
a1.set_xlim(-45, 45); a1.set_ylim(-12, 12); a1.set_xlabel('y (um)'); a1.set_ylabel('z (um)')
a1.set_title('what a camera sees:  $|E|^2$')
plt.colorbar(im1, ax=a1, fraction=0.03, label='$I/I_0$')
a1.annotate('nothing.', (0, 0), color='w', fontsize=13, ha='center')
im2 = a2.imshow(np.angle(st_E_atoms), extent=ext, origin='lower', cmap='inferno')
a2.set_xlim(-45, 45); a2.set_ylim(-12, 12); a2.set_xlabel('y (um)'); a2.set_ylabel('z (um)')
a2.set_title('where the information hides:  $\\arg E$')
plt.colorbar(im2, ax=a2, fraction=0.03, label='phase (rad)')
fig.suptitle('Step 2 — a pure phase object is invisible', fontsize=12)
plt.tight_layout(); plt.show()
print(f"max deviation of |E|^2 from 1 over the whole field: {np.max(np.abs(np.abs(st_E_atoms)**2 - 1)):.1e}"
      f"   -> exactly blank, not approximately blank.")

# %% [cell 70: markdown]
# #### Step 3 — split the field: a big carrier and a small wavelet at 90°
#
# **Theory.** Write $e^{i\varphi} = 1 + w$, with $w \equiv e^{i\varphi}-1$. The "1" is the
# **carrier** — the probe as if no atoms were there. The **wavelet** $w$ carries all the atomic
# information; for small phase, $w \approx i\varphi$: it points at $90°$ to the carrier in the
# complex plane (exactly: $\arg w = 90° + \varphi/2$, $|w| = 2\sin(\varphi/2)$).
#
# **What is happening.** This decomposition *is* the reason for invisibility. Interference changes
# intensity through the cross term $2\,\mathrm{Re}(1^*\!\cdot w)$, and a wavelet at $\sim90°$ has
# almost no real part — the phasor diagram below, drawn **from the actual centre-pixel values of the
# simulated field**, shows the tip-to-tail sum landing right back on the unit circle:
# $|1+w| = 1$. The two arrows are there, carrying 0.2 rad of signal, and they cannot beat against
# each other. Every mode in Steps 6–13 is one way of re-arranging these two arrows.
#
# **What the code does.** the `- 1` inside `fft2(np.exp(1j*phi*prof) - 1)` — the carrier is split
# off *analytically* and handled exactly from here on.

# %% [cell 71: code]
# ---- STEP 3: the wavelet map + the phasor from actual simulated values -------
st_w = st_E_atoms - 1                                  # the wavelet, exactly as the code forms it
_c = Ngrid//2
st_wc = st_w[_c-1:_c+1, _c-1:_c+1].mean()              # actual value at the cloud centre

fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.4, 4.1), gridspec_kw=dict(width_ratios=[1.5, 1]))
im = a1.imshow(np.abs(st_w), extent=ext, origin='lower', cmap='inferno')
a1.set_xlim(-45, 45); a1.set_ylim(-12, 12); a1.set_xlabel('y (um)'); a1.set_ylabel('z (um)')
a1.set_title('the wavelet\'s size:  $|w| = 2\\sin(\\varphi/2)$')
plt.colorbar(im, ax=a1, fraction=0.03, label='$|w|$')

# phasor: carrier, then wavelet tip-to-tail, then their sum -- actual numbers
a2.add_patch(plt.Circle((0, 0), 1.0, fill=False, ls=':', color='gray', lw=1))
a2.annotate('', xy=(1, 0), xytext=(0, 0),
            arrowprops=dict(arrowstyle='-|>', color='#1f5fa8', lw=2.4))
a2.annotate('carrier  1', (0.45, -0.09), color='#1f5fa8', fontsize=10, ha='center')
a2.annotate('', xy=(1+st_wc.real, st_wc.imag), xytext=(1, 0),
            arrowprops=dict(arrowstyle='-|>', color='#c4161c', lw=2.4))
a2.annotate(f'wavelet  $w$\n(|w| = {abs(st_wc):.3f},  at {np.degrees(np.angle(st_wc)):.0f}°)',
            (1.10, 0.10), color='#c4161c', fontsize=9)
a2.annotate('', xy=(1+st_wc.real, st_wc.imag), xytext=(0, 0),
            arrowprops=dict(arrowstyle='-|>', color='k', lw=1.6, ls='--'))
a2.annotate(f'sum:  $|1+w|$ = {abs(1+st_wc):.4f}\nstill ON the unit circle', (0.18, 0.24),
            fontsize=9)
a2.set_xlim(-0.15, 1.45); a2.set_ylim(-0.22, 0.55); a2.set_aspect('equal')
a2.set_xlabel('Re E'); a2.set_ylabel('Im E'); a2.grid(alpha=0.2)
a2.set_title('phasor at the cloud centre (actual values)')
fig.suptitle('Step 3 — carrier + wavelet: 90° apart, so no interference contrast', fontsize=12)
plt.tight_layout(); plt.show()

# %% [cell 72: markdown]
# #### Step 4 — the Fourier plane: the one place the two arrows sit apart
#
# **Theory.** The first lens performs a spatial Fourier transform. The carrier — a plane wave —
# focuses to a delta function at the centre (DC); the wavelet, cloud-shaped in real space, spreads
# into a halo of width $\sim 1/(\text{cloud size})$: the tight $z$ direction (1.5 µm) scatters wide,
# the long $y$ direction (25 µm) scatters narrow. Only what lies inside the NA circle survives.
#
# **What is happening.** In real space, carrier and wavelet overlap everywhere — you cannot touch
# one without the other. *Here* they separate: a physical object of the right size placed at this
# plane (a phase dot, a dark stop) acts on the carrier alone. This plane is the handle for PCI and
# DGI. The right panel shows the aperture doing its (destructive) share: the outer halo — the
# sharpest spatial detail — is gone before any mode-specific optics even begins.
#
# **What the code does.** `fft2(...) * pupil` — and the carrier's delta is never on the grid at all:
# the code carries it as the exact analytic term outside the transform, so the Fourier-plane
# element later acts on it *exactly* (the one idealisation shared by §§7–18: a real dot or stop has
# finite size and would also clip the innermost halo).

# %% [cell 73: code]
# ---- STEP 4: the Fourier plane, before and after the aperture ----------------
st_F = np.fft.fftshift(np.fft.fft2(st_w))
st_P2 = np.abs(st_F)**2 / np.abs(st_F).max()**2
st_fax = np.fft.fftshift(np.fft.fftfreq(Ngrid, FOV/Ngrid)) * 1e-6
st_fNA = NA/lam * 1e-6
st_pup = np.fft.fftshift(pupil)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.4, 4.6))
for a, dat, ttl in [(a1, st_P2, 'before the aperture'),
                    (a2, st_P2*st_pup, 'after the aperture:  outer halo gone')]:
    im = a.imshow(np.log10(dat + 1e-12), extent=[st_fax[0], st_fax[-1], st_fax[0], st_fax[-1]],
                  origin='lower', cmap='magma', vmin=-8, vmax=0)
    th_c = np.linspace(0, 2*np.pi, 200)
    a.plot(st_fNA*np.cos(th_c), st_fNA*np.sin(th_c), 'c-', lw=1.5)
    a.set_xlim(-0.55, 0.55); a.set_ylim(-0.55, 0.55)
    a.set_xlabel('$f_y$ (1/um)'); a.set_ylabel('$f_z$ (1/um)'); a.set_title(ttl)
a1.plot(0, 0, 'w+', ms=13, mew=2.5)
a1.annotate('carrier: analytic delta at DC\n(the big arrow of Step 3, parked here)',
            (0.04, 0.33), color='w', fontsize=8.5)
a1.annotate('the wavelet\'s halo', (-0.5, -0.45), color='w', fontsize=9)
a1.annotate('NA circle', (st_fNA*0.6, st_fNA*1.06), color='c', fontsize=9)
plt.colorbar(im, ax=a2, fraction=0.045, label='$\\log_{10}$ power (norm.)')
fig.suptitle('Step 4 — carrier and wavelet finally separate in space', fontsize=12)
plt.tight_layout(); plt.show()
print(f"scattered power inside the NA circle: "
      f"{(st_P2*st_pup).sum()/st_P2.sum():.3f}  -> the aperture keeps ~71% of the wavelet;"
      f" this loss is the blur factor of Sec. 7 and is paid by every mode alike.")

# %% [cell 74: markdown]
# #### Step 5 — through a *perfect* imaging system: still invisible
#
# **Theory.** With nothing placed at the Fourier plane, the second lens reassembles
# $E = 1 + w_\mathrm{NA}$ (the wavelet now band-limited). Since the surviving wavelet is still
# essentially at 90° to the carrier, the cross term is still $\approx 0$: the intensity is flat up
# to second-order crumbs, $|1+w_\mathrm{NA}|^2 \approx 1 + \mathcal{O}(\varphi^2) + \text{clipping ripples}$.
#
# **What is happening.** This is the dead end the story has been building to: a *perfect,
# aberration-free, lossless* 4f system faithfully delivers the phase object to its image plane —
# still as a phase object. The colour bar below spans ±6% and the cloud barely registers. Nothing in
# the optics so far has converted signal into intensity, because nothing has touched the *angle
# between the two arrows*. That, and only that, is what the four modes now do — two of them by
# grabbing the carrier at the Fourier plane (Steps 6–9), two by splitting the field into
# polarizations that rotate against each other (Steps 10–13).
#
# **What the code does.** `sim_image(axis, phi, mode)` with the fall-through branch `E = 1 + Esc` —
# the "clear" reference the noisy frames of §7 are normalised against.

# %% [cell 75: code]
# ---- STEP 5: the image plane with no trick installed --------------------------
st_I_clear, _ = sim_image(0, st_phi0, 'clear')

fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.4, 3.4), gridspec_kw=dict(width_ratios=[1.5, 1]))
im = a1.imshow(st_I_clear, extent=ext, origin='lower', cmap='inferno', vmin=0.94, vmax=1.06)
a1.set_xlim(-45, 45); a1.set_ylim(-12, 12); a1.set_xlabel('y (um)'); a1.set_ylabel('z (um)')
a1.set_title('image-plane intensity, no Fourier-plane element')
plt.colorbar(im, ax=a1, fraction=0.03, label='$I/I_0$')
a2.plot(gax*1e6, st_I_clear[:, Ngrid//2], 'C0', lw=1.8, label='cut along z')
a2.axhline(1, color='gray', ls=':', lw=1)
a2.set_xlim(-12, 12); a2.set_ylim(0.98, 1.02); a2.set_xlabel('z (um)'); a2.set_ylabel('$I/I_0$')
a2.set_title('centre cut: crumbs of order $\\varphi^2$'); a2.legend(fontsize=8.5); a2.grid(alpha=0.25)
fig.suptitle('Step 5 — the dead end: a perfect lens system still shows (almost) nothing', fontsize=12)
plt.tight_layout(); plt.show()
print(f"peak deviation from flatness: {np.max(np.abs(st_I_clear-1)):.3f} of I0 -- second-order only, "
      f"while the phase signal sitting in the field is {st_phi0:.3f} rad. Time for the tricks.")

# %% [cell 76: markdown]
# ### 19.2 PCI — rotate the carrier (Steps 6–7)
#
# #### Step 6 — the phase plate turns the carrier by 90°, and the cloud appears
#
# **Theory.** A small dot at the Fourier-plane centre retards (only) the carrier by $\pi/2$ and
# transmits $t_p = 0.95$ of its amplitude: $1 \mapsto i\,t_p$. Nothing happens to the wavelet. Now
# the two arrows are (nearly) **parallel**, and the cross term switches on:
# $$I = |i\,t_p + w|^2 = t_p^2 + 2t_p\sin\varphi + 2(1-\cos\varphi) \approx t_p^2 + 2\,t_p\,\varphi .$$
#
# **What is happening.** Compare the two phasor panels below — same actual centre-pixel values as
# Step 3, before and after the plate. The wavelet has not changed at all; the carrier has been
# *rotated under it*. Tip-to-tail, the sum now pokes **outside** the circle of radius $t_p$: the
# cloud becomes a bright object on the $t_p^2$ background, linear in $\varphi$. The third panel is
# the image plane: the same field as Step 5, transformed from blank to a condensate portrait by one
# rotation.
#
# **What the code does.** `E = t_p*np.exp(1j*theta) + Esc` — the entire optical trick is that one
# analytic carrier coefficient.

# %% [cell 77: code]
# ---- STEP 6: phasor before / after the plate + the image that appears ---------
st_I_pci, _ = sim_image(0, st_phi0, 'PCI')
st_carrier_pci = t_p*np.exp(1j*theta)

fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(14.6, 4.3), gridspec_kw=dict(width_ratios=[1, 1, 1.5]))
for a, car, lab, ttl in [(a1, 1+0j, 'carrier  1', 'before the plate:  $|1+w|^2 = 1$'),
                         (a2, st_carrier_pci, 'carrier  $i\\,t_p$', 'after:  $|i t_p + w|^2$ = %.2f' %
                          abs(st_carrier_pci+st_wc)**2)]:
    a.add_patch(plt.Circle((0, 0), 1.0, fill=False, ls=':', color='gray', lw=1))
    a.add_patch(plt.Circle((0, 0), t_p, fill=False, ls=':', color='#b8b8b8', lw=1))
    a.annotate('', xy=(car.real, car.imag), xytext=(0, 0),
               arrowprops=dict(arrowstyle='-|>', color='#1f5fa8', lw=2.4))
    a.annotate(lab, (car.real*0.5+0.06, car.imag*0.5-0.02), color='#1f5fa8', fontsize=9)
    tip = car + st_wc
    a.annotate('', xy=(tip.real, tip.imag), xytext=(car.real, car.imag),
               arrowprops=dict(arrowstyle='-|>', color='#c4161c', lw=2.4))
    a.annotate('$w$ (unchanged)', (tip.real+0.05, tip.imag-0.03), color='#c4161c', fontsize=9)
    a.annotate('', xy=(tip.real, tip.imag), xytext=(0, 0),
               arrowprops=dict(arrowstyle='-|>', color='k', lw=1.6, ls='--'))
    a.set_xlim(-0.45, 1.35); a.set_ylim(-0.15, 1.35); a.set_aspect('equal')
    a.set_xlabel('Re E'); a.set_ylabel('Im E'); a.grid(alpha=0.2); a.set_title(ttl, fontsize=10)
a2.annotate('arrows now parallel:\ninterference ON', (0.30, 0.30), fontsize=9)
im = a3.imshow(st_I_pci, extent=ext, origin='lower', cmap='inferno', vmin=0.85, vmax=1.20)
a3.set_xlim(-45, 45); a3.set_ylim(-12, 12); a3.set_xlabel('y (um)'); a3.set_ylabel('z (um)')
a3.set_title('the image plane, one rotation later')
plt.colorbar(im, ax=a3, fraction=0.03, label='$I/I_0$')
fig.suptitle('Step 6 — PCI: rotate the carrier under the wavelet', fontsize=12)
plt.tight_layout(); plt.show()
print(f"centre of the simulated image: I/I0 = {st_I_pci[Ngrid//2-1:Ngrid//2+1, Ngrid//2-1:Ngrid//2+1].mean():.3f}"
      f"   (infinite-NA phasor value {abs(st_carrier_pci+st_wc)**2:.3f}; the gap is Step 4's aperture loss)")

# %% [cell 78: markdown]
# #### Step 7 — onto the camera: pixels, photons, noise
#
# **Theory.** The camera integrates the image over 15×15-sample bins (one 2.93 µm pixel each), each
# bin collecting $N_\mathrm{ph} \approx 3800$ detected photons at this power — then nature rolls the
# dice: Poisson shot noise on every pixel plus 7 e⁻ of Gaussian read noise (§7.4).
#
# **What is happening.** The left panel is the same field as Step 6 after binning — what the
# experiment would record with a noiseless detector. The right panel is *one actual noisy frame*:
# the condensate is comfortably visible in a single shot, riding on the granular $t_p^2$ background —
# the shot noise of the bright carrier is the price PCI pays for linearity. The lineout shows the
# noisy pixels scattered around the noiseless profile, which is precisely the SNR ≈ 15/pixel of §8.
#
# **What the code does.** `to_camera(I_img, P_mW)` — binning, `rng.poisson`, `rng.normal`, in three
# lines.

# %% [cell 79: code]
# ---- STEP 7: PCI, the recorded frame ------------------------------------------
st_cam_pci, st_ideal_pci = to_camera(st_I_pci, st_P)
st_nb = Ngrid//15; st_ycam = (np.arange(st_nb)-st_nb//2+0.5)*pix_obj*1e6; st_mid = st_nb//2

fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(14.6, 3.5))
im1 = a1.imshow(st_ideal_pci, extent=ext, origin='lower', cmap='inferno', vmin=0.85, vmax=1.20)
a1.set_title('binned, noiseless'); plt.colorbar(im1, ax=a1, fraction=0.03)
im2 = a2.imshow(st_cam_pci, extent=ext, origin='lower', cmap='inferno', vmin=0.85, vmax=1.20)
a2.set_title('ONE noisy frame (what you actually get)'); plt.colorbar(im2, ax=a2, fraction=0.03, label='$I/I_0$')
for a in (a1, a2):
    a.set_xlim(-45, 45); a.set_ylim(-12, 12); a.set_xlabel('y (um)'); a.set_ylabel('z (um)')
a3.plot(st_ycam, st_ideal_pci[st_mid], 'C0', lw=2, label='noiseless')
a3.plot(st_ycam, st_cam_pci[st_mid], 'C0.', ms=4, alpha=0.75, label='noisy pixels')
a3.axhline(bg_plate, color='gray', ls=':', lw=1, label='$t_p^2$ background')
a3.set_xlim(-45, 45); a3.set_xlabel('y (um)'); a3.set_ylabel('$I/I_0$')
a3.set_title('lineout along the cigar'); a3.legend(fontsize=8); a3.grid(alpha=0.25)
fig.suptitle('Step 7 — PCI recorded: linear signal on a bright, shot-noisy background', fontsize=12)
plt.tight_layout(); plt.show()

# %% [cell 80: markdown]
# ### 19.3 DGI — delete the carrier (Steps 8–9)
#
# #### Step 8 — the dark stop removes the big arrow; only the wavelet flies on
#
# **Theory.** Instead of rotating the carrier, absorb it: an OD-4 stop leaves
# $1 \mapsto \varepsilon = 10^{-\mathrm{OD}/2} = 0.01$. The image field is essentially the wavelet
# alone:
# $$I = |\varepsilon + w|^2 \approx |w|^2 = 4\sin^2(\varphi/2) \approx \varphi^2 ,$$
# **quadratic** in the signal, on a background of the leak, $\varepsilon^2 = 10^{-4}$.
#
# **What is happening.** In the phasor picture the blue arrow is simply *gone* — what remains is the
# red wavelet, whose squared length is now the intensity directly. No interference, no local
# oscillator: the cloud glows by its own scattered light on a black field, the classic dark-ground
# image of Andrews et al. (1996). The catch is in the numbers: $|w|^2 = 0.041$ — the signal is
# *self-referenced* and therefore small-squared. Note the colour bar relative to Step 6.
#
# **What the code does.** `E = 10**(-OD/2) + Esc` — again a single analytic carrier coefficient;
# same wavelet, same propagation, same everything else.

# %% [cell 81: code]
# ---- STEP 8: carrier deleted, the wavelet's own light ------------------------
st_I_dgi, _ = sim_image(0, st_phi0, 'DGI')
st_eps = 10**(-4.0/2)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.4, 4.1), gridspec_kw=dict(width_ratios=[1, 1.5]))
a1.add_patch(plt.Circle((0, 0), 1.0, fill=False, ls=':', color='gray', lw=1))
a1.annotate('', xy=(st_eps, 0), xytext=(0, 0),
            arrowprops=dict(arrowstyle='-|>', color='#1f5fa8', lw=2.4))
a1.annotate('carrier: $\\varepsilon$ = 0.01\n(all but deleted)', (0.05, -0.13), color='#1f5fa8', fontsize=9)
tip = st_eps + st_wc
a1.annotate('', xy=(tip.real, tip.imag), xytext=(st_eps, 0),
            arrowprops=dict(arrowstyle='-|>', color='#c4161c', lw=2.4))
a1.annotate('$w$: now the whole field\n$|w|^2$ = %.3f' % abs(st_wc)**2, (0.09, 0.14), color='#c4161c', fontsize=9)
a1.set_xlim(-0.30, 1.10); a1.set_ylim(-0.30, 0.60); a1.set_aspect('equal')
a1.set_xlabel('Re E'); a1.set_ylabel('Im E'); a1.grid(alpha=0.2)
a1.set_title('the phasor after the stop')
im = a2.imshow(st_I_dgi, extent=ext, origin='lower', cmap='inferno', vmin=0, vmax=0.035)
a2.set_xlim(-45, 45); a2.set_ylim(-12, 12); a2.set_xlabel('y (um)'); a2.set_ylabel('z (um)')
a2.set_title('the image plane: self-luminous cloud on black')
plt.colorbar(im, ax=a2, fraction=0.03, label='$I/I_0$')
fig.suptitle('Step 8 — DGI: delete the carrier, keep the wavelet', fontsize=12)
plt.tight_layout(); plt.show()
print(f"peak of the dark-ground image: {st_I_dgi.max():.4f} of I0 (vs 4 sin^2(phi/2) = "
      f"{4*np.sin(st_phi0/2)**2:.4f} before aperture loss); background floor eps^2 = {st_eps**2:.0e}.")

# %% [cell 82: markdown]
# #### Step 9 — the dark frame on the camera
#
# **Theory.** The same detection step as Step 7, but the arithmetic is now brutally different: the
# peak intensity is 0.03 of $I_0$, so the brightest pixel collects ~100 detected photons instead of
# ~4800 — shot noise of ~10 e⁻ against a 7 e⁻ read floor.
#
# **What is happening.** One noisy dark-ground frame. The cloud is still detectable in a single shot
# (SNR ≈ 6/pixel, §18.4) but visibly grainier than PCI's — the quadratic transfer starved it of
# photons *before* the noise was ever rolled. In exchange: essentially zero background, so nothing
# of the probe's brightness or flicker contaminates the frame. This trade — signal linearity versus
# background darkness — is the entire PCI-vs-DGI decision of §9, now visible as two camera frames
# of the same cloud.
#
# **What the code does.** identical `to_camera(...)` call; only the input image changed.

# %% [cell 83: code]
# ---- STEP 9: DGI, the recorded frame ------------------------------------------
st_cam_dgi, st_ideal_dgi = to_camera(st_I_dgi, st_P)

fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(14.6, 3.5))
im1 = a1.imshow(st_ideal_dgi, extent=ext, origin='lower', cmap='inferno', vmin=-0.005, vmax=0.03)
a1.set_title('binned, noiseless'); plt.colorbar(im1, ax=a1, fraction=0.03)
im2 = a2.imshow(st_cam_dgi, extent=ext, origin='lower', cmap='inferno', vmin=-0.005, vmax=0.03)
a2.set_title('ONE noisy frame'); plt.colorbar(im2, ax=a2, fraction=0.03, label='$I/I_0$')
for a in (a1, a2):
    a.set_xlim(-45, 45); a.set_ylim(-12, 12); a.set_xlabel('y (um)'); a.set_ylabel('z (um)')
a3.plot(st_ycam, st_ideal_dgi[st_mid], 'C3', lw=2, label='noiseless')
a3.plot(st_ycam, st_cam_dgi[st_mid], 'C3.', ms=4, alpha=0.75, label='noisy pixels')
a3.axhline(0, color='gray', ls=':', lw=1)
a3.set_xlim(-45, 45); a3.set_xlabel('y (um)'); a3.set_ylabel('$I/I_0$')
a3.set_title('lineout: small signal, dark floor'); a3.legend(fontsize=8); a3.grid(alpha=0.25)
fig.suptitle('Step 9 — DGI recorded: fewer photons, but nothing behind them', fontsize=12)
plt.tight_layout(); plt.show()

# %% [cell 94: markdown]
# ---
# **The story in one paragraph.** The atoms hand the light a phase map and nothing else (Step 1);
# phase is invisible to any camera (Step 2) because the information rides in a small wavelet at 90°
# to a big carrier (Step 3); the two separate only at the Fourier plane (Step 4), and a perfect lens
# system alone reunites them as invisibly as they arrived (Step 5). PCI rotates the carrier under
# the wavelet and buys a *linear*, bright-background image (Steps 6–7); DGI deletes the carrier and
# buys a *quadratic*, dark-background one (Steps 8–9); Faraday sidesteps the Fourier plane entirely
# by letting the two circular halves of the beam phase against each other — read through a crossed
# polarizer it is DGI's quadratic twin (Steps 10–11), read through two balanced ports it is PCI's
# linear twin with built-in immunity to the probe's own noise (Steps 12–13). All four spend the
# identical destruction budget; the filmstrip of Step 14 is that budget being spent. The tables of
# §17.3 and §18.4 are the accounting; this section is the film.

# %% [cell 95: markdown]
# ---
# ## References
#
# The apparatus and BEC operating point used throughout are taken from the group's own thesis (K24);
# the imaging techniques and destruction-budget arguments draw on the dispersive- and
# Faraday-imaging literature below. (This consolidates and extends the short citation list in the
# introduction at the top of the notebook.)
#
# **Apparatus**
#
# 1. Kucera, DPhil thesis, University of Oxford (2024) [K24] -- imaging arm (§3.1.3), 401 nm
#    transition (§3.2.2), BEC operating point (§6.3.1).
#
# **Dispersive imaging of BECs (PCI / DGI)**
#
# 2. W. Ketterle, D. S. Durfee, D. M. Stamper-Kurn, "Making, probing and understanding
#    Bose-Einstein condensates," in *Proceedings of the International School of Physics "Enrico
#    Fermi," Course CXL* (IOS Press, 1999), arXiv:cond-mat/9904034 -- dispersive-imaging and
#    phase-contrast review used throughout §5-§8.
# 3. M. R. Andrews, M.-O. Mewes, N. J. van Druten, D. S. Durfee, D. M. Kurn, W. Ketterle, "Direct,
#    Nondestructive Observation of a Bose Condensate," Science **273**, 84 (1996) -- first
#    non-destructive dark-ground observation of a BEC; motivates the DGI mode of §7.1.
# 4. R. Meppelink, R. A. Rozendaal, S. B. Koller, J. M. Vogels, P. van der Straten,
#    "Thermodynamics of Bose-Einstein-condensed clouds using phase-contrast imaging," Phys. Rev. A
#    **81**, 053632 (2010) -- primary reference for the PCI transfer-curve convention used in
#    §5.3 and §7.2.
# 5. F. Böttcher, J.-N. Schmidt, M. Wenzel, J. Hertkorn, M. Guo, T. Langen, T. Pfau, "Transient
#    Supersolid Properties in an Array of Dipolar Quantum Droplets," Phys. Rev. X **9**, 011051
#    (2019) -- single-shot, in-situ far-detuned phase-contrast imaging of a strongly dipolar (Dy)
#    condensate at the analogous 421 nm line; cited in §12 as the "right-substance" precedent
#    for the revised destruction model.
#
# **Faraday (polarization-rotation) imaging -- new in §17**
#
# 6. M. Gajdacz, P. L. Pedersen, T. Mørch, A. J. Hilliard, J. Arlt, J. F. Sherson,
#    "Non-destructive Faraday imaging of dynamically controlled ultracold atoms," Rev. Sci.
#    Instrum. **84**, 083105 (2013) -- the dark-field (crossed-polarizer) configuration of
#    §17.2(a); reports up to 2000 images of the same cloud.
# 7. F. Kaminski, N. S. Kampel, M. P. H. Steenstrup, A. Griesmaier, E. S. Polzik, J. H. Müller,
#    "In-situ dual-port polarization contrast imaging of Faraday rotation in a high optical depth
#    ultracold $^{87}$Rb atomic ensemble," Eur. Phys. J. D **66**, 227 (2012) -- the balanced
#    dual-port configuration of §17.2(b).
# 8. M. A. Kristensen, M. Gajdacz, P. L. Pedersen, C. Klempt, J. F. Sherson, J. J. Arlt, A. J.
#    Hilliard, "Sub-atom shot noise Faraday imaging of ultracold atom clouds," J. Phys. B **50**,
#    034004 (2017) -- quantitative shot-noise-limited precision analysis behind the SNR comparison
#    of §17.3.
#
# **Erbium apparatus context**
#
# 9. K. Aikawa, A. Frisch, M. Mark, S. Baier, A. Rietzler, R. Grimm, F. Ferlaino, "Bose-Einstein
#    Condensation of Erbium," Phys. Rev. Lett. **108**, 210401 (2012) -- first quantum-degenerate
#    erbium gas ($^{168}$Er); established the 401 nm cooling/imaging transition shared with the
#    $^{166}$Er apparatus modelled here.
#
# ---
# *Note on $\kappa_F$ (§17.1): the vector-coupling coefficient for the 401 nm $^{166}$Er
# stretched-state transition was not found in the literature searched for this notebook and is left
# as an explicit placeholder ($\kappa_F=1$, an idealised bound) pending either an atomic-structure
# calculation or a lab measurement -- the same status as the `QE_cam`/`read_e` placeholders of
# §2.*
