# 08 Shot Noise
#
# Exported from: 1 calculations revised 2  multishot  6  extended.ipynb
# This file is a mechanical notebook-section export for refactoring.
# Keep physics, numerical algorithms, parameter values, and execution order equivalent to the notebook.


# %% [cell 21: markdown]
# ## 8. Photon shot-noise SNR
#
# ### 8.1 Idealised single-shot SNR (original)
# $\mathrm{SNR}_\mathrm{shot}=2t_p\varphi\sqrt{N_\mathrm{phot}}$ at the centre pixel — the
# photon-shot-noise limit at $\mathrm{QE}=1$, infinite resolution, linear regime.
#
# ### 8.2 Realistic detection [rev]
# Four corrections, all now applied consistently:
# * **Peak intensity.** Detected photons use the on-axis probe intensity (§6.1), the
#   same convention as scattering — the cloud and the pixels imaging it both sit at the
#   beam centre. (An earlier draft used the area average here, a factor-$\sqrt2$ pessimism, now removed.)
# * **Plate-dimmed background.** The processed ratio references the no-atom, plate-in
#   image at level $t_p^2$, so the signal is $I_\mathrm{img}-t_p^2$ and the local shot
#   noise scales with the actual level $I_\mathrm{img}$, not $1$.
# * **Full-model sensitivity & a phase-regime guard.** Contrast uses the full transfer
#   $I=|t_p e^{i\theta}+e^{i\varphi}-1|^2$, not the linear tangent, and SNR is **only
#   reported where $\varphi<\pi$** (single-valued); beyond that the inversion is
#   ambiguous and a linear SNR is meaningless (this is why the along-cigar axis at
#   1.5 GHz, $\varphi_y=4.25$ rad, returns NaN rather than a spurious large number).
# * **Axis-aware NA blur** $b_\mathrm{axis}$ from §7.3.
#
# ### 8.3 Per-pixel vs per-resolution-element, from the simulation
# Because the cloud is unresolved across the cigar, the honest per-resolution-element
# SNR is obtained by **integrating the real Fourier-optics image** over a resolution
# element, not by an analytic $b\times\mathrm{bin}^2$ formula (which would wrongly treat
# NA-lost light as recoverable). Both are reported.
#
# ### 8.4 Accumulated-SNR invariance
# $\mathrm{SNR}_\mathrm{total}=\mathrm{SNR}_\mathrm{shot}\sqrt{N}$ with $N\le N_\mathrm{max}$.
# Because $\mathrm{SNR}_\mathrm{shot}\propto\sqrt{P\tau}/\delta$ and $N_\mathrm{max}\propto\delta^2/(P\tau)$,
# the **photon-limited** total SNR is independent of $(\Delta,P,\tau)$ — set only by the
# destruction budget. Read noise breaks this invariance (it favours fewer, brighter shots).

# %% [cell 22: code]
def I_full(phi):
    return np.abs(t_p*np.exp(1j*theta) + np.exp(1j*phi) - 1)**2   # incident-I0 units; no-atom = t_p^2

def regime_of(phi):
    return 'linear' if phi < 0.5 else ('Meppelink' if phi < np.pi else 'phase-wrapped')

def SNR_shot_ideal(Delta_Hz, P_mW, axis=0, tau_s=None):
    # original idealised reference: linear contrast, QE=1, infinite NA
    phi = phi_peak(Delta_Hz, n_col[axis])
    Nph = N_phot_pix(P_mW, tau_s, QE=QE_ideal)
    return 2*t_p*phi*np.sqrt(Nph)

def SNR_pixel(Delta_Hz, P_mW, axis=0, tau_s=None):
    # realistic single-pixel SNR at cloud centre (full model, peak intensity, plate bg, axis blur)
    phi = phi_peak(Delta_Hz, n_col[axis])
    if phi >= np.pi:
        return np.nan                       # ambiguous; not a single-valued measurement
    Iimg = bg_plate + blur_axis[axis]*(I_full(phi) - bg_plate)
    Nph  = N_phot_pix(P_mW, tau_s)
    return abs(Iimg - bg_plate)*Nph/np.sqrt(abs(Iimg)*Nph + read_e**2)

def SNR_reselem_sim(Delta_Hz, P_mW, axis=0, tau_s=None, half=1):
    # trustworthy per-resolution-element SNR by integrating the real simulated image
    phi = phi_peak(Delta_Hz, n_col[axis])
    if phi >= np.pi:
        return np.nan
    Iimg, _ = sim_image(axis, phi, 'PCI')
    nb = (Ngrid//15)*15
    binned = Iimg[:nb,:nb].reshape(nb//15,15,nb//15,15).mean(axis=(1,3))
    Nd = N_phot_pix(P_mW, tau_s)
    m = binned.shape[0]//2
    blk = binned[m-half:m+half+1, m-half:m+half+1]
    S = (blk - bg_plate).sum()*Nd
    noise = np.sqrt(blk.sum()*Nd + blk.size*read_e**2)
    return S/noise

print("Across-cigar (x), reference detuning 1.5 GHz, regime:", regime_of(phi_peak(1.5e9, n_col[0])))
for P in [2.0, 5.0, 10.0]:
    print(f"  P={P:4.1f} mW:  ideal(QE=1)={SNR_shot_ideal(1.5e9,P):5.1f}   "
          f"realistic/pixel={SNR_pixel(1.5e9,P):4.1f}   "
          f"realistic/res-elem(sim)={SNR_reselem_sim(1.5e9,P):4.1f}")
dgi_peak = sim_image(0, phi_peak(1.5e9, n_col[0]), 'DGI')[0].max()
_dgi5 = dgi_peak*N_phot_pix(5.0)
_dgi_stat = 'quantitative (> 3x read noise)' if _dgi5 > 3*read_e else ('above the read-noise floor' if _dgi5 > read_e else 'below the read-noise floor')
print(f"\nDGI single-shot of the BEC (5 mW): peak signal ~{dgi_peak:.3f} of I0 = "
      f"{_dgi5:.1f} e- vs {read_e} e- read -> {_dgi_stat}")

# --- numerical demonstration of the accumulated-SNR invariance (photon-limited) ---
print("\n[8.4 check] photon-limited total SNR = SNR_shot(QE=1) x sqrt(Nmax) over a full budget:")
print("   detuning   SNR_shot   Nmax(30%)   accumulated = SNR_shot*sqrt(Nmax)")
for Dg in [0.75, 1.0, 1.5, 2.0, 3.0]:
    ss = SNR_shot_ideal(Dg*1e9, 2.0); nm = Nmax_loss(Dg*1e9, 2.0, 0.30)
    print(f"   {Dg:4.2f} GHz   {ss:6.1f}     {nm:7.0f}      {ss*np.sqrt(nm):8.1f}")
print("   -> ~constant: total photon-limited SNR is set by the destruction budget, not by (Delta,P,tau)")

# %% [cell 23: markdown]
# ## 9. Maximum non-destructive shots — threshold family [revised]
#
# The shot budget follows from the loss model $N_0(s)=N_0e^{-\eta N_\gamma s}$ (§6.3).
# For a tolerated condensate-loss fraction $f$, $N_\mathrm{max}=-\ln(1-f)/(\eta N_\gamma)$.
# Several thresholds are physically motivated; we evaluate all of them and let the plots
# decide which is most reasonable.
#
# | threshold | meaning | equivalent loss |
# |---|---|---|
# | $\varepsilon=\mu$ (original) | heating $=\mu$; with $E_\mathrm{rec}\gg\mu$ this is $\mu/E_\mathrm{rec}=13\%$ loss | 13% |
# | 30% loss | the protocols-doc §4.5 acceptance test ("<30% loss over ~100 shots") | 30% |
# | 50% loss | half-life of the condensate | 50% |
# | $\varepsilon=k_BT_\mathrm{cloud}$ | thermal-pedestal ceiling (different observable) | — |
# | 1 photon/atom | classic dispersive-imaging bar (Ketterle/Andrews) | — |
#
# Collisional secondaries enter through $\eta$; the original used $\eta=1$.

# %% [cell 24: code]
# Nmax_loss and Nmax_heat are defined in section 6.3 (the loss model); applied here.
D, P = 1.5e9, 2.0
print("At reference point (1.5 GHz, 2 mW), peak-intensity convention, eta = %.1f:" % eta_coll)
print(f"  eps=mu (=13% loss)     N_max = {Nmax_loss(D,P,T_mu/T_rec):5.0f}")
print(f"  30% loss               N_max = {Nmax_loss(D,P,0.30):5.0f}")
print(f"  50% loss               N_max = {Nmax_loss(D,P,0.50):5.0f}")
print(f"  eps=kT_cloud           N_max = {Nmax_heat(D,P,T_cloud):5.0f}")
print(f"  1 photon/atom (eta=1)  N_max = {1.0/N_scatt(D,P):5.0f}")
print(f"\nFor comparison, the ORIGINAL notebook (area-avg intensity, eta=1, eps=mu): N_max = {(T_mu/T_rec)/(N_scatt(D,P)/2):.0f}")
print(f"  -> the two corrections (peak intensity x2, eta={eta_coll}) reduce eps=mu N_max from ~134 to ~{Nmax_loss(D,P,T_mu/T_rec):.0f}")
print(f"  -> but the recommended 30%-loss criterion gives N_max = {Nmax_loss(D,P,0.30):.0f}, still clearing the 100-shot goal")

# %% [cell 25: markdown]
# ### Figure — shot budget vs detuning for the candidate thresholds

# %% [cell 26: code]
Dg = np.linspace(0.5, 4.0, 400)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 5.0),
                             gridspec_kw=dict(width_ratios=[1.5,1], wspace=0.3))
curves = [
    (Nmax_loss(Dg*1e9, 2.0, T_mu/T_rec), r'$\varepsilon=\mu$ (=13% loss)', 'tab:blue'),
    (Nmax_loss(Dg*1e9, 2.0, 0.30),       '30% loss (protocols 4.5)',      'tab:red'),
    (Nmax_loss(Dg*1e9, 2.0, 0.50),       '50% loss',                       'tab:orange'),
    (Nmax_heat(Dg*1e9, 2.0, T_cloud),    r'$k_B T_{cloud}$',               'tab:green'),
    (1.0/N_scatt(Dg*1e9, 2.0),           '1 photon/atom',                  'tab:gray'),
]
for y, lab, col in curves:
    a1.semilogy(Dg, y, color=col, lw=2.2, label=lab)
a1.fill_between(Dg, Nmax_loss(Dg*1e9,2.0,0.30,eta=2.0), Nmax_loss(Dg*1e9,2.0,0.30,eta=1.0),
                color='tab:red', alpha=0.15, label=r'30% curve, $\eta\in[1,2]$')
a1.axhline(100, color='k', ls='--', lw=0.8); a1.axvline(1.5, color='k', ls=':', lw=0.8)
a1.annotate('Andrews 1996 ~100 shots (Na, larger mu)', (0.55, 110), fontsize=8)
a1.set_xlabel('Detuning (GHz)'); a1.set_ylabel('N_max  (P=2 mW, tau=40 us)')
a1.set_ylim(3, 5e3); a1.grid(alpha=0.3, which='both'); a1.legend(fontsize=8, loc='upper left')
a1.set_title('(a) Shot budget vs detuning, by threshold')

labels = ['eps=mu\n(13% loss)','30%\nloss','50%\nloss','kT_cloud','1 phot\n/atom']
vals = [Nmax_loss(1.5e9,2.0,T_mu/T_rec), Nmax_loss(1.5e9,2.0,0.30), Nmax_loss(1.5e9,2.0,0.50),
        Nmax_heat(1.5e9,2.0,T_cloud), 1.0/N_scatt(1.5e9,2.0)]
cols = ['tab:blue','tab:red','tab:orange','tab:green','tab:gray']
bars = a2.bar(labels, vals, color=cols, alpha=0.85)
for b,v in zip(bars, vals): a2.text(b.get_x()+b.get_width()/2, v*1.05, f'{v:.0f}', ha='center', fontsize=9)
a2.set_yscale('log'); a2.set_ylim(5, 2e3); a2.set_ylabel('N_max'); a2.grid(alpha=0.3, axis='y', which='both')
a2.set_title('(b) At (1.5 GHz, 2 mW)')
fig.suptitle('Destruction budget: E_rec/kB=359 nK >> mu/kB=48 nK  =>  loss model, not heating', y=1.02)
plt.tight_layout(); plt.show()

# %% [cell 33: markdown]
# ### 12.1 Pulse-duration sweep at the reference detuning
#
# At $\Delta=1.5$ GHz: as $\tau$ grows, PCI SNR rises, the DGI signal crosses the read-noise floor
# (quantitative for $\tau\gtrsim40$–$50\,\mu$s at 3.5–5 mW), and $N_\mathrm{max}$ falls along a curve of
# constant $N_\mathrm{max}\tau$ (right axis of panel b).

# %% [cell 34: code]
# ===== CELL F: tau sweep — SNR rises, DGI clears the read-noise floor, N_max falls at fixed N_max*tau =====
tau_us = np.array([20,30,40,50,65,80,100,120,150,200,250,300])
Dref, Pref, Pdgi = 1.5e9, 2.0, 5.0
snr_res = np.array([SNR_reselem_sim(Dref,Pref,0,t*1e-6) for t in tau_us])
snr_pix = np.array([SNR_pixel(Dref,Pref,0,t*1e-6)       for t in tau_us])
dgi_e   = np.array([0.016*N_phot_pix(Pdgi,t*1e-6)        for t in tau_us])
nm_real = np.array([Nmax_heating(Dref,Pref,0.30,t*1e-6)  for t in tau_us])
nm_opt  = np.array([Nmax_cleanloss(Dref,Pref,0.30,t*1e-6) for t in tau_us])

fig,(axA,axB)=plt.subplots(1,2,figsize=(13.2,4.9))
tcap=tau_max_now*1e6

# (a) per-shot SNR (left axis) and DGI signal (right axis)
l1,=axA.plot(tau_us,snr_res,'o-',color='#1f5fa8',lw=2,label='PCI SNR / resolution element')
l2,=axA.plot(tau_us,snr_pix,'s--',color='#1f5fa8',lw=1.5,alpha=0.55,label='PCI SNR / pixel')
axA.set_xlabel('pulse duration  τ (µs)'); axA.set_ylabel('PCI SNR',color='#1f5fa8')
axA.tick_params(axis='y',labelcolor='#1f5fa8'); axA.set_ylim(0,max(snr_res)*1.15)
axR=axA.twinx()
l3,=axR.plot(tau_us,dgi_e,'^-',color='#c4161c',lw=2,label='DGI peak signal (5 mW)')
axR.axhline(read_e,color='#c4161c',ls=':',lw=1.3); axR.axhline(3*read_e,color='#c4161c',ls=':',lw=1,alpha=0.5)
axR.set_ylabel('DGI peak signal  (e⁻)',color='#c4161c'); axR.tick_params(axis='y',labelcolor='#c4161c')
axR.set_ylim(0,max(dgi_e)*1.1)
axR.text(tau_us[-1],read_e+3,'read-noise floor (7 e⁻)',ha='right',va='bottom',color='#c4161c',fontsize=8.5)
axR.text(tau_us[-1],3*read_e+3,'3× read → DGI quantitative',ha='right',va='bottom',color='#c4161c',fontsize=8,alpha=0.8)
axA.axvspan(tcap,tau_us[-1],color='gray',alpha=0.12)
axA.text(tcap+4,axA.get_ylim()[1]*0.5,'τ > τ_max\n(in-trap blur)',fontsize=8.5,color='gray')
axA.set_title(f'(a) Per-shot signal vs τ   (Δ = {Dref/1e9:.1f} GHz)')
axA.axvline(40,color='0.35',ls='--',lw=1.3,alpha=0.85)
axA.text(41,max(snr_res)*0.30,'40 µs\noperating\npoint',fontsize=8,color='0.3',va='center')
axA.legend(handles=[l1,l2,l3],fontsize=8.6,loc='upper left'); axA.grid(alpha=0.25)

# (b) shot budget vs tau, with N_max*tau (constant) on twin axis
axB.semilogy(tau_us,nm_opt,'o-',color='#8a8a8a',lw=2,label='clean-loss (optimistic)')
axB.semilogy(tau_us,nm_real,'o-',color='#c4161c',lw=2,label='heating + reabsorption (realistic)')
axB.fill_between(tau_us,nm_real,nm_opt,color='#c4161c',alpha=0.08)
axB.set_xlabel('pulse duration  τ (µs)'); axB.set_ylabel('N_max  (30% condensate loss)')
axB.set_ylim(2,300); axB.grid(alpha=0.25,which='both')
axB2=axB.twinx()
axB2.plot(tau_us,nm_real*tau_us,'k--',lw=1.6)
axB2.set_ylabel('N_max × τ  (µs)   [dashed = constant]'); axB2.set_ylim(0,1.6*np.max(nm_real*tau_us))
axB.axvspan(tcap,tau_us[-1],color='gray',alpha=0.12)
for tx,lab in [(40,'40'),(80,'80')]:
    i=int(np.where(tau_us==tx)[0][0]); axB.annotate(f'{nm_real[i]:.0f}',(tx,nm_real[i]),textcoords='offset points',xytext=(0,7),ha='center',fontsize=8.5,color='#c4161c')
axB.axvline(40,color='0.35',ls='--',lw=1.3,alpha=0.85)
axB.set_title('(b) Shot budget vs τ : fewer shots, fixed total integration')
axB.legend(fontsize=8.6,loc='upper right')
plt.tight_layout(); plt.savefig('fig_tausweep.png',dpi=140,bbox_inches='tight',facecolor='white'); plt.show()
