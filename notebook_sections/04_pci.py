# 04 PCI
#
# Exported from: 1 calculations revised 2  multishot  6  extended.ipynb
# This file is a mechanical notebook-section export for refactoring.
# Keep physics, numerical algorithms, parameter values, and execution order equivalent to the notebook.


# %% [cell 15: markdown]
# ## 7. Image formation — phase to camera intensity
#
# This section makes the imaging chain explicit and quantitative: how the column
# phase $\varphi(\mathbf r_\perp)$ becomes the intensity recorded on the camera, for
# each of the three modes, and what the finite numerical aperture and the detector
# do to the measured contrast.
#
# ### 7.1 Decomposition at the Fourier plane
# The field just after the cloud is $E_0\,e^{i\varphi(\mathbf r_\perp)}$ (absorption
# negligible, §5.1). Writing $e^{i\varphi}=1+(e^{i\varphi}-1)$ separates the
# **un-diffracted carrier** "1" — which focuses to the $d_\mathrm{DC}=2.44\lambda f_1/D=6.1\,\mu$m
# DC spot — from the **scattered field** $(e^{i\varphi}-1)$, which spreads across the
# Fourier plane on the scale set by the cloud size ($\gtrsim$ mm). The Fourier-plane
# optic acts only on the DC spot:
#
# * **PCI** — phase dot multiplies the carrier by $t_p e^{i\theta}$:
#   $I/I_0=|t_p e^{i\theta}+(e^{i\varphi}-1)|^2$.
# * **DGI** — opaque stop removes the carrier ($\times10^{-\mathrm{OD}/2}\!\to\!0$):
#   $I/I_0=|e^{i\varphi}-1|^2=4\sin^2(\varphi/2)$.
# * **RAI / clear** — carrier untouched: $I/I_0=|e^{i\varphi}|^2=1$ in the dispersive limit
#   (no absorption signal; this is why RAI needs near-resonant light, not the dispersive probe).
#
# ### 7.2 Transfer curves
# **Two normalisation conventions, stated explicitly** (a frequent source of factor-of-$t_p$
# confusion). Expanding the full PCI expression for small $\varphi$ at $\theta=\pi/2$ gives
# $I/I_0 = t_p^2 + 2(1-\cos\varphi) + 2t_p\sin\varphi \approx t_p^2 + 2t_p\varphi$, i.e.
# relative to the **incident** intensity $I_0$ the linear PCI signal sits on the plate
# pedestal $t_p^2$ with slope $2t_p$. In practice each shot is divided by a no-atom image
# *taken with the plate in beam* (level $t_p^2$), so the **processed ratio** is
# $I/I_{0,\mathrm{plate}} \approx 1 + (2/t_p)\,\varphi$ — slope $2/t_p\approx2.1$, intercept 1.
# The two differ only by the constant $t_p^2$ and the $t_p$ in the slope ($\le$10% for
# $t_p=0.95$); §8 and the simulation use the incident-$I_0$ form ($t_p^2+2t_p\varphi$)
# consistently. The curve below plots the incident-$I_0$ convention, so the linear
# approximation is the genuine tangent to the full curve at $\varphi=0$ (both start at $t_p^2$).
# DGI is intrinsically quadratic and sign-blind (it cannot distinguish $+\varphi$ from
# $-\varphi$). The sign of $\theta$ flips PCI between bright-on-bright and dark-on-bright
# (the protocols' "dark cloud" diagnostic).

# %% [cell 16: code]
phi = np.linspace(0, 2.6, 600)
I_pci_p = np.abs(t_p*np.exp(1j*np.pi/2) + np.exp(1j*phi) - 1)**2
I_pci_m = np.abs(t_p*np.exp(-1j*np.pi/2) + np.exp(1j*phi) - 1)**2
I_lin   = t_p**2 + 2*t_p*phi          # genuine tangent: intercept t_p^2 (plate pedestal), slope 2 t_p
I_dgi   = np.abs(np.exp(1j*phi) - 1)**2
phx = phi_peak(1.5e9, n_col[0])

fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.axvspan(0, 0.5, color='C1', alpha=0.07, label=r'linear regime $\varphi<0.5$')
ax.axhline(t_p**2, color='gray', lw=0.8, ls=':')
ax.annotate(r'plate pedestal $t_p^2$', (1.7, t_p**2-0.18), color='gray', fontsize=8)
ax.plot(phi, I_pci_p, 'C0', lw=2, label=r'PCI full $\theta=+\pi/2$')
ax.plot(phi, I_pci_m, 'C0:', lw=2, label=r'PCI full $\theta=-\pi/2$')
ax.plot(phi, I_lin, 'C1--', lw=1.6, label=r'PCI linear $t_p^2+2t_p\varphi$')
ax.plot(phi, I_dgi, 'C3', lw=2, label=r'DGI $4\sin^2(\varphi/2)$')
ax.axvline(phx, color='gray', ls='--', lw=1)
ax.annotate(r'$\varphi_x$ @ 1.5 GHz', (phx+0.04, 3.6), color='gray', fontsize=9)
ax.set_xlabel(r'phase shift $\varphi$ (rad)'); ax.set_ylabel(r'$I/I_0$')
ax.set_title('Phase-to-intensity transfer (incident-$I_0$ convention)'); ax.set_ylim(0, 4.3)
ax.legend(fontsize=8.5); ax.grid(alpha=0.25); plt.tight_layout(); plt.show()

# %% [cell 17: markdown]
# ### 7.3 Fourier-optics simulation through the real arm
# We propagate the actual TF column-phase profile through the system: forward FFT to
# the Fourier plane, apply (i) the finite-NA pupil of L1 ($\mathrm{NA}=D/2f_1=0.08$) and
# (ii) the Fourier-plane optic on the carrier, inverse FFT to the image plane.
#
# The **blur (contrast-dilution) factor** for imaging along a given axis is a pure
# geometry/NA property: the ratio of the NA-limited centre contrast to the ideal
# (infinite-NA) contrast, evaluated at a small test phase so it is independent of the
# operating $\varphi$. It depends on the **two transverse radii in the image plane**:
#
# * across-cigar (image $x$, along $\hat x$): plane carries $(y,z)=(24.8,1.49)\,\mu$m — one
#   resolved, one not.
# * along-cigar (image along $\hat y$): plane carries $(x,z)=(1.19,1.49)\,\mu$m — **both**
#   below the $3.06\,\mu$m Rayleigh limit, so along-cigar imaging is *more* blurred, not
#   less. It buys a large column density (hence large $\varphi$) at the cost of
#   resolving nothing transversely; it is an atom-number probe, not a structure probe.
#
# The simulation reference background is the plate-only level $t_p^2$ (the no-atom image
# with the phase plate in beam), so the signal is $I_\mathrm{img}-t_p^2$, not $I_\mathrm{img}-1$.

# %% [cell 18: code]
from non_destructive_image import simulate_fourier_image

Ngrid, FOV = 1024, 100e-6
dgrid = FOV/Ngrid
gax = (np.arange(Ngrid)-Ngrid//2)*dgrid
GA, GB = np.meshgrid(gax, gax)
NA = (D_probe/2)/f1
fx = np.fft.fftfreq(Ngrid, dgrid); FX, FY = np.meshgrid(fx, fx)
pupil = (np.sqrt(FX**2+FY**2) <= NA/lam).astype(float)
bg_plate = t_p**2                       # no-atom, plate-in reference level

def _tf_profile(Ra, Rb):
    return np.maximum(0, 1 - GA**2/Ra**2 - GB**2/Rb**2)**1.5

def sim_image(axis, phi_peak_val, mode='PCI', OD=4.0):
    plane = [i for i in range(3) if i != axis]
    prof = _tf_profile(R[plane[0]], R[plane[1]])
    object_field = np.exp(1j*phi_peak_val*prof)
    if mode == 'PCI':   ref = t_p*np.exp(1j*theta)
    elif mode == 'DGI': ref = 10**(-OD/2)
    else:               ref = 1
    return simulate_fourier_image(object_field, pupil, ref), prof

def blur_for_axis(axis, phi_test=0.1):
    Iimg, _ = sim_image(axis, phi_test, 'PCI')
    ideal = np.abs(t_p*np.exp(1j*theta) + (np.exp(1j*phi_test)-1))**2  # infinite-NA, uniform peak
    return (Iimg.max()-bg_plate)/(ideal-bg_plate)

blur_axis = {a: blur_for_axis(a) for a in range(3)}
blur_factor = blur_axis[0]   # across-cigar, used by the reference figures below
print(f"NA = {NA:.3f}, Rayleigh dx = {0.61*lam/NA*1e6:.2f} um")
for a, nm in zip(range(3), ['across (image y,z)', 'along  (image x,z)', 'across (image x,y)']):
    pl = [i for i in range(3) if i != a]
    print(f"  axis {('xyz')[a]} {nm}: transverse R = "
          f"({R[pl[0]]*1e6:.1f}, {R[pl[1]]*1e6:.1f}) um  ->  blur = {blur_axis[a]:.2f}")

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

# %% [cell 35: markdown]
# ### 12.2 $(\Delta, P, \tau)$ optimisation
#
# Maximise usable non-destructive frames at a target per-shot SNR, subject to $\varphi_x<0.5$ (linear
# PCI), $\tau\le\tau_\mathrm{max}$, and the route's power ceiling. Because accumulated SNR is invariant,
# *higher* detuning always buys more frames — but at smaller, less-measurable $\varphi$ and outside the
# AOM-cascade range, so the practical $\Delta$ is set by calibration and phase measurability
# ($\varphi_x\sim0.13$–$0.20$), not by sensitivity.

# %% [cell 36: code]
# ===== CELL H: (Delta, P, tau) optimisation — operating map + invariance frontier + recommendation =====
def _pci_block(Delta_Hz, axis=0, half=1):
    Iimg,_=sim_image(axis, phi_peak(Delta_Hz, n_col[axis]), 'PCI')
    nb=(Ngrid//15)*15; b=Iimg[:nb,:nb].reshape(nb//15,15,nb//15,15).mean(axis=(1,3))
    m=b.shape[0]//2; k=b[m-half:m+half+1,m-half:m+half+1]
    return (k-bg_plate).sum(), k.sum(), k.size
def SNRres_fast(block,P_mW,tau_s):
    S,Sig,n=block; Nd=N_phot_pix(P_mW,tau_s); return S*Nd/np.sqrt(Sig*Nd+n*read_e**2)

fig,(axA,axB)=plt.subplots(1,2,figsize=(13.6,5.4),gridspec_kw=dict(wspace=0.30))

# ---------- (a) operating map at the recommended detuning ----------
Dmap=1.5e9; _blk=_pci_block(Dmap)
Pax=np.linspace(0.5,10,150); Tax=np.linspace(15,tau_max_now*1e6,150); PP,TT=np.meshgrid(Pax,Tax)
NMmap=np.vectorize(lambda P,tu: Nmax_heating(Dmap,P,0.30,tu*1e-6))(PP,TT)
SNmap=np.vectorize(lambda P,tu: SNRres_fast(_blk,P,tu*1e-6))(PP,TT)
DGmap=np.vectorize(lambda P,tu: 0.016*N_phot_pix(P,tu*1e-6))(PP,TT)
im=axA.contourf(PP,TT,np.log10(np.clip(NMmap,1,None)),levels=np.linspace(0,2.4,25),cmap='viridis')
cb=plt.colorbar(im,ax=axA,pad=0.02,ticks=[0,1,2]); cb.set_label('usable frames  N_max  (log)'); cb.set_ticklabels(['1','10','100'])
cs=axA.contour(PP,TT,SNmap,levels=[8,12,20],colors='white',linewidths=1.2)
axA.clabel(cs,fmt='SNR/res=%d',fontsize=8,colors='white',manual=[(2.0,18),(2.5,33),(3.0,76)])
cd=axA.contour(PP,TT,DGmap,levels=[3*read_e],colors='#ff5555',linewidths=2.2,linestyles='dashed')
axA.clabel(cd,fmt='DGI quantitative',fontsize=8.5,colors='#ff3333',manual=[(6.2,28)])
mk=[('balanced\n3.5 mW · 40 µs',3.5,40,(16,16)),('dynamics\n3.5 · 15',3.5,15,(16,2)),('hi-SNR\n5 · 80',5.0,80,(16,2))]
_bb=dict(boxstyle='round,pad=0.25',fc='white',ec='0.5',alpha=0.85)
for lab,P,tu,off in mk:
    axA.plot(P,tu,'o',color='k',ms=10,mfc='white',mew=1.8,zorder=6)
    axA.annotate(lab,(P,tu),xytext=off,textcoords='offset points',fontsize=8.4,fontweight='bold',
                 color='#202020',ha='left',va='center',bbox=_bb,zorder=7)
axA.set_xlabel('probe power  P (mW)'); axA.set_ylabel('pulse duration  τ (µs)')
axA.set_title(f'(a) Operating map at Δ = {Dmap/1e9:.1f} GHz\nframes (colour) · PCI SNR/res (white) · DGI floor (red)')

# ---------- (b) accumulated-SNR invariance: the frames-vs-SNR curve is the same for every detuning ----------
Dscan=np.array([1.0,1.5,2.0,3.0,4.0])*1e9; _bk={D:_pci_block(D) for D in Dscan}
tt=np.linspace(15,tau_max_now*1e6,60)
curves=[]
for D in Dscan:
    sn=np.array([SNRres_fast(_bk[D],3.5,t*1e-6) for t in tt]); nm=np.array([Nmax_heating(D,3.5,0.30,t*1e-6) for t in tt])
    curves.append((sn,nm))
# common SNR grid -> band (min/max N_max across detunings) to show they coincide
sgrid=np.linspace(3,28,120); band=[]
for s in sgrid:
    vals=[np.interp(s,c[0],c[1],left=np.nan,right=np.nan) for c in curves]
    vals=[v for v in vals if np.isfinite(v)]; band.append((min(vals),max(vals)) if vals else (np.nan,np.nan))
band=np.array(band)
axB.fill_between(sgrid,band[:,0],band[:,1],color='tab:purple',alpha=0.25,label='Δ = 1–4 GHz (all curves)')
axB.plot(curves[1][0],curves[1][1],'-',color='tab:purple',lw=2.4,label='Δ = 1.5 GHz (representative)')
axB.set_yscale('log'); axB.set_xlim(3,28); axB.set_ylim(2,300)
axB.axvline(8,color='gray',ls=':',lw=1.2); axB.annotate('SNR/res = 8',(8.3,2.6),fontsize=8.5,color='gray')
axB.annotate('the band is razor-thin:\nframes vs per-shot SNR is the SAME\nfor every detuning — accumulated\nSNR is set by the destruction\nbudget, not by (Δ, P, τ)',
             xy=(12,np.interp(12,curves[1][0],curves[1][1])),xytext=(15,70),fontsize=8.6,color='#4a2d6b',
             arrowprops=dict(arrowstyle='->',color='#4a2d6b',lw=1.2))
axB.set_xlabel('per-shot PCI SNR / resolution element'); axB.set_ylabel('usable frames  N_max  (realistic)')
axB.set_title('(b) Why τ only re-slices a fixed budget  (P = 3.5 mW)')
axB.legend(fontsize=8.6,loc='upper right'); axB.grid(alpha=0.25,which='both')
plt.tight_layout(); plt.savefig('fig_optimisation.png',dpi=140,bbox_inches='tight',facecolor='white'); plt.show()

# ---------- recommended operating points ----------
print("RECOMMENDED OPERATING POINTS  (across-cigar x; realistic heating+reabs budget; 30% condensate loss)")
hdr=f"{'mode':<24}{'Δ/GHz':>6}{'P/mW':>6}{'τ/µs':>6}{'phi_x':>7}{'N_max':>6}{'SNR/res':>8}{'SNR/pix':>8}{'DGI/e-':>7}"
print(hdr); print('-'*len(hdr))
for lab,Dg,P,tus in [("balanced (Route A)",1.5,3.5,40),("dynamics: many frames",1.5,3.5,15),
                     ("hi-SNR / quant. DGI",1.5,5.0,80),("higher-Δ, more frames",2.5,5.0,40),
                     ("max range (beat-lock)",3.0,8.0,40)]:
    D=Dg*1e9; ts=tus*1e-6; nm=Nmax_heating(D,P,0.30,ts); dgi=0.016*N_phot_pix(P,ts)
    print(f"{lab:<24}{Dg:>6.1f}{P:>6.1f}{tus:>6d}{phi_peak(D,n_col[0]):>7.3f}{nm:>6.0f}"
          f"{SNRres_fast(_pci_block(D),P,ts):>8.1f}{SNR_pixel(D,P,0,ts):>8.1f}{dgi:>7.0f}")
print("\nRule of thumb: τ set by the read-noise floor (DGI quantitative for τ >~ 40-50 µs at 3.5-5 mW);")
print("Δ set by calibration range + keeping phi_x measurable (~0.13-0.20); P by the chosen source route.")
