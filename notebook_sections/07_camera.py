# 07 Camera
#
# Exported from: 1 calculations revised 2  multishot  6  extended.ipynb
# This file is a mechanical notebook-section export for refactoring.
# Keep physics, numerical algorithms, parameter values, and execution order equivalent to the notebook.


# %% [cell 20: code]
Mag = f2/f1; pix_obj = pix_cam/Mag

def N_phot_pix(P_mW, tau_s=None, QE=None):
    # Detected photons per camera pixel at I/I0 = 1, using the probe intensity AT the
    # cloud (beam centre). Consistent with intensity_at_atoms used for scattering.
    tau_set = 100e-6
    if tau_s is None: tau_s = tau_set
    if QE is None: QE = QE_cam
    return intensity_at_atoms(P_mW) * pix_obj**2 * tau_s * QE / E_phot

I_img_pci, prof0 = sim_image(0, phx, 'PCI')
I_img_dgi, _     = sim_image(0, phi_peak(1.5e9, n_col[0]), 'DGI')

def to_camera(Iratio, P_mW, QE=None):
    nb = (Ngrid//15)*15
    binned = Iratio[:nb,:nb].reshape(nb//15,15,nb//15,15).mean(axis=(1,3))
    Nd = N_phot_pix(P_mW, QE=QE)
    counts = rng.poisson(np.clip(binned,0,None)*Nd) + rng.normal(0, read_e, binned.shape)
    return counts/Nd, binned

cam_pci, ideal_pci = to_camera(I_img_pci, 2.0)
cam_dgi, ideal_dgi = to_camera(I_img_dgi, 5.0)
ext = [-FOV/2*1e6, FOV/2*1e6, -FOV/2*1e6, FOV/2*1e6]

fig, axs = plt.subplots(2, 2, figsize=(11, 6.4))
im0 = axs[0,0].imshow(cam_pci, extent=ext, origin='lower', cmap='inferno', vmin=0.85, vmax=1.15)
axs[0,0].set_title('PCI camera frame (2 mW)'); plt.colorbar(im0, ax=axs[0,0], fraction=0.03, label='I/I0')
im1 = axs[0,1].imshow(cam_dgi, extent=ext, origin='lower', cmap='inferno', vmin=-0.01, vmax=0.05)
axs[0,1].set_title('DGI camera frame (5 mW)'); plt.colorbar(im1, ax=axs[0,1], fraction=0.03, label='I/I0')
for a in axs[0]:
    a.set_xlim(-45,45); a.set_ylim(-12,12); a.set_xlabel('y (um)'); a.set_ylabel('z (um)')

nb = Ngrid//15; ycam = (np.arange(nb)-nb//2+0.5)*pix_obj*1e6; mid = nb//2
axs[1,0].plot(ycam, ideal_pci[mid], 'C0', lw=2, label='noiseless (NA-filtered)')
axs[1,0].plot(ycam, cam_pci[mid], 'C0.', ms=4, alpha=0.7, label='single shot + noise')
axs[1,0].axhline(bg_plate, color='gray', ls=':', lw=1, label=r'plate background $t_p^2$')
axs[1,0].set_xlim(-45,45); axs[1,0].set_xlabel('y (um)'); axs[1,0].set_ylabel('I/I0')
axs[1,0].set_title('PCI lineout along cigar'); axs[1,0].legend(fontsize=8); axs[1,0].grid(alpha=0.25)

zf = gax*1e6
axs[1,1].plot(zf, I_img_pci[:, Ngrid//2], 'C0', lw=1.8, label='NA = 0.08')
axs[1,1].plot(zf, bg_plate+2*t_p*phx*np.maximum(0,1-(gax)**2/R[2]**2)**1.5, 'k--', lw=1, label='infinite NA')
axs[1,1].set_xlim(-8,8); axs[1,1].set_xlabel('z (um)'); axs[1,1].set_ylabel('I/I0')
axs[1,1].set_title('Across-cigar cut (unresolved)'); axs[1,1].legend(fontsize=8); axs[1,1].grid(alpha=0.25)
plt.tight_layout(); plt.show()

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
