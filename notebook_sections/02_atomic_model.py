# 02 Atomic Model
#
# Exported from: 1 calculations revised 2  multishot  6  extended.ipynb
# This file is a mechanical notebook-section export for refactoring.
# Keep physics, numerical algorithms, parameter values, and execution order equivalent to the notebook.


# %% [cell 7: markdown]
# ## 4. Condensate profile (Thomas-Fermi)
# $\mu=\tfrac{\hbar\bar\omega}{2}(15N_0a_s/a_\mathrm{ho})^{2/5}$,
# $n_\mathrm{peak}=\mu m/(4\pi\hbar^2a_s)$, $R_i=\sqrt{2\mu/m\omega_i^2}$,
# $\tilde n_i^{(\mathrm{peak})}=\tfrac43 n_\mathrm{peak}R_i$.
#
# **[rev] §4.3 condensate-fraction check.** $T_\mathrm{cloud}=200$ nK is consistent with a
# condensate only if $T<T_c\simeq0.94\,\hbar\bar\omega N_\mathrm{tot}^{1/3}/k_B$. Solving
# $N_0=N_\mathrm{tot}[1-(T/T_c)^3]$ self-consistently fixes $N_\mathrm{tot}$ and the
# condensate fraction, to be checked against K24 §6.3.1.

# %% [cell 8: code]
from non_destructive_image import build_thomas_fermi_state

tf_initial = build_thomas_fermi_state(
    atom_number=N0,
    scattering_length=a_s,
    trap_frequencies_hz=trap_Hz,
    atomic_mass=m,
    hbar=hbar,
    boltzmann_constant=kB,
)
omega = tf_initial.trap_angular_frequencies
omega_bar = tf_initial.geometric_mean_frequency
a_ho = tf_initial.harmonic_oscillator_length
mu = tf_initial.chemical_potential
T_mu = tf_initial.chemical_potential_temperature
n_peak = tf_initial.peak_density
R = tf_initial.radii
n_col = tf_initial.column_density
N_check = tf_initial.atom_number_check

print(f"omega_bar/2pi = {omega_bar/(2*np.pi):.1f} Hz   a_ho = {a_ho*1e6:.3f} um")
print(f"mu/kB    = {T_mu*1e9:.1f} nK")
print(f"n_peak   = {n_peak*1e-6:.3e} cm^-3")
print(f"R(x,y,z) = ({R[0]*1e6:.2f}, {R[1]*1e6:.2f}, {R[2]*1e6:.2f}) um")
print(f"N check  = {N_check:.2e}  (input {N0:.2e})")
for ax_, nc in zip('xyz', n_col):
    print(f"   ncol along {ax_}: {nc*1e-4:.3e} cm^-2")

from scipy.optimize import brentq
def _frac_residual(Ntot):
    Tc = 0.94*hbar*omega_bar/kB * Ntot**(1/3)
    f  = 1-(T_cloud/Tc)**3 if Tc > T_cloud else 0.0
    return Ntot*f - N0                  # condensed atoms must equal N0
N_tot_sc = brentq(_frac_residual, 3e4, 1e7)
Tc_sc = 0.94*hbar*omega_bar/kB * N_tot_sc**(1/3)
print(f"\n[rev] self-consistent: N_tot = {N_tot_sc:.2e}, T_c = {Tc_sc*1e9:.0f} nK, "
      f"condensate fraction = {N0/N_tot_sc:.2f}")
print(f"[rev] => bimodal fit runs against a large thermal pedestal (~{1-N0/N_tot_sc:.0%}); cross-check K24 6.3.1")
print(f"[rev] mu/E_rec = {T_mu/T_rec:.3f}  ->  one scattered photon delivers {T_rec/T_mu:.1f}x mu")
