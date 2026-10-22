"""Fig. 1 of Paper C: fundamental-mode eigenvalue of the validated Galerkin solver (exact l=0 self-gravity closure,
surface boundary term) across uniform-Gamma1 polytropes, against the Pekeris threshold Gamma1 = 4/3.
Uses the Newtonian solver of the supplementary technical report (../technical_report_phase2/ancillary/linear_perturbation).
Writes ../fig1_pekeris_validation.png and prints the table.  Runtime: ~20 s."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'technical_report_phase2', 'ancillary', 'linear_perturbation'))
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.linalg import eigh
from galerkin_fem import assemble
from validate_polytrope import build_polytrope, G

def count_nodes(xi):
    s = np.sign(xi[2:-2]); s = s[s != 0]
    return int(np.sum(np.abs(np.diff(s)) > 0))

def fundamental(Gamma, K=1e10, rho_c=1e3, N=300):
    sol = build_polytrope(K, Gamma, rho_c=rho_c, rmax=5e8)
    r_raw, rho_raw = sol.t, np.maximum(sol.y[1], 1e-10)
    r = np.concatenate([[0.0], np.linspace(r_raw[1], r_raw[-1]*0.9995, N-1)])
    rho0 = np.interp(r, r_raw, rho_raw); rho0[0] = rho_raw[0]
    P0 = K*rho0**Gamma; Gamma1 = np.full(len(r), Gamma)
    M, D = assemble(r, rho0, P0, Gamma1, G=G, include_self_gravity=True)
    ev, vec = eigh(M[1:, 1:], D[1:, 1:])
    for idx in np.argsort(ev)[:25]:
        if count_nodes(vec[:, idx]) == 0:
            return ev[idx]
    return np.nan

gam = np.array([1.24, 1.2667, 1.30, 1.3167, 4/3, 1.3667, 1.40, 1.4333, 1.4667, 1.50, 1.55, 1.60, 5/3])
w2 = np.array([fundamental(g) for g in gam])
for g, w in zip(gam, w2): print(f"Gamma1={g:.4f}  fundamental omega^2 = {w:+.5e} s^-2  {'unstable' if w < 0 else 'stable'}")
fig, ax = plt.subplots(figsize=(6.2, 4.2))
m = gam < 4/3 - 1e-9
ax.plot(gam[m], w2[m], 'o-', c='#b02a18', label=r'unstable ($\omega^2<0$)')
ax.plot(gam[~m], w2[~m], 'o-', c='#1f5fb4', label=r'stable ($\omega^2>0$)')
ax.axhline(0, ls='--', c='gray', lw=0.8); ax.axvline(4/3, ls=':', c='gray', lw=0.8); ax.text(1.35, 0.9*np.nanmax(w2), r'$\Gamma_1=4/3$')
ax.set_xlabel(r'adiabatic index $\Gamma_1$'); ax.set_ylabel(r'fundamental-mode $\omega^2$  [s$^{-2}$]')
ax.set_title('Galerkin solver versus the classical Pekeris criterion'); ax.grid(alpha=.3); ax.legend(loc='lower right', fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(HERE, '..', 'fig1_pekeris_validation.png'), dpi=200)
