"""Eigenvalues of the discrete Jacobian of the Lagrangian scheme (lag_gr.Lag), Paper C Sec. 5.7 'Validation against the linear results'.
Usage: python3 lag_eigs.py N [N ...]      (e.g. 100 200 400 800; ~1 min each, N=800 ~5 min)
Prints per N: (i) unstable growth rate at the fiducial rho_c=5e13 kg/m^3; (ii) the lowest discrete sigma^2 at the turning
point rho_c = rho_max (from fold_step1.json) and the second radial mode omega_2 there; (iii) the stable-branch fundamental
frequency at rho_c=1e13 kg/m^3 and the extreme real parts of its Jacobian eigenvalues (should vanish: no damping, no growth)."""
import sys, json, os
import numpy as np, background as bg
bg.mu_e = 2.0
from lag_gr import Lag
from lag_jac import jac_eigs
HERE = os.path.dirname(os.path.abspath(__file__))
rm = json.load(open(os.path.join(HERE, 'fold_step1.json')))['rho_max']
for N in [int(a) for a in sys.argv[1:]] or [100, 200, 400]:
    ev, _ = jac_eigs(Lag(5e13, N=N, pert=None))
    print(f"N={N}  fiducial 5e13: growth rate = {ev.real.max():.4f} s^-1", flush=True)
    ev, _ = jac_eigs(Lag(rm, N=N, pert=None))
    ok = np.abs(ev) > 1e-6; om2 = -(ev**2).real
    idx = np.where(ok)[0][np.argsort(np.abs(om2[ok]))[:4]]
    pos = np.sort(om2[ok][om2[ok] > 0])
    print(f"N={N}  turning point: lowest sigma^2 = {om2[idx[0]]:.5f} s^-2 ; next modes omega^2 = {pos[pos > 10][:2]} (omega_2 = {np.sqrt(pos[pos > 10][0]):.2f} rad/s)", flush=True)
    ev, _ = jac_eigs(Lag(1e13, N=N, pert=None))
    im = np.abs(ev.imag); im = np.sort(im[im > 1e-3])
    print(f"N={N}  stable 1e13: fundamental omega = {im[0]:.5f} rad/s ; max Re(eigenvalue) = {ev.real.max():+.2e} s^-1, min = {ev.real.min():+.2e} s^-1", flush=True)
