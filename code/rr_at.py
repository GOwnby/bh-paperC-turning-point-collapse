"""Rayleigh-Ritz (n = 5) sigma^2 and, for sigma^2 > 0, the fundamental angular frequency, at the central densities given on the command line
[kg/m^3].  Used for the 'stable branch' value quoted in Sec. 5.8 (rho_c = 1e10 g/cm^3 = 1e13 kg/m^3: omega = 2.802 rad/s).
usage: python rr_at.py 1e13 5e13"""
import sys, numpy as np
from scipy.linalg import eigh
import background as bg
from metric import build_metric
from gr_rayleigh_ritz import build_RR_matrices
for a in sys.argv[1:]:
    rc = float(a)
    bgd = build_metric(rc)
    A, B, _, _ = build_RR_matrices(bgd, 5)
    s2 = eigh(A, B)[0][0]
    msg = f'rho_c = {rc:.4e}: M = {bgd["M"]/bg.M_sun:.5f} Msun, sigma^2 = {s2:+.5f} s^-2'
    msg += f', omega = {np.sqrt(s2):.4f} rad/s' if s2 > 0 else f', growth rate = {np.sqrt(-s2):.4f} 1/s'
    print(msg, flush=True)
