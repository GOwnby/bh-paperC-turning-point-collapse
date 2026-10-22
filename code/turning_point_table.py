"""Turning-point test of Paper C (Sec. 3.6): for each equation of state, the central density at the mass maximum dM/drho_c = 0
(maximiser tolerance 1e-8 in log10 rho_c) and the zero of the 5-function Rayleigh-Ritz sigma^2 (12000 radial points,
brentq tolerance 1e7 kg/m^3 = 4e-7 relative).  Also: grid dependence of sigma^2 at rho_max and at rho_c=5e13, the basis-size
dependence, and the regression test of the 'first implementation' (first-law energy density off, 8 pi G term off / sign flipped).
Usage: python3 turning_point_table.py [table|grid|first]   (default: all; ~15-25 min)"""
import sys, functools
import numpy as np, background as bg, tov_solve, metric, gr_rayleigh_ritz as rr
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar, brentq
cases = {'ideal (mu_e=2)': (0, 2.0), 'He4': (2, 2.0), 'C12': (6, 2.0), 'O16': (8, 2.0), 'Fe56': (26, 56/26)}
orig_solve = tov_solve.solve_tov
def sigma2(rc, nb=5, n=12000):
    metric.solve_tov = functools.partial(orig_solve, n_eval=n)
    return eigh(*rr.build_RR_matrices(metric.build_metric(rc), nb)[:2])[0][0]
def mass(lr, n=12000):
    metric.solve_tov = functools.partial(orig_solve, n_eval=n); return -metric.build_metric(10**lr)['M']
what = sys.argv[1:] or ['table', 'grid', 'first']
if 'table' in what:
    print("EOS | rho_c(dM/drho=0) [1e13 kg/m3] | rho_c(sigma^2=0) [1e13] | difference [%] | M_max [Msun]")
    for name, (Z, mue) in cases.items():
        bg.Zlat = float(Z); bg.mu_e = mue
        res = minimize_scalar(mass, bounds=(13.2, 13.6), method='bounded', options={'xatol': 1e-8}); rp = 10**res.x
        rz = brentq(sigma2, 0.9*rp, 1.1*rp, xtol=1e7)
        print(f"{name} | {rp/1e13:.5f} | {rz/1e13:.5f} | {100*(rz-rp)/rp:.4f} | {-res.fun/bg.M_sun:.5f}", flush=True)
if 'grid' in what:
    bg.Zlat = 0.0; bg.mu_e = 2.0
    rm = 10**minimize_scalar(mass, bounds=(13.2, 13.6), method='bounded', options={'xatol': 1e-8}).x
    print("sigma^2(rho_max), n=5:", {n: round(float(sigma2(rm, 5, n)), 5) for n in (3000, 12000, 48000)}, flush=True)
    print("sigma^2(rho_max) vs basis size (12000 pts):", {nb: round(float(sigma2(rm, nb)), 5) for nb in (3, 4, 5, 6, 7)}, flush=True)
    print("sigma^2(5e13), n=5, 3000/6000/12000/24000 pts:", [round(float(sigma2(5e13, 5, n)), 3) for n in (3000, 6000, 12000, 24000)], flush=True)
    z3 = brentq(lambda rc: sigma2(rc, 5, 3000), 0.9*rm, 1.1*rm, xtol=1e7)
    print("zero crossing with 3000 points: %.5e ; relative offset from rho_max = %.3f %%" % (z3, 100*(z3 - rm)/rm), flush=True)
if 'first' in what:
    bg.Zlat = 0.0; bg.mu_e = 2.0
    for eps_flag, f8 in ((False, 0.0), (True, 0.0), (False, 1.0), (True, 1.0), (True, -1.0)):
        bg.use_eps = eps_flag; rr.INCLUDE_8PI = f8
        z = brentq(lambda rc: sigma2(rc), 0.8e13, 6e13, xtol=1e8)
        pk = minimize_scalar(mass, bounds=(13.0, 14.0), method='bounded', options={'xatol': 1e-6})
        print(f"use_eps={eps_flag}  8piG-term factor={f8:+.0f}: sigma^2=0 at {z:.4e}  ; mass peak at {10**pk.x:.4e} (M_max={-pk.fun/bg.M_sun:.4f}) ; sigma^2(5e13)={sigma2(5e13):.3f}", flush=True)
