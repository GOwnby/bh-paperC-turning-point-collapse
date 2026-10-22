"""Test (j): size of the time variation of Gamma^2 along the shells (Gamma is NOT constant in time, contrary to the frozen-Gamma
assumption of evolve2.py) compared with the relativistic correction 2m/R.  Backs the statement in derivations.md / the report that
delta(Gamma^2) ~ 1e-4 ... 1e-3 during the runs, comparable to 2m/R ~ 2.6e-3.
Report's Newtonian initial data, GR dynamics, N = 400, inviscid, T = 50 L/c (< first edge crossing at ~59), CFL-limited fixed step.
Quantities: max over t and interfaces (enclosing <= 99.9 % of the mass) of |Gamma^2(t) - Gamma^2(0)|; max |U|; max 2m/R over the same interfaces."""
import json, numpy as np
from harness import *

OUT = {}
print('# Test (j): delta(Gamma^2) along the shells, report data, N = 400, T = 50 L/c')
for kind in ('pure', 'neu'):
    star = make_star(kind, 'newton', 400)
    star.start()
    Mb = np.concatenate([[0.0], np.cumsum(star.dM)]); f = Mb / Mb[-1]
    sel = f[1:] <= 0.999
    G2_0 = (star.geo['G'][1:]**2).copy()
    dt0 = star.dt_cfl(star.geo, 0.4)
    nsteps = int(np.ceil(50.0 / dt0)); dt = 50.0 / nsteps
    dG2 = 0.0; umax = 0.0; two_m_over_R = 0.0; dG2_med = 0.0
    for n in range(nsteps):
        star.step(dt)
        g = star.geo
        d = np.abs(g['G'][1:]**2 - G2_0)[sel]
        dG2 = max(dG2, float(d.max())); dG2_med = max(dG2_med, float(np.median(d)))
        umax = max(umax, float(np.max(np.abs(star.U[1:][sel]))))
        two_m_over_R = max(two_m_over_R, float(np.max(2.0 * g['m'][1:][sel] / star.R[1:][sel])))
    OUT[kind] = dict(N=400, T=50.0, nsteps=nsteps, max_dGamma2=dG2, max_of_median_dGamma2=dG2_med, umax=umax, max_2m_over_R=two_m_over_R)
    print(f"  {kind}: steps={nsteps}  max|dGamma^2|={dG2:.2e}  (max over t of the median over interfaces {dG2_med:.2e})  max|U|={umax:.3f} c  max 2m/R={two_m_over_R:.2e}")
json.dump(OUT, open('../results/test_j_gamma_size.json', 'w'), indent=1)
ok = all(1e-5 <= v['max_dGamma2'] <= 1e-2 for v in OUT.values())
print('PASS/FAIL (max |dGamma^2| between 1e-5 and 1e-2, i.e. a non-negligible fraction of 2m/R) ->', 'PASS' if ok else 'FAIL')
