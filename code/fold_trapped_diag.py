"""Diagnostics at the first trapped shell of the Lagrangian collapse from the fold (Paper C Sec. 5.7, 'Through the model's wall').
Usage: python3 fold_trapped_diag.py N [y_handoff=0.1] [cfl=0.4]      (N=200: ~1 min; 400: ~2 min; 800: ~5 min; 1600: ~12 min)
Same hand-off as fold_inhom.py / fold_to_BH.py; integrates (RK4) until max_r 2Gm/(r c^2) >= 1 and reports, linearly
interpolated to max compactness = 1: time, m_AH/M, central rest-mass density amplification (zone average), radii of the shells
enclosing 5,20,50,80,95,99 % of the rest mass and of the surface relative to their initial radii, the Hamiltonian-constraint
violation (max over the run while lambda_h>0.2, and at the trapped shell), t(lambda_h=0.5), t(0.15) and the time elapsed
between lambda_h=0.5 and the trapped shell.  Appends one JSON line to fold_trapped_N<N>_cfl<cfl>.json."""
import sys, json, time, os
import numpy as np
import background as bg; bg.mu_e = 2.0
from fold_inhom import handoff, sc, rm
from lag_gr import Lag
G, c = bg.G, bg.c
N = int(sys.argv[1]); YH = float(sys.argv[2]) if len(sys.argv) > 2 else 0.1; CFL = float(sys.argv[3]) if len(sys.argv) > 3 else 0.4
th, state = handoff(sc, yh=YH)
L = Lag(rm, N=N, pert=(state[0], state[1]))
y = L.y0.copy(); t = 0.0; k = 0; t0 = time.time(); jh = np.searchsorted(L.mu, 0.5*L.mu_tot)
cross = {}; viol_max = 0.0; prev = None
def snapshot(y, t):
    r, U, m = L.unpack(y); Gam, rho, P, e = L.zone(r, U, m)
    comp = 2*G*m[1:]/(r[1:]*c**2); jm = int(np.argmax(comp)) + 1
    V = 4*np.pi/3*(r[1:]**3 - r[:-1]**3)
    viol = np.max(np.abs(np.diff(m) - e*V/c**2))/m[-1]
    return dict(t=t, lam=r[jh]/L.r0[jh], mc=comp.max(), mAH=m[jm]/m[-1], rhoc=rho[0]/L.rho_c, viol=viol, r=r.copy())
while True:
    s = snapshot(y, t)
    for v in (0.5, 0.15):
        if v not in cross and s['lam'] < v and prev is not None:
            cross[v] = float(np.interp(v, [s['lam'], prev['lam']], [s['t'], prev['t']]))
    if s['lam'] > 0.2: viol_max = max(viol_max, s['viol'])
    if s['mc'] >= 1.0 and prev is not None:
        w = (1.0 - prev['mc'])/(s['mc'] - prev['mc']); lerp = lambda a, b: a + w*(b - a)
        mu_frac = L.mu/L.mu_tot; fr = {}
        rr = lambda sn: sn['r']/np.where(L.r0 > 0, L.r0, 1)
        for q in (0.05, 0.2, 0.5, 0.8, 0.95, 0.99):
            j = np.searchsorted(mu_frac, q); fr[str(q)] = float(lerp(rr(prev)[j], rr(s)[j]))
        out = dict(N=N, y_handoff=YH, cfl=CFL, t_AH=lerp(prev['t'], s['t']), m_AH_over_M=lerp(prev['mAH'], s['mAH']), rhoc_amp=lerp(prev['rhoc'], s['rhoc']),
                   shell_radius_ratio=fr, surface_ratio=float(lerp(rr(prev)[-1], rr(s)[-1])), viol_max_lam_gt_0p2=viol_max, viol_at_AH=lerp(prev['viol'], s['viol']),
                   t_lam05=cross.get(0.5), t_lam015=cross.get(0.15), dt_lam05_to_AH=lerp(prev['t'], s['t']) - cross.get(0.5), steps=k, wall=time.time() - t0)
        print(json.dumps(out), flush=True)
        open(f'fold_trapped_N{N}_cfl{CFL}.json', 'a').write(json.dumps(out) + '\n'); break
    if np.isnan(y).any() or t > 40 or np.min(np.diff(s['r'])) <= 0:
        print("aborted", t, s['mc']); break
    prev = s
    dt = L.dt_cfl(y, CFL)
    k1 = L.f(t, y); k2 = L.f(t + dt/2, y + dt/2*k1); k3 = L.f(t + dt/2, y + dt/2*k2); k4 = L.f(t + dt, y + dt*k3)
    y = y + dt/6*(k1 + 2*k2 + 2*k3 + k4); t += dt; k += 1
