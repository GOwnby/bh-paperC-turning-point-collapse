"""Test (d): sensitivity to the treatment of the vacuum edge (outer boundary).
 variants (all N = 400 unless stated; report initial data, GR dynamics; T = 60 or loss of ordering):
   baseline         free surface: p = 0 outside the last zone, no floor, no artificial viscosity
   floor f          zone pressures clipped from below at p_floor = f * p_c   (the 'surface pressure floor' of the report)
   AV               von Neumann-Richtmyer viscosity (q = 2 rho dU^2, compression)
   last zone x w    mass of the outermost zone multiplied by w
   exterior p_ext   a constant exterior pressure p_ext = f p_c acting on the last interface
 and the first-step kick: max |U_j| over the 5 outermost interfaces after ONE time step, in units of g dt (g = local gravity)."""
import sys, json, numpy as np
from harness import *
from lag import Lag

OUT = {}
def one(kind, N, label, T=60.0, relax=None, **kw):
    S = Setup(kind)
    d = S.newton_data(N)
    wlast = kw.pop('wlast', 1.0)
    dM = d['dM'].copy(); dM[-1] *= wlast
    star = S.make(dict(d, dM=dM), **{k: v for k, v in kw.items() if k not in ('floor', 'pext')})
    star.start()
    pc = star.geo['p'][0]
    if 'floor' in kw: star.p_floor = kw['floor'] * pc
    if 'pext' in kw: star.p_ext = kw['pext'] * pc
    # first-step kick
    star2 = S.make(dict(d, dM=dM), **{k: v for k, v in kw.items() if k not in ('floor', 'pext')})
    star2.start()
    if 'floor' in kw: star2.p_floor = kw['floor'] * pc
    if 'pext' in kw: star2.p_ext = kw['pext'] * pc
    star2.start()
    dt = star2.dt_cfl(star2.geo, 0.4)
    g = star2.geo['m'][-1] / star2.R[-1]**2
    star2.step(dt)
    kick = np.max(np.abs(star2.U[-5:])) / (g * dt)
    r = evolve(star, T, nout=6)
    sn = r['snaps']
    row = dict(label=label, N=N, status=r['status'], t_end=sn[-1]['t'], t_cross=r['t_cross'], Rpmin=min(s['Rpmin'] for s in sn),
               kmin=sn[-1]['kmin'], fout=sn[-1]['mass_out_at_min'], drho30=None, drho_end=sn[-1]['rho_c'] / sn[0]['rho_c'] - 1, kick=float(kick),
               Rs_end=sn[-1]['Rs'])
    for s in sn:
        if abs(s['t'] - 30.0) < 1e-6: row['drho30'] = s['rho_c'] / sn[0]['rho_c'] - 1
    print(f"  {label:28s} {r['status']:9s} t_end={row['t_end']:7.3f} t_cross={'-' if row['t_cross'] is None else '%.3f' % row['t_cross']:>7s} Rp_min={row['Rpmin']:+.3f} (zone {row['kmin']}/{N}, outer mass frac {row['fout']:.1e})  "
          f"drho_c(30)={row['drho30'] if row['drho30'] is None else round(row['drho30'], 5)}  first-step kick = {kick:.3f} g dt")
    return row

for kind in ('pure', 'neu'):
    print(f'\n## {kind}: N = 400, report initial data')
    rows = []
    rows.append(one(kind, 400, 'baseline'))
    for f in (1e-14, 1e-12, 1e-10, 1e-8, 1e-6):
        rows.append(one(kind, 400, f'floor {f:.0e} p_c', floor=f))
    for f in (1e-12, 1e-8):
        rows.append(one(kind, 400, f'exterior p_ext {f:.0e} p_c', pext=f))
    rows.append(one(kind, 400, 'AV q=2 rho dU^2', cq2=2.0))
    rows.append(one(kind, 400, 'AV quad 2 + linear 0.3', cq2=2.0, cl=0.3))
    for w in (0.25, 0.5, 2.0, 4.0):
        rows.append(one(kind, 400, f'last zone mass x {w:g}', wlast=w))
    OUT[kind] = rows
json.dump(OUT, open('../results/test_d_edge.json', 'w'), indent=1, default=float)
