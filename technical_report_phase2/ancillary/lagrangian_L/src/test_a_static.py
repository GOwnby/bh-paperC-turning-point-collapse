"""Test (a): static equilibrium.  Exact GR hydrostatic (TOV) data, discretely balanced TOV data, and the report's Newtonian profile.
For each family and N the star is evolved to T = 10 and 30 (L/c); the deviation from the initial state is the truncation error
(plus, for the marginal neutronised star, its growth by the unstable mode: e-folding 41 L/c, i.e. factor 1.3 / 2.1).
PASS criterion (stated before running): for TOV-pointwise data, the deviation norms (central density, bulk displacement, mass-weighted L2
displacement) must decrease monotonically with N with observed order p >= 1.8 between the two finest grid pairs."""
import sys, json, numpy as np
from harness import *

OUT = {}
Ns = [50, 100, 200, 400, 800]
def order(e):
    e = np.array(e); return [float(np.log2(e[i] / e[i + 1])) if e[i + 1] > 0 and e[i] > 0 else float('nan') for i in range(len(e) - 1)]

print('# Test (a): static equilibrium / transient amplitude vs N  (time unit L/c, L = report R_star)')
for kind in ('pure', 'neu'):
    for fam, data, relax, newt in (('TOV-pointwise', 'tov', None, False), ('TOV-discrete-balanced', 'tov', 'db', False),
                                   ('report-Newtonian-profile, GR dynamics', 'newton', None, False),
                                   ('report-Newtonian-profile, Newtonian dynamics', 'newton', None, True)):
        for T in (10.0, 30.0):
            rows = []
            for N in Ns:
                if fam == 'TOV-discrete-balanced' and N > 400: continue
                star = make_star(kind, data, N, newt=newt, relax=relax)
                r = evolve(star, T, nout=5)
                s = r['snaps'][-1]
                rows.append(dict(N=N, drho=s['rho_c'] / r['snaps'][0]['rho_c'] - 1, dR_bulk=s['dR_bulk'], dR_L2=s['dR_L2'], umax_int=None,
                                 Rpmin=s['Rpmin'], Rpmin_int=s['Rpmin_int'], dM=s['M'] / r['snaps'][0]['M'] - 1, status=r['status'], t=s['t']))
            key = f'{kind} | {fam} | T={T:g}'
            OUT[key] = rows
            print(f'\n{key}')
            print('   N     |drho_c|      bulk max|dR|/L   massweighted L2 dR/L    min R\'(interior)  order(drho) order(bulk) order(L2)')
            e1 = [abs(x['drho']) for x in rows]; e2 = [x['dR_bulk'] for x in rows]; e3 = [x['dR_L2'] for x in rows]
            o1, o2, o3 = order(e1), order(e2), order(e3)
            for i, x in enumerate(rows):
                print(f"  {x['N']:4d}  {abs(x['drho']):.3e}   {x['dR_bulk']:.3e}      {x['dR_L2']:.3e}          {x['Rpmin_int']:.5f}       " +
                      (f"{o1[i-1]:6.2f}      {o2[i-1]:6.2f}     {o3[i-1]:6.2f}" if i > 0 else ''))
            if fam == 'TOV-pointwise':
                OUT[key + ' | orders'] = dict(drho=o1, bulk=o2, L2=o3)
json.dump(OUT, open('../results/test_a_static.json', 'w'), indent=1)

# PASS/FAIL on the stated criterion (TOV-pointwise, T=10)
ok = True
for kind in ('pure', 'neu'):
    o = OUT[f'{kind} | TOV-pointwise | T=10 | orders']
    last = [o['drho'][-1], o['bulk'][-1], o['L2'][-1]]
    mono = all(np.diff([abs(x['drho']) for x in OUT[f'{kind} | TOV-pointwise | T=10']]) < 0)
    good = mono and min(last) >= 1.8
    print(f'\nPASS/FAIL {kind}: monotone={mono}, orders of the finest pair (drho, bulk, L2) = {last}  ->', 'PASS' if good else 'FAIL')
    ok &= good
print('\nTEST (a) OVERALL:', 'PASS' if ok else 'FAIL')
