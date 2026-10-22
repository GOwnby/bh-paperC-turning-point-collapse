"""Test (h): collapse metrics that can be converged.  Report initial data, GR dynamics, von Neumann-Richtmyer AV (q = 2 rho dU^2), which removes the
edge shock-crossing (see test (g)) so that the evolution continues until the CORE collapses.
Metrics: first time at which rho_c/rho_c(0) reaches 1.5, 2, 3, 5, 10 (log-linear interpolation of the output series, 0.5 L/c spacing), and t_cross (loss of ordering).
Pure Chandrasekhar: run to T = 400 L/c; neutronised: to the core crossing (T <= 120)."""
import sys, json, numpy as np
from harness import *
from test_c_converge import run_case, line

OUT = {}
THR = (1.5, 2.0, 3.0, 5.0, 10.0, 30.0, 100.0, 300.0)
def hits(series):
    t = np.array([s['t'] for s in series]); r = np.array([1 + s['drho'] for s in series])
    out = {}
    for th in THR:
        k = np.where(r >= th)[0]
        if len(k) == 0 or k[0] == 0: out[th] = None; continue
        k = k[0]
        out[th] = float(np.interp(np.log(th), [np.log(r[k-1]), np.log(r[k])], [t[k-1], t[k]]))
    return out
def do(kind, T, nout, N, ds, adaptive=False):
    d = dict(N=N, dt_scale=ds, cq2=2.0)
    s = run_case(kind, 'newton', N, T=T, nout=nout, cq2=2.0, dt_scale=ds, adaptive=adaptive)
    h = hits(s['series']); s['hits'] = {str(k): v for k, v in h.items()}
    OUT[f'{kind}|AV{"ad" if adaptive else ""}|' + (f'dt{ds:g}|' if ds != 1.0 else '') + f'N{N}'] = s
    print(line(d, s)); print('      first time rho_c/rho_c0 >= ' + ', '.join(f'{k:g}: ' + ('-' if v is None else f'{v:.3f}') for k, v in h.items()),
                            f'   max drho = {max(x["drho"] for x in s["series"]):.3f}', flush=True)
print('\n## neu: report data, AV q=2, T = 120.  (i) fixed step dt = T/n from the CFL bound of the INITIAL state (x factor)')
for N, ds in ((400, 1.0), (400, 0.5), (400, 0.25)): do('neu', 120.0, 240, N, ds)
print('## neu: (ii) the step re-limited every step by the CFL bound of the CURRENT state (dt = 0.4 x factor x dt_CFL(t)); stops when rho_c = 300 rho_c0')
for N, ds in ((200, 1.0), (400, 1.0), (800, 1.0), (1600, 1.0), (400, 0.5)): do('neu', 120.0, 240, N, ds, adaptive=True)
print('\n## pure: report data, AV q=2, T = 300.  (i) fixed step')
for N, ds in ((200, 1.0), (400, 1.0), (800, 1.0), (400, 0.5)): do('pure', 300.0, 600, N, ds)
print('## pure: (ii) adaptive step')
for N, ds in ((400, 1.0),): do('pure', 300.0, 600, N, ds, adaptive=True)
json.dump(OUT, open('../results/test_h_collapse.json', 'w'), indent=1, default=float)
print('done')
