"""Test (g): extension of test (c) to T = 150 L/c (about 3.6 e-foldings of the unstable mode) for the report's initial data, both EOS.
Purpose: determine whether the neutronised configuration crosses at all (test (c) stops at 60), and when; and whether the crossing time converges.
Inviscid scheme N = 200..1600 (stops at the first loss of ordering); the viscous continuation is test (h).  Also the exact (discretely balanced) TOV
equilibria to T = 150 as the control: R' = 1 to rounding, no crossing."""
import sys, json, numpy as np
from harness import *
from test_c_converge import run_case, line

T = 150.0
OUT = {}
for kind in ('pure', 'neu'):
    print(f'\n## {kind}: report data, GR dynamics, T = {T:g}')
    for N in (200, 400, 800, 1600):
        d = dict(N=N, dt_scale=1.0, cq2=0.0)
        s = run_case(kind, 'newton', N, T=T, nout=15); OUT[f'{kind}|inviscid|N{N}'] = s; print(line(d, s), flush=True)
    for N in (200, 400):
        d = dict(N=N, dt_scale=1.0, cq2=0.0)
        s = run_case(kind, 'tov', N, T=T, nout=15, relax='db'); OUT[f'{kind}|TOVdb|N{N}'] = s; print('TOV-db' + line(d, s), flush=True)
json.dump(OUT, open('../results/test_g_long.json', 'w'), indent=1, default=float)
print('done')
