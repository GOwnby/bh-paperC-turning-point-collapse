"""Test (c): convergence of the quantities of interest of the Lagrangian evolution up to t = 60 L/c or loss of ordering.
Initial-data families:  'report' = the report's Newtonian profile (R_j = r_j, m_j) evolved with the GR (Misner-Sharp) equations, which is what the
report's runs did;  'TOV-db' = discretely balanced exact GR equilibrium;  'TOV-pt' = pointwise-sampled TOV.
Shell-crossing criterion: R_{j+1} - R_j <= 0 for some j (loss of ordering of two interface radii) <=> R' = dR/dr <= 0 with label r = initial radius.
Quantities: rho_c(t)/rho_c(0)-1 (centre zone), min_j R' over all zones and over the interior (zones inside 99.9% of the mass), mass fraction outside the
min-R' zone, Lagrangian radii at 50/90/99% of the mass, surface radius, total mass.
The step is fixed per run, dt = T/n with CFL <= 0.4 initially; 'x0.5' means half that step."""
import sys, json, numpy as np
from harness import *

TOUT = 60.0
OUT = {}
def line(d, s):
    return (f"   N={d['N']:5d} dt*{d['dt_scale']:3.1f} {'AV' if d['cq2']>0 else '  '} {s['status']:9s} t_end={s['tend']:7.3f}  t_cross={s['t_cross'] if s['t_cross'] is None else round(s['t_cross'],3)}"
            f"  Rp_min={s['Rpmin']:.4f} (zone {s['kmin']}/{d['N']}, outer mass frac {s['fout']:.2e})  Rp_int={s['Rpmin_int']:.5f}  steps {s['nsteps']}  {s['wall']:.0f}s")

def run_case(kind, data, N, dt_scale=1.0, cq2=0.0, relax=None, T=TOUT, nout=6, adaptive=False, **kw):
    star = make_star(kind, data, N, relax=relax, cq2=cq2, **kw)
    r = (evolve_adaptive if adaptive else evolve)(star, T, nout=nout, dt_scale=dt_scale)
    sn = r['snaps']
    last = sn[-1]
    summ = dict(status=r['status'], tend=last['t'], t_cross=r['t_cross'], Rpmin=min(s['Rpmin'] for s in sn), kmin=sn[-1]['kmin'], fout=sn[-1]['mass_out_at_min'],
                Rpmin_int=min(s['Rpmin_int'] for s in sn), nsteps=r['nsteps'], wall=r['wall'], dt=r['dt'],
                series=[dict(t=s['t'], drho=s['rho_c'] / sn[0]['rho_c'] - 1, Rs=s['Rs'], rl=s['rl'], M=s['M'], Rpmin=s['Rpmin'], Rpmin_int=s['Rpmin_int'], qmax=s['q_max']) for s in sn])
    return summ

def table_at(cases, key_t, kind, field='drho'):
    out = {}
    for k, c in cases.items():
        for s in c['series']:
            if abs(s['t'] - key_t) < 1e-6 or (abs(s['t'] - key_t) < 0.03 * key_t and False):
                out[k] = s[field]
    return out

if __name__ == '__main__':
    kinds = sys.argv[1].split(',') if len(sys.argv) > 1 else ['pure', 'neu']
    Ns = [100, 200, 400, 800, 1600]
    for kind in kinds:
        for fam, data, relax in (('report', 'newton', None), ('TOV-db', 'tov', 'db'), ('TOV-pt', 'tov', None)):
            print(f'\n## {kind} / {fam}: inviscid scheme, N = {Ns}')
            cases = {}
            for N in Ns:
                if fam != 'report' and N > 800: continue
                d = dict(N=N, dt_scale=1.0, cq2=0.0)
                s = run_case(kind, data, N, relax=relax)
                cases[f'N{N}'] = s; print(line(d, s), flush=True)
            OUT[f'{kind}|{fam}|inviscid'] = cases
            if fam == 'report':
                # time-step halving
                for N in (200, 400):
                    d = dict(N=N, dt_scale=0.5, cq2=0.0)
                    s = run_case(kind, data, N, dt_scale=0.5)
                    OUT[f'{kind}|{fam}|dt0.5|N{N}'] = s; print(line(d, s), flush=True)
                # artificial viscosity on
                print(f'## {kind} / {fam}: von Neumann-Richtmyer viscosity q = 2 rho (dU)^2 (compression only)')
                for N in (200, 400, 800, 1600):
                    d = dict(N=N, dt_scale=1.0, cq2=2.0)
                    s = run_case(kind, data, N, cq2=2.0)
                    OUT[f'{kind}|{fam}|AV|N{N}'] = s; print(line(d, s), flush=True)
    json.dump(OUT, open('../results/test_c_converge_%s.json' % ('_'.join(kinds)), 'w'), indent=1, default=float)
