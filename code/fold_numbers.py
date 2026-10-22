"""Derived numbers of the fold-passage sections of Paper C (Secs. 5.6, 5.7, 5.9), computed from the exact-GR coefficients in
fold_step1.json / fold_comp.json and the accretion rates; no new physics.  Usage: python3 fold_numbers.py  (~1 s)
 (1) beta, alpha of the scaling for Mdot = 1e-9 ... 1e-5 Msun/yr (table `fold`, columns beta and alpha) and 3.41167 beta;
 (2) dissipation examples of Sec. 5.7: Mdot = 1e-7, nu = 10 s^-1 and Mdot = 1e-5, nu = 100 s^-1: nu~, t_collapse, Delta M in the overdamped regime;
 (3) composition table `foldcomp` (Sec. 5.9): beta and t_coll = 3.41167 beta at Mdot = 1e-7 for the five equations of state, relative changes of
     kappa, a, t and alpha with respect to the ideal gas.
"""
import json, numpy as np
yr = 3.15576e7
s_star = 3.41167; tau_star = 2.338107
d = json.load(open('fold_step1.json')); p = d['pts'][:3]
kap = np.mean([q['kappa'] for q in p]); a = np.mean([q['a'] for q in p])
print('exact-GR fold coefficients (mean of y = 0.002, 0.005, 0.01): kappa = %.4f s^-2, a = %.5f Msun (paper 18.733, 4.027e-3)' % (kap, a))
print('rho_max = %.5e kg/m^3, M_max = %.5f Msun' % (d['rho_max'], d['Mmax']))
bet = lambda Md, kap=kap, a=a: (4*a/(kap**2*Md/yr))**0.2
print('\n(1) Mdot[Msun/yr]  beta[s]  alpha  3.41167*beta[s]')
for Md in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
    b = bet(Md); print('    %.0e  %7.2f  %.2e  %8.1f' % (Md, b, 2/(kap*b**2), s_star*b))
print('\n(2) dissipation examples (overdamped Airy law t = 2.338 beta nu~^(2/3), Delta M = Mdot t)')
for Md, nu in ((1e-7, 10.0), (1e-5, 100.0)):
    b = bet(Md); nt = nu*b; t = tau_star*b*nt**(2/3); dM = Md/yr*t
    dM2 = tau_star*(Md/yr)**(2/3)*(4*a/kap**2)**(1/3)*nu**(2/3)
    print('    Mdot = %.0e, nu = %.0f: beta = %.1f s, nu~ = %.0f, t_collapse = %.3g s, Delta M = %.2e Msun (closed form %.2e)' % (Md, nu, b, nt, t, dM, dM2))
print('\n(3) composition (Sec. 5.9), Mdot = 1e-7 Msun/yr')
fc = json.load(open('fold_comp.json'))
rows = fc.items() if isinstance(fc, dict) else [(r['name'], r) for r in fc]
base = None
for k, v in rows:
    b = bet(1e-7, v['kappa'], v['a']); t = s_star*b; al = 2/(v['kappa']*b**2)
    if base is None: base = (v['kappa'], v['a'], t, al)
    print('    %-8s rho_max = %.3e g/cc  Mmax = %.4f  kappa = %.2f  a = %.3f e-3  beta = %.1f s  t_coll = %.0f s   d(kappa) %+5.1f%%  d(a) %+5.1f%%  d(t) %+5.1f%%  d(alpha) %+5.1f%%' % (
        k, v['rho_max']/1e3, v['Mmax'], v['kappa'], v['a']*1e3, b, t, 100*(v['kappa']/base[0]-1), 100*(v['a']/base[1]-1), 100*(t/base[2]-1), 100*(al/base[3]-1)))
