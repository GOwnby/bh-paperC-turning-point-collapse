"""Test (f) [diagnostic, not a convergence test]: the early evolution of the report's initial data is a LINEAR response about the GR equilibrium with the same zone contents.
The report's data (Newtonian profile) = discrete GR equilibrium R_eq (Newton solve at the same zone masses dM) + displacement d0 = R_report - R_eq.
Linear theory:  d^2 dR/dt^2 = A dR, A = diag(a) J (J = dF/dR at R_eq) => dR(t) = V diag(cos(omega_n t) | cosh(sigma t)) V^-1 d0.
We compare the central density from this prediction with the nonlinear evolution (N = 200 and 400, inviscid), and report the share of the unstable mode.
NOTE: the displacement d0 is NOT small (it changes the central density by tens of per cent because the star is nearly neutrally stable), so the linear prediction is only an order-of-magnitude guide; its purpose is to quantify how far the report's data are from the GR equilibrium with the same zone masses."""
import sys, json, numpy as np
from scipy.linalg import eig, inv
from harness import *
from lag import relax_equilibrium, fd_jacobian

OUT = {}
print('# initial acceleration of the report data and of the pointwise TOV data in the staggered scheme: (dU/dt)/g with g = m/R^2, interfaces enclosing 1 %..99 % of the mass')
for kind in ('pure', 'neu'):
    for data in ('newton', 'tov'):
        for N in (200, 400, 800):
            st = make_star(kind, data, N); st.start()
            F, a = st.accel(st.geo)
            g = st.geo['m'][1:] / st.R[1:]**2
            Mb = np.cumsum(st.dM) / st.dM.sum()
            sel = (Mb > 0.01) & (Mb < 0.99)
            rel = F[sel] / g[sel]
            print(f'  {kind} {data:6s} N={N:4d}: median (dU/dt)/g = {np.median(rel):+.5f}, max |.| = {np.max(np.abs(rel)):.5f}')
            OUT[f'acc {kind} {data} {N}'] = dict(median=float(np.median(rel)), maxabs=float(np.max(np.abs(rel))))
for kind in ('pure', 'neu'):
    for N in (200, 400):
        star = make_star(kind, 'newton', N)
        R_rep = star.R.copy()
        star.start(); a0 = star.a.copy()
        eqs = make_star(kind, 'newton', N)
        res, it = relax_equilibrium(eqs)
        R_eq = eqs.R.copy()
        d0 = (R_rep - R_eq)[1:]
        x = R_eq[1:].copy()
        J = fd_jacobian(eqs, x)
        F, a = eqs.force_R(x)
        A = a[:, None] * J
        w, V = eig(A); w = w.real; V = V.real
        order = np.argsort(-w)          # most unstable first (largest lambda)
        w, V = w[order], V[:, order]
        q0 = np.linalg.solve(V, d0)
        def rho_c_of(dR):
            Rn = R_eq.copy(); Rn[1:] += dR
            geo = eqs.geometry(Rn, np.zeros(N + 1), Gguess=np.ones(N + 1))
            return geo['e'][0]
        e_eq = rho_c_of(np.zeros(N))
        times = [5.0, 10.0, 20.0, 30.0, 40.0]
        pred = []; pred_no_unstable = []
        for t in times:
            sw = np.sqrt(np.abs(w))
            qt = np.where(w > 0, q0 * np.cosh(sw * t), q0 * np.cos(sw * t))
            pred.append(rho_c_of(V @ qt) / e_eq)
            qt2 = qt.copy(); qt2[w > 0] = 0.0
            pred_no_unstable.append(rho_c_of(V @ qt2) / e_eq)
        # nonlinear evolution (inviscid) of the report data: rho_c relative to the report data's own initial value
        sn = evolve(make_star(kind, 'newton', N), 40.0, nout=8)['snaps']
        evo = {round(s['t'], 6): s['rho_c'] / sn[0]['rho_c'] for s in sn}
        e_rep0 = rho_c_of(d0) / e_eq                     # central density of the report data relative to the GR equilibrium with the same zones
        starT = make_star(kind, 'tov', N, relax='db'); starT.start(); e_tov = starT.geo['e'][0]
        eq_vs_tov = e_eq / e_tov - 1                     # central density of that GR equilibrium relative to the TOV star of rho_c = 5e12
        cs_ = (2.998e8 / star.S.L)
        print(f'\n{kind} N={N}: Newton residual {res:.1e}; number of unstable modes {int(np.sum(w>0))}; lowest omega^2 = {-w[0]*cs_**2:+.3f} s^-2 (e-folding {1/np.sqrt(w[0]):.1f} L/c if unstable); '
              f'central density of the report data / GR equilibrium with the same zone masses - 1 = {e_rep0-1:+.5f}; that equilibrium vs the TOV star (5e12) = {eq_vs_tov:+.5f}; '
              f'|q0| of the unstable mode / norm(d0) = {abs(q0[0])*np.linalg.norm(V[:,0])/np.linalg.norm(d0):.3f}')
        print('   t      rho_c/rho_c(0)-1: nonlinear evolution | linear mode-sum | linear without unstable mode')
        for t, p, pn in zip(times, pred, pred_no_unstable):
            ev = evo.get(t, None)
            print(f'  {t:5.1f}   {ev/1-1 if ev else float("nan"):+.5f} | {p/e_rep0-1:+.5f} | {pn/e_rep0-1:+.5f}')
        OUT[f'{kind} {N}'] = dict(times=times, evolution=[evo.get(t) for t in times], linear=[p / e_rep0 for p in pred], linear_no_unstable=[p / e_rep0 for p in pred_no_unstable], e_rep0=e_rep0, eq_vs_tov=eq_vs_tov, w2_low=float(-w[0] * cs_**2),
                                  n_unstable=int(np.sum(w > 0)), lam0=float(w[0]), q0_frac=float(abs(q0[0]) * np.linalg.norm(V[:, 0]) / np.linalg.norm(d0)))
json.dump(OUT, open('../results/test_f_forced.json', 'w'), indent=1, default=float)
