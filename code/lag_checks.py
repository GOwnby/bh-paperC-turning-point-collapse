"""Auxiliary checks on the Lagrangian collapse from the fold, Paper C Sec. 5.7 (the statements that are not tabulated in
Table `lag`).  Sub-commands (run in this directory):
  python3 lag_checks.py resid                 static residual of the initial discrete equilibrium, as a fraction of local gravity,
                                               at several mass fractions and at the outermost nodes (N = 100 ... 1600);
                                               mass of the outermost zone.   (Sec. 5.7 'well-balancing': ~0.5 of local gravity at the outermost node)
  python3 lag_checks.py shells N              t(lambda=0.5), t(0.15) and the plunge for the shells enclosing 2,10,25,50,75,90,100 % of the
                                               rest mass (Sec. 5.7: 'the choice of shell, up to 90 % of the mass, changes the plunge by only 0.0522-0.0541 s'),
                                               and the homologous model (rescaled/raw) from the same hand-off.   (N=200: ~1.5 min)
  python3 lag_checks.py velnorm N mode        mode in {base, div, mul}: hand-off velocity read as dr/dtau (base), dr/dt (div:
                                               lambda_dot r0/alpha), or multiplied by alpha (mul); prints t(lambda_h=0.5).
                                               (Sec. 5.7: reading it as dr/dt shifts t(0.5) by -0.11 %)     (N=200: ~1 min each)
  python3 lag_checks.py stages                fraction of t(lambda_h=0.5) spent at lambda > 0.8 etc., for the homologous model and for the
                                               shipped Lagrangian time series fold_inhom_N<N>_y0.1.npy (N=400,800)  (Sec. 5.7: 'about 85 %')
"""
import sys, io, contextlib, os
import numpy as np
import background as bg; bg.mu_e = 2.0
G, c = bg.G, bg.c
from lag_gr import Lag

def _load():
    with contextlib.redirect_stdout(io.StringIO()):
        import fold_inhom as fi
    return fi

def cmd_resid():
    fi = _load(); rc = fi.rm
    fr = [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.999]
    print('static residual S_j / g_local (g = G m_j/r_j^2) at mass fractions', fr, 'and |S|/g at the last 3 nodes')
    for N in (100, 200, 400, 800, 1600):
        L = Lag(rc, N=N, pert=None); r, U, m = L.unpack(L.y0); g = G*m[1:]/r[1:]**2
        out = [L.S[int(np.argmin(abs(L.mu[1:]/L.mu_tot - f)))]/g[int(np.argmin(abs(L.mu[1:]/L.mu_tot - f)))] for f in fr]
        print(N, ' '.join('%.1e' % x for x in out), '| last3', ' '.join('%.2f' % (L.S[j]/g[j]) for j in (N-3, N-2, N-1)),
              ' mass of outer zone/M %.1e' % (L.dmu[-1]/L.mu_tot), flush=True)

def cmd_shells(N):
    fi = _load(); sc, rm = fi.sc, fi.rm
    th, state = fi.handoff(sc, yh=0.1)
    L = Lag(rm, N=N, pert=(state[0], state[1]))
    fr = [0.02, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]
    idx = [int(np.argmin(abs(L.mu/L.mu_tot - f))) for f in fr]
    print('mass-fraction nodes', [(f, i, round(L.mu[i]/L.mu_tot, 4)) for f, i in zip(fr, idx)])
    y = L.y0.copy(); t = 0.0; cross = {}
    lams = lambda y: np.array([L.unpack(y)[0][i]/L.r0[i] for i in idx])
    prev = lams(y)
    while True:
        dt = L.dt_cfl(y); k1 = L.f(t, y); k2 = L.f(t+dt/2, y+dt/2*k1); k3 = L.f(t+dt/2, y+dt/2*k2); k4 = L.f(t+dt, y+dt*k3)
        yn = y+dt/6*(k1+2*k2+2*k3+k4); cur = lams(yn)
        for v in (0.5, 0.15):
            for k, f in enumerate(fr):
                if (v, f) not in cross and cur[k] < v <= prev[k]: cross[(v, f)] = t+dt*(prev[k]-v)/(prev[k]-cur[k])
        prev = cur; y = yn; t += dt
        r, U, m = L.unpack(y)
        if np.max(2*G*m[1:]/(r[1:]*c**2)) >= 1.0 or t > 10: break
    print('N', N, 't_AH', t)
    for f in fr:
        a = cross.get((0.5, f)); b = cross.get((0.15, f))
        print('  mass frac %.2f: t(0.5)=%.4f t(0.15)=%s plunge=%s' % (f, a, '%.4f' % b if b else None, '%.4f' % (b-a) if b else None))
    for scl in (sc, 1.0):
        t5, t15 = fi.model_from(state, scl); print('  homologous scale %.3f t(0.5)=%.4f t(0.15)=%.4f plunge=%.4f' % (scl, t5, t15, t15-t5))

def cmd_velnorm(N, mode):
    fi = _load(); sc, rm = fi.sc, fi.rm
    th, state = fi.handoff(sc, yh=0.1); lam, ld = state
    L = Lag(rm, N=N, pert=(lam, ld))
    r = L.r0*lam
    r0 = L.r0.copy(); U0 = np.zeros_like(r0); m0 = L.solve_m(r0, U0)
    Gam, rho, P, e = L.zone(r0, U0, m0); h = (e+P)/rho
    alpha_s = (1-2*G*m0[-1]/(r0[-1]*c**2))/Gam[-1]
    az = alpha_s*c**2/h; aj = np.empty(N+1); aj[1:N] = 0.5*(az[1:]+az[:-1]); aj[N] = alpha_s; aj[0] = az[0]
    print('alpha_s %.5f alpha_center %.5f' % (alpha_s, aj[1]))
    if mode == 'div': U = ld*L.r0/aj
    elif mode == 'mul': U = ld*L.r0*aj
    else: U = ld*L.r0
    U[0] = 0
    m = L.solve_m(r, U); L.y0 = np.concatenate([r[1:], U[1:], m[1:]])
    jh = np.searchsorted(L.mu, 0.5*L.mu_tot)
    y = L.y0.copy(); t = 0.0; lp = lam; t5 = None
    while True:
        dt = L.dt_cfl(y, 0.8); k1 = L.f(t, y); k2 = L.f(t+dt/2, y+dt/2*k1); k3 = L.f(t+dt/2, y+dt/2*k2); k4 = L.f(t+dt, y+dt*k3)
        y = y+dt/6*(k1+2*k2+2*k3+k4); t += dt
        l = L.unpack(y)[0][jh]/L.r0[jh]
        if l < 0.5: t5 = t-dt*(0.5-l)/(lp-l); break
        lp = l
    print('N=%d mode=%s t(0.5)=%.5f' % (N, mode, t5))

def cmd_stages():
    fi = _load(); sc = fi.sc
    from fold_dyn import F
    from scipy.integrate import solve_ivp
    th, state = fi.handoff(sc, yh=0.1)
    levels = [0.9, 0.8, 0.7, 0.6, 0.5, 0.3, 0.15]
    def hom(scale):
        evs = [(lambda t, y, v=v: y[0]-v) for v in levels]
        s = solve_ivp(lambda t, y: [y[1], scale*F(y[0])], (0, 10), list(state), method='DOP853', events=evs, rtol=1e-11, atol=1e-14, max_step=0.01)
        return np.array([s.t_events[i][0] for i in range(len(levels))])
    def lag(N):
        rows = np.load(f'fold_inhom_N{N}_y0.1.npy'); t, l = rows[:, 0], rows[:, 1]
        o = []
        for v in levels:
            i = np.argmax(l < v); o.append(np.interp(v, [l[i], l[i-1]], [t[i], t[i-1]]))
        return np.array(o)
    print('levels   ', levels)
    res = {'hom rescaled': hom(sc), 'hom raw': hom(1.0)}
    for N in (400, 800): res['lag N=%d' % N] = lag(N)
    for k, v in res.items(): print('%-13s' % k, np.round(v, 4), ' t(0.8)/t(0.5) = %.3f' % (v[1]/v[4]))

if __name__ == '__main__':
    a = sys.argv[1:]
    if a[0] == 'resid': cmd_resid()
    elif a[0] == 'shells': cmd_shells(int(a[1]))
    elif a[0] == 'velnorm': cmd_velnorm(int(a[1]), a[2])
    elif a[0] == 'stages': cmd_stages()
    else: print(__doc__)
