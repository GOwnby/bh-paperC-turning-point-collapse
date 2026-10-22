"""(iii) Reproduction and diagnosis of the re-creation ancillary/lagrangian_envelope.py (cubic-spline derivatives, RK4, dt = 1e-3).

1. Fast exact replica: the cubic spline derivative at the knots is a LINEAR operator, D = CubicSpline(r, I)(r, 1) (N x N), so
   derivs(y) = D @ y reproduces CubicSpline(r, y)(r, 1) to rounding error; the right-hand side, the RK4 loop and the stopping rule
   are otherwise the original ones.  We verify the replica against the ORIGINAL right-hand side (executed from the original source text)
   to rounding error, and the first 300 RK4 steps of both.
2. Reproduce the quoted crossings (pure: t = 14.78 at shell 0; neutronised: t = 2.19 at the outermost shell).
3. Exhibit the discretisation error: (a) residual acceleration of the equilibrium initial data (should be 0), (b) node-number dependence
   of the crossing time/shell, (c) removal of the smooth 'weight' that switches the pressure force off where rho+p < 1e-6 max(rho),
   (d) the free-fall estimate of the outermost node, (e) the same initial data in the converged staggered scheme.
"""
import sys, os, io, time, contextlib, json, re
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import cumulative_trapezoid

HERE = os.path.dirname(os.path.abspath(__file__)); REF = os.path.join(HERE, '..', 'ref')
sys.path.insert(0, REF)
with contextlib.redirect_stdout(io.StringIO()):
    import eos as E, initial_profile as IPp, initial_profile_neu as IPn, eos_neutronized as EN


class SplineScheme:
    def __init__(self, kind, N=300, hole=0.01, mask=True, a_power=2):
        self.kind, self.N, self.mask, self.a_power = kind, N, mask, a_power
        IP = IPn if kind == 'neu' else IPp
        self.P_of_rho = EN.P_of_rho if kind == 'neu' else IPp.P_of_rho
        sol = IP.build_star(5e12, rmax=2e6, n=8000)
        r_raw, m_raw = sol.t, sol.y[0]
        self.r_surf = r_surf = r_raw[-1]
        r_grid = np.linspace(hole * r_surf, r_surf, N)
        m_grid = np.interp(r_grid, r_raw, m_raw)
        self.L = L = r_surf
        self.r_lab = r_grid / L
        self.m0 = (m_grid * E.G / E.c**2) / L
        self.Gamma_shell = np.sqrt(np.maximum(1 - 2 * self.m0 / self.r_lab, 1e-10))
        self.Gscale = E.G * L**2 / E.c**2
        self.D = CubicSpline(self.r_lab, np.eye(N))(self.r_lab, 1)       # derivative operator (not-a-knot, as the original default)
        self.M_Msun = m_raw[-1] / E.Msun

    def p_dimless(self, rho_d):
        return self.Gscale * (self.P_of_rho(np.maximum(rho_d, 0) / self.Gscale) / E.c**2)

    def rhs(self, t, y):
        N = self.N; D = self.D
        R, m, U = y[0:N], y[N:2 * N], y[2 * N:3 * N]
        Rp = D @ R; mp = D @ m
        rho = np.maximum(mp / (4 * np.pi * np.maximum(R, 1e-12)**2 * np.where(np.abs(Rp) > 1e-12, Rp, 1e-12)), 0)
        p_d = self.p_dimless(rho)
        dp_dr = D @ p_d
        denom = rho + p_d
        threshold = 1e-6 * np.max(rho)
        denom_safe = np.maximum(denom, 1e-300)
        w = denom**2 / (denom**2 + threshold**2) if self.mask else np.ones(N)
        integrand = w * (-dp_dr / denom_safe)
        a = np.clip(np.exp(cumulative_trapezoid(integrand[::-1], self.r_lab[::-1], initial=0.0)[::-1]), 1e-6, 1e2)
        X = Rp / np.maximum(self.Gamma_shell, 1e-12)
        pgrad = w * (self.Gamma_shell / denom_safe) * (dp_dr / np.where(np.abs(X) > 1e-12, X, 1e-12))
        dU = -a**self.a_power * ((m + 4 * np.pi * R**3 * p_d) / np.maximum(R, 1e-12)**2 + pgrad)
        return np.concatenate([a * U, -4 * np.pi * a * p_d * R**2 * U, dU]), rho, Rp, a, pgrad

    def initial_y(self):
        return np.concatenate([self.r_lab.copy(), self.m0.copy(), np.zeros(self.N)])

    def run(self, dt=1e-3, tmax=30.0, log_every=500, verbose=False, stop=True):
        N = self.N
        y = self.initial_y(); t = 0.0; n = 0
        _, rho0, Rp0, _, _ = self.rhs(0.0, y)
        hist = []
        status = 'tmax'
        while t < tmax:
            k1, rho, Rp, a, pg = self.rhs(t, y)
            if n % log_every == 0:
                hist.append((t, rho[0] / rho0[0], float(Rp.min()), int(Rp.argmin())))
                if verbose: print(f'  t={t:8.3f}  rho_c/rho_c0={rho[0]/rho0[0]:.5f}  min R\'={Rp.min():.4f} (shell {Rp.argmin()})', flush=True)
            if (np.min(Rp) < 0 or not np.all(np.isfinite(y))) and stop:
                status = 'crossing'; break
            k2 = self.rhs(t + dt / 2, y + dt / 2 * k1)[0]; k3 = self.rhs(t + dt / 2, y + dt / 2 * k2)[0]; k4 = self.rhs(t + dt, y + dt * k3)[0]
            y = y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4); t += dt; n += 1
        i = int(np.argmin(Rp))
        return dict(status=status, t=t, shell=i, r_over_R=float(self.r_lab[i]), rho_over_rhoc=float(rho[i] / rho0[0]), drho_c=float(rho[0] / rho0[0] - 1),
                    hist=hist, n=n)


def original_rhs(kind):
    """execute the head of the ORIGINAL lagrangian_envelope.py (everything before the time loop) and return its namespace"""
    src = open(os.path.join(REF, 'lagrangian_envelope.py')).read()
    head = src.split("y = np.concatenate([r_lab.copy()")[0]
    ns = {'__name__': 'orig'}
    old = sys.argv; sys.argv = ['lagrangian_envelope.py', kind, '1e-3']
    cwd = os.getcwd(); os.chdir(REF)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(head, 'lagrangian_envelope_head', 'exec'), ns)
    finally:
        sys.argv = old; os.chdir(cwd)
    return ns


if __name__ == '__main__':
    out = {}
    print('# 1. replica vs original right-hand side')
    for kind in ('pure', 'neu'):
        sch = SplineScheme(kind)
        ns = original_rhs(kind)
        y = sch.initial_y()
        f_new = sch.rhs(0.0, y)[0]; f_old = ns['rhs'](0.0, y)[0]
        rng = np.random.default_rng(0)
        y2 = y * (1 + 1e-4 * rng.standard_normal(y.size)); y2[2 * sch.N:] = 1e-4 * rng.standard_normal(sch.N)
        g_new = sch.rhs(0.0, y2)[0]; g_old = ns['rhs'](0.0, y2)[0]
        d1 = np.max(np.abs(f_new - f_old)) / np.max(np.abs(f_old)); d2 = np.max(np.abs(g_new - g_old)) / np.max(np.abs(g_old))
        # RK4 steps
        ya = y.copy(); yb = y.copy(); dt = 1e-3; t = 0.0
        for n in range(300):
            ka = [None] * 4
            ka[0] = sch.rhs(t, ya)[0]; ka[1] = sch.rhs(t + dt / 2, ya + dt / 2 * ka[0])[0]; ka[2] = sch.rhs(t + dt / 2, ya + dt / 2 * ka[1])[0]; ka[3] = sch.rhs(t + dt, ya + dt * ka[2])[0]
            ya = ya + dt / 6 * (ka[0] + 2 * ka[1] + 2 * ka[2] + ka[3])
            kb = [None] * 4
            kb[0] = ns['rhs'](t, yb)[0]; kb[1] = ns['rhs'](t + dt / 2, yb + dt / 2 * kb[0])[0]; kb[2] = ns['rhs'](t + dt / 2, yb + dt / 2 * kb[1])[0]; kb[3] = ns['rhs'](t + dt, yb + dt * kb[2])[0]
            yb = yb + dt / 6 * (kb[0] + 2 * kb[1] + 2 * kb[2] + kb[3]); t += dt
        d3 = np.max(np.abs(ya - yb))
        print(f'  {kind}: rhs(initial) rel.diff {d1:.1e}; rhs(perturbed) rel.diff {d2:.1e}; after 300 RK4 steps max|y_replica - y_original| = {d3:.1e}   {"PASS" if max(d1,d2)<1e-8 and d3<1e-9 else "FAIL"}')
        out[f'replica_{kind}'] = dict(d_rhs0=d1, d_rhs1=d2, d_300=d3)

    print('\n# 2. reproduction of the re-creation (dt = 1e-3, N = 300, R'' < 0 stop)')
    for kind, tmax in (('neu', 6.0), ('pure', 20.0)):
        sch = SplineScheme(kind)
        t0 = time.time()
        r = sch.run(dt=1e-3, tmax=tmax, log_every=500)
        print(f'  [{kind}] {r["status"]} at t = {r["t"]:.3f} (L/c) = {r["t"]*sch.L/E.c*1e3:.2f} ms; shell {r["shell"]} of {sch.N} (r/R = {r["r_over_R"]:.3f}, rho/rho_c = {r["rho_over_rhoc"]:.3e}); '
              f'central density changed by {100*r["drho_c"]:+.2f} %; wall {time.time()-t0:.0f} s')
        amp = max(abs(h[1] - 1) for h in r['hist'])
        print(f'      max |rho_c/rho_c0 - 1| over the logged history: {amp:.3f}')
        out[f'repro_{kind}'] = dict(t=r['t'], shell=r['shell'], r_over_R=r['r_over_R'], drho_c=r['drho_c'], max_amp=amp, hist=r['hist'][::max(1, len(r['hist']) // 60)])

    print('\n# 3a. residual acceleration of the equilibrium initial data (exact equilibrium: dU/dt = 0 up to discretisation error)')
    for kind in ('pure', 'neu'):
        sch = SplineScheme(kind)
        f, rho, Rp, a, pg = sch.rhs(0.0, sch.initial_y())
        dU = f[2 * sch.N:]
        g = sch.m0 / sch.r_lab**2
        rel = dU / g
        idx = [0, 1, 2, 5, 10, 30, 100, 200, 280, 290, 295, 297, 298, 299]
        print(f'  {kind}: (dU/dt)/g at nodes {idx}:')
        print('     ', np.array2string(rel[idx], precision=3))
        print(f'      max |dU/dt| = {np.max(np.abs(dU)):.3e} (g at the surface = {g[-1]:.3e})')
        d = rho + sch.p_dimless(rho); w = d**2 / (d**2 + (1e-6 * rho.max())**2)
        print(f'      number of nodes with weight w < 0.99: {int(np.sum(w<0.99))} (indices {np.where(w<0.99)[0].tolist()[:8]}), rho/rho_c there: {np.array2string(rho[w<0.99][:4]/rho.max(), precision=2)}')
        out[f'residual_{kind}'] = dict(rel=rel[idx].tolist(), idx=idx)

    print('\n# 3b. node-number dependence of the spline scheme (dt scaled as 1e-3*300/N)')
    rows = []
    for kind, tmax in (('neu', 6.0), ('pure', 20.0)):
        for N in (150, 300, 600):
            sch = SplineScheme(kind, N=N)
            dt = 1e-3 * 300 / N
            t0 = time.time()
            r = sch.run(dt=dt, tmax=tmax, log_every=10**9)
            print(f'  [{kind}] N={N}: {r["status"]} t = {r["t"]:.3f}, shell {r["shell"]} (r/R = {r["r_over_R"]:.3f}), drho_c = {100*r["drho_c"]:+.2f} %  ({time.time()-t0:.0f} s)')
            rows.append(dict(kind=kind, N=N, t=r['t'], shell=r['shell'], r_over_R=r['r_over_R'], status=r['status']))
    out['node_dependence'] = rows

    print('\n# 3c. spline scheme with the pressure-force weight removed (w = 1) and with a^2 -> a (consistent proper-time form)')
    rows = []
    for kind, tmax in (('neu', 6.0), ('pure', 20.0)):
        for label, kw in (('no weight mask', dict(mask=False)), ('a^1 instead of a^2', dict(a_power=1)), ('no mask and a^1', dict(mask=False, a_power=1))):
            sch = SplineScheme(kind, **kw)
            r = sch.run(dt=1e-3, tmax=tmax, log_every=10**9)
            print(f'  [{kind}] {label}: {r["status"]} t = {r["t"]:.3f}, shell {r["shell"]} (r/R = {r["r_over_R"]:.3f}), drho_c = {100*r["drho_c"]:+.2f} %')
            rows.append(dict(kind=kind, variant=label, t=r['t'], shell=r['shell'], status=r['status'], drho_c=r['drho_c']))
    out['variants'] = rows

    print('\n# 3d. free-fall estimate for the outermost node when its pressure force is switched off')
    for kind in ('pure', 'neu'):
        sch = SplineScheme(kind)
        dr = sch.r_lab[-1] - sch.r_lab[-2]
        g = sch.m0[-1] / sch.r_lab[-1]**2
        print(f'  {kind}: spacing {dr:.5f}, g = {g:.4e}  ->  sqrt(2 dr/g) = {np.sqrt(2*dr/g):.3f} L/c  (crossing of the outermost pair if the outer node alone falls freely)')
        out[f'freefall_{kind}'] = float(np.sqrt(2 * dr / g))
    json.dump(out, open(os.path.join(HERE, '..', 'results', 'repro_spline.json'), 'w'), indent=1, default=float)
