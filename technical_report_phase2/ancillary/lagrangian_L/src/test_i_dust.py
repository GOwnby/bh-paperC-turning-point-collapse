"""Test (i): comparison with EXACT general-relativistic solutions (no free parameters, no fitting).

Pressureless dust is the one case of the Misner-Sharp system with a closed-form solution for arbitrary initial density profile
(Oppenheimer-Snyder 1939 for a homogeneous ball; Lemaitre 1933 - Tolman 1934 - Bondi 1947 for an inhomogeneous ball): with p = 0 every shell
is an independent marginally-bound/bound geodesic,  m(r) is constant along the shell and, for a shell starting at rest at areal radius R0 with
enclosed mass M = m(R0), the proper-time evolution is the cycloid
       R = (R0/2)(1 + cos eta),      tau = (1/2) sqrt(R0^3/(2M)) (eta + sin eta)        (G = c = 1),
valid until the shell reaches R = 0 (eta = pi) or crosses another shell.  The scheme of lag.py is run with a dust "equation of state"
(p = 0, e = rho_N) and a FIXED step; the shells are placed with the DISCRETE constraint m_j = sum e V and the DISCRETE Gamma_j^2 = 1 - 2 m_j/R_j,
so the discrete initial data are exact.  Everything that is then compared (R_j(t), U_j(t), and m_j(t) which must stay constant) is produced by the
scheme's evolution equations (1)-(3), (5) of lag.py, including the relativistic Gamma-dependence, without any pressure support.
Strong field: initial compactness 2M/R0 = 0.6 (homogeneous) and 0.5 (inhomogeneous), i.e. far from the Newtonian regime.

PASS criterion (fixed before running): the max relative error of R_j(t) over all shells and over t <= 0.8 t_collapse(innermost shell)
decreases monotonically with N (with dt = t_ref/(4N)) with observed order >= 1.8 between the two finest grids, and m_j drifts by < 1e-3 of m_j
on the finest grid.
"""
import sys, json
import numpy as np
from lag import Lag

FOUR_PI = 4.0 * np.pi


class Dust:
    """barotropic 'equation of state' p = 0, e = rho_N"""
    def thermo(self, rN):
        rN = np.asarray(rN, float)
        return rN.copy(), np.zeros_like(rN), np.zeros_like(rN)


def build(N, ek_fun, rb=1.0):
    """zones k=0..N-1 on a uniform interface grid R_j = rb j/N; energy density e_k = ek_fun(r_centre)."""
    R = rb * np.arange(N + 1) / N
    V = FOUR_PI / 3.0 * (R[1:]**3 - R[:-1]**3)
    e = ek_fun(0.5 * (R[1:] + R[:-1]))
    m = np.concatenate([[0.0], np.cumsum(e * V)])                      # discrete constraint (2), exact
    G = np.ones(N + 1); G[1:] = np.sqrt(1.0 - 2.0 * m[1:] / R[1:])      # discrete definition of Gamma_j at rest
    Gz = 0.5 * (G[:-1] + G[1:])
    dM = e * V / Gz                                                    # so that rho_N Gamma_k / ... reproduces e_k exactly
    return R, dM, m


def eta_of_tau(tau, R0, M):
    """invert tau = (1/2) sqrt(R0^3/(2M)) (eta + sin eta) by Newton/bisection (monotone for eta in [0, pi])"""
    A = 0.5 * np.sqrt(R0**3 / (2.0 * M))
    y = tau / A                                                         # = eta + sin eta
    eta = np.minimum(y / 2.0, np.pi)                                    # starting point (exact at 0 and pi)
    for _ in range(100):
        f = eta + np.sin(eta) - y
        fp = 1.0 + np.cos(eta)
        step = f / np.maximum(fp, 1e-12)
        eta = np.clip(eta - step, 0.0, np.pi)
        if np.max(np.abs(f)) < 1e-15: break
    return eta


def exact(R0, M, tau):
    """exact R and U = dR/dtau for the cycloid"""
    eta = eta_of_tau(tau, R0, M)
    R = 0.5 * R0 * (1.0 + np.cos(eta))
    A = 0.5 * np.sqrt(R0**3 / (2.0 * M))
    U = -0.5 * R0 * np.sin(eta) / (A * (1.0 + np.cos(eta)))
    return R, U


def run_case(name, ek_fun, Ns, frac=0.8, **kw):
    res = []
    for N in Ns:
        R0, dM, m0 = build(N, ek_fun)
        # collapse time of the innermost shell (exact)
        M1, Rb1 = m0[1], R0[1]
        tcol = 0.5 * np.sqrt(Rb1**3 / (2.0 * M1)) * np.pi
        # all shells: tau_collapse_j = (pi/2) sqrt(R_j^3/(2 m_j)); the inhomogeneous case has the smallest at the centre (e decreasing)
        tcol_all = 0.5 * np.pi * np.sqrt(R0[1:]**3 / (2.0 * m0[1:]))
        tmin = tcol_all.min()
        T = frac * tmin
        nsteps = 4 * N
        dt = T / nsteps
        star = Lag(Dust(), dM, R0, **kw)
        star.start()
        eR = 0.0; eU = 0.0; dm = 0.0; mworst = 0.0
        for n in range(nsteps):
            star.step(dt)
            R_ex, U_ex = exact(R0[1:], m0[1:], star.t)
            eR = max(eR, float(np.max(np.abs(star.R[1:] - R_ex) / R_ex)))
            eU = max(eU, float(np.max(np.abs(star.U[1:] - U_ex) / np.maximum(np.abs(U_ex), 1e-3))))
            geo = star.geometry(star.R, star.U)
            dm = max(dm, float(np.max(np.abs(geo['m'][1:] / m0[1:] - 1.0))))
        res.append(dict(N=N, nsteps=nsteps, T=T, t_over_tcol_min=T / tmin, err_R=eR, err_U=eU, drift_m=dm, status='ok'))
    return res


def orders(vals):
    return [float(np.log2(vals[i] / vals[i + 1])) if vals[i + 1] > 0 and vals[i] > 0 else float('nan') for i in range(len(vals) - 1)]


if __name__ == '__main__':
    Ns = [25, 50, 100, 200, 400]
    ok = True
    OUT = {}
    # (i) homogeneous ball (Oppenheimer-Snyder): 2M/R0 = 0.6   -> e0 = 3 M / (4 pi R0^3), M = 0.3 R0  (rb = 1)
    e_hom = 3.0 * 0.3 / FOUR_PI
    # (ii) inhomogeneous ball (LTB): e(r) = e_c (1 - 0.5 (r/rb)^2), chosen so that 2 m(rb)/rb = 0.5
    # m(rb) = 4 pi e_c (rb^3/3 - 0.5 rb^3/5) = 4 pi e_c rb^3 (0.3333 - 0.1) -> e_c = 0.25 / (4 pi * 0.233333)
    e_c = 0.25 / (FOUR_PI * (1.0 / 3.0 - 0.1))
    cases = [('homogeneous (Oppenheimer-Snyder), 2M/R0 = 0.6', lambda r: e_hom * np.ones_like(r)),
             ('inhomogeneous (Lemaitre-Tolman-Bondi), 2M/R0(surface) = 0.5, e_c/e_surf = 2', lambda r: e_c * (1.0 - 0.5 * r**2))]
    print('# Test (i): exact dust solutions (Oppenheimer-Snyder / LTB cycloids); fixed dt = T/(4N), T = 0.8 t_collapse(innermost shell)')
    for name, f in cases:
        rows = run_case(name, f, Ns)
        OUT[name] = rows
        oR = orders([r['err_R'] for r in rows]); oU = orders([r['err_U'] for r in rows])
        print(f'\n{name}\n   N    nsteps  max rel err R(t)   order   max rel err U(t)   order   max |m/m0-1|')
        for i, r in enumerate(rows):
            print(f"  {r['N']:4d}  {r['nsteps']:5d}   {r['err_R']:.3e}      " + (f"{oR[i-1]:5.2f}" if i else '     ') +
                  f"    {r['err_U']:.3e}      " + (f"{oU[i-1]:5.2f}" if i else '     ') + f"    {r['drift_m']:.2e}")
        mono = all(np.diff([r['err_R'] for r in rows]) < 0)
        good = mono and min(oR[-2:]) >= 1.8 and rows[-1]['drift_m'] < 1e-3
        print(f"PASS/FAIL {name}: monotone={mono}, orders of the two finest pairs = {oR[-2:]}, finest drift in m = {rows[-1]['drift_m']:.1e} ->", 'PASS' if good else 'FAIL')
        OUT[name + ' | orders'] = dict(R=oR, U=oU)
        ok &= good

    # (iii) why the Gamma coupling is iterated to convergence: N = 50 homogeneous ball, steps = 100 ... 1600, fixed count of 3 iterations vs converged
    print('\n(iii) homogeneous ball, N = 50: max rel err of R(t) versus the number of time steps, fixed 3 iterations of the rho_N-m-Gamma coupling vs iteration to 1e-13')
    e_h = lambda r: e_hom * np.ones_like(r)
    def err_steps(ns, **kw):
        R0, dM, m0 = build(50, e_h)
        tmin = (0.5 * np.pi * np.sqrt(R0[1:]**3 / (2.0 * m0[1:]))).min(); T = 0.8 * tmin; dt = T / ns
        star = Lag(Dust(), dM, R0, **kw); star.start(); e = 0.0
        for n in range(ns):
            star.step(dt); Rx, _ = exact(R0[1:], m0[1:], star.t)
            e = max(e, float(np.max(np.abs(star.R[1:] - Rx) / Rx)))
        return e
    steps = [100, 200, 400, 800, 1600]
    e3 = [err_steps(n, gtol=0.0, n_iter=3) for n in steps]; ec = [err_steps(n) for n in steps]
    o3, oc = orders(e3), orders(ec)
    print('  steps   3 iterations   order    converged     order')
    for i, n in enumerate(steps):
        print(f"  {n:5d}   {e3[i]:.3e}   " + (f"{o3[i-1]:5.2f}" if i else '     ') + f"    {ec[i]:.3e}   " + (f"{oc[i-1]:5.2f}" if i else ''))
    OUT['iteration count'] = dict(steps=steps, err_fixed3=e3, err_converged=ec, order_fixed3=o3, order_converged=oc)
    good_c = min(oc[-2:]) >= 1.9
    print('PASS/FAIL converged iteration is second order in dt (orders of the two finest pairs >= 1.9) ->', 'PASS' if good_c else 'FAIL')
    ok &= good_c
    json.dump(OUT, open('../results/test_i_dust.json', 'w'), indent=1)
    print('\nTEST (i) OVERALL:', 'PASS' if ok else 'FAIL')
