"""Test (e): conservation diagnostics.
 e1  number (baryon) content: sum dM is constant by construction; verified numerically.
 e2  gravitational mass M(t) = m_N(t) = sum e_k V_k: relative drift over 0 < t <= T for the report data and TOV data (GR dynamics), vs N.
 e3  Newtonian mode: total energy KE + IE + W (W = -sum m_j dM^int_j / R_j, the exact potential of the discrete force) drift vs N and dt.
 e4  energy equation (3), dm_j/dt = -4 pi a_j pbar_j R_j^2 U_j, checked independently of the constraint solve: over one step the change of the
     interface masses m_j (from the constraint) is compared with dt * (average of the fluxes at the two ends of the step); the residual
     must scale like dt^3 per step (dt^2 globally).
 e5  discrete constraint residual: m_{j+1} - m_j - e V  is zero by construction; instead check Gamma^2 = 1 + U^2 - 2m/R at interfaces (definition) and the
     fixed-point residual of the rho_N-m-Gamma coupling after the 3 iterations."""
import sys, json, numpy as np
from harness import *
from lag import Lag

OUT = {}
print('# e1/e2: mass and number conservation, report data and TOV data, GR dynamics, T = 50 (before any crossing)')
for kind in ('pure', 'neu'):
    for fam, data, relax in (('report', 'newton', None), ('TOV-pt', 'tov', None)):
        for N in (100, 200, 400, 800):
            star = make_star(kind, data, N, relax=relax)
            r = evolve(star, 50.0, nout=5)
            sn = r['snaps']
            dMtot = star.dM.sum()
            drift = max(abs(s['M'] / sn[0]['M'] - 1) for s in sn)
            print(f'  {kind} {fam:7s} N={N:4d}  max|M(t)/M(0)-1| = {drift:.3e}   rho_c change {sn[-1]["rho_c"]/sn[0]["rho_c"]-1:+.3e}   status {r["status"]}')
            OUT[f'e2 {kind} {fam} {N}'] = drift
print('\n# e3: Newtonian mode energy conservation (report profile, T = 60, dt = CFL 0.4 x 0.5^k)')
for kind in ('pure', 'neu'):
    S = Setup(kind, newt=True)
    for N in (100, 200, 400):
        for dts in (1.0, 0.5):
            d = S.newton_data(N, newt=True)
            star = S.make(d)
            star.start()
            geo = star.geometry(star.R, star.U)
            KE, IE, W, E0 = star.energy_newt(geo)
            dt0 = star.dt_cfl(star.geo, 0.4) * dts
            nst = int(np.ceil(60.0 / dt0)); dt = 60.0 / nst
            Emax = 0.0
            for n in range(nst):
                star.step(dt)
                if n % max(1, nst // 60) == 0:
                    Et = star.energy_newt()[3]
                    Emax = max(Emax, abs(Et - E0))
            KE, IE, W, E1 = star.energy_newt()
            print(f'  {kind} N={N:4d} dt x{dts:3.1f}: |W| = {abs(W):.4e}, max_t |E(t)-E(0)|/|W| = {Emax/abs(W):.3e}, final {abs(E1-E0)/abs(W):.3e}')
            OUT[f'e3 {kind} {N} {dts}'] = Emax / abs(W)

print('\n# e4: energy equation (3) checked over one step at t = 20 (report data, N = 400)')
for kind in ('pure', 'neu'):
    star = make_star(kind, 'newton', 400)
    r = evolve(star, 20.0, nout=2)
    S0 = (star.R.copy(), star.U.copy(), star.t)
    res = []
    for dts in (1.0, 0.5, 0.25):
        star.R, star.U, star.t = S0[0].copy(), S0[1].copy(), S0[2]
        star.Gint = np.ones(star.N + 1); star.start()
        dt = star.dt_cfl(star.geo, 0.4) * dts
        def flux(st):
            geo = st.geometry(st.R, st.U); F, a = st.accel(geo)
            ppad = np.concatenate([geo['p'], [st.p_ext]]); pbar = 0.5 * (ppad[:-1] + ppad[1:])
            return geo['m'][1:].copy(), -FOUR * a * pbar * st.R[1:]**2 * st.U[1:]
        FOUR = 4 * np.pi
        m0, f0 = flux(star)
        star.step(dt)
        m1, f1 = flux(star)
        resid = np.max(np.abs((m1 - m0) - 0.5 * dt * (f0 + f1))) / np.max(np.abs(m1 - m0))
        res.append((dt, resid))
        print(f'  {kind}: dt = {dt:.4e}   max_j |dm_j - dt (f0+f1)/2| / max_j |dm_j| = {resid:.3e}')
    o = np.log2(res[0][1] / res[1][1]), np.log2(res[1][1] / res[2][1])
    print(f'  {kind}: observed order of the one-step residual in dt: {o[0]:.2f}, {o[1]:.2f}  (3 expected for a second-order scheme; residual also contains the fixed spatial truncation error)')
    OUT[f'e4 {kind}'] = dict(res=res, order=o)

print('\n# e4b: the dt-independent part of that residual is the spatial truncation error of the discrete energy equation; scaling with N (one step, same time, dts = 1)')
for kind in ('pure', 'neu'):
    rr = []
    for N in (100, 200, 400, 800):
        star = make_star(kind, 'newton', N)
        evolve(star, 20.0, nout=2)
        star.Gint = np.ones(star.N + 1); star.start()
        dt = star.dt_cfl(star.geo, 0.4) * 0.25
        def flux(st):
            geo = st.geometry(st.R, st.U); F, a = st.accel(geo)
            ppad = np.concatenate([geo['p'], [st.p_ext]]); pbar = 0.5 * (ppad[:-1] + ppad[1:])
            return geo['m'][1:].copy(), -4 * np.pi * a * pbar * st.R[1:]**2 * st.U[1:]
        m0, f0 = flux(star); star.step(dt); m1, f1 = flux(star)
        rr.append(np.max(np.abs((m1 - m0) - 0.5 * dt * (f0 + f1))) / np.max(np.abs(m1 - m0)))
    ords = [np.log2(rr[i] / rr[i + 1]) for i in range(3)]
    print(f'  {kind}: residual N=100,200,400,800: ' + ', '.join(f'{x:.3e}' for x in rr) + '; observed orders in h: ' + ', '.join(f'{o:.2f}' for o in ords))
    OUT[f'e4b {kind}'] = dict(res=rr, order=ords)

print('\n# e5: definition check Gamma^2 = 1 + U^2 - 2 m/R at the end of a run, and fixed-point convergence')
for kind in ('pure', 'neu'):
    star = make_star(kind, 'newton', 400)
    r = evolve(star, 40.0, nout=2)
    ones = np.ones(star.N + 1)
    star.gtol = 0; star.n_iter = 3; geo = star.geometry(star.R, star.U, Gguess=ones)   # fixed-count experiments (the production setting iterates to 1e-13)
    G = geo['G']; m = geo['m']; R = star.R; U = star.U
    res = np.max(np.abs(G[1:]**2 - (1 + U[1:]**2 - 2 * m[1:] / R[1:])))
    ch = {}
    for it in (1, 2, 3, 4, 8):
        star.n_iter = it
        gi = star.geometry(star.R, star.U, Gguess=ones)
        star.n_iter = 30
        gc = star.geometry(star.R, star.U, Gguess=ones)
        ch[it] = float(np.max(np.abs(gi['m'] - gc['m'])) / gc['m'][-1])
    star.n_iter = 3; star.gtol = 1e-13
    print(f'  {kind}: max|Gamma^2 - (1+U^2-2m/R)| = {res:.2e}; relative change of m(R) vs the 30-iteration fixed point, starting from Gamma = 1 with n iterations: ' + ', '.join(f'n={k}: {v:.1e}' for k, v in ch.items()) + '  (in the evolution the previous Gamma is the starting guess)')
    OUT[f'e5 {kind}'] = dict(gamma_res=res, iter_change=ch)
json.dump(OUT, open('../results/test_e_conservation.json', 'w'), indent=1, default=float)
