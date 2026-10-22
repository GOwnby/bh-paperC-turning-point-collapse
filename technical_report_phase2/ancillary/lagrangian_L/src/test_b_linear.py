"""Test (b): linear oscillation / instability.
 b1  spectrum of the linearised semi-discrete GR scheme (Jacobian eigenvalues) vs N, Richardson-extrapolated, against
       (i) the independent GR Sturm-Liouville solver (Chandrasekhar 1964 / Gondek et al. 1997; linear_gr.py),
       (ii) in Newtonian mode, the report's Galerkin solver (ancillary/linear_perturbation/galerkin_fem.py, called as the report does).
 b2  Newtonian polytrope P = K rho^(4/3 + eps): scheme (Newtonian mode) vs report Galerkin, and ->0 as eps->0 (Pekeris).
 b3  time-domain: small-amplitude oscillation of the stable rho_c=1e11 star (pure and neutronised): frequency of the evolution vs eigenvalue.
 b4  time-domain: seeded growth of the unstable mode of the rho_c=5e12 neutronised star vs the linear e-folding rate.
omega^2 in s^-2 (physical units, time = proper time of the surface shell)."""
import sys, json, numpy as np
from scipy.optimize import curve_fit
from harness import *
import linear_gr as LG
import report_galerkin as RG
import eosmodel as M, structure as S
from lag import Lag, relax_equilibrium, linear_modes

OUT = {}
SEC = sys.argv[1] if len(sys.argv) > 1 else 'b1b2b3b4'
def rich(vals, p=2.0):
    """Richardson extrapolation from the last two values of a halving sequence"""
    return vals[-1] + (vals[-1] - vals[-2]) / (2**p - 1)

if 'b1' in SEC:
    print('# b1: GR linear spectrum of the scheme vs independent Sturm-Liouville solver')
    configs = [('neu', 5e12), ('pure', 5e12), ('neu', 1e11), ('pure', 1e11)]
    for kind, rc in configs:
        Sx = Setup(kind, rho_c=rc)
        s = (C_SI / Sx.L)**2
        res = {}
        for N in (50, 100, 200, 400):
            star = make_star(kind, 'tov', N, rho_c=rc, relax='db')
            om2, vec, imag = linear_modes(star, 4)
            res[N] = om2 * s
        arr = np.array([res[N] for N in (50, 100, 200, 400)])
        lim = [rich(arr[:, i]) for i in range(4)]
        wsl = {n: LG.modes(Sx.eos, Sx.tov, n=n, nmodes=4)[0] * s for n in (500, 1000, 2000)}
        sl = np.array([wsl[n] for n in (500, 1000, 2000)])
        sllim = [rich(sl[:, i]) for i in range(4)]
        print(f'\n{kind} rho_c={rc:.0e}  (L/c = {Sx.L/C_SI*1e3:.3f} ms)   omega^2 [s^-2]')
        for N in (50, 100, 200, 400):
            print(f'   scheme N={N:4d}:', np.array2string(res[N], precision=5, floatmode='fixed'))
        print('   scheme Richardson (p=2):', np.array2string(np.array(lim), precision=5, floatmode='fixed'))
        for n in (500, 1000, 2000):
            print(f'   SL-FEM n={n:4d}:', np.array2string(wsl[n], precision=5, floatmode='fixed'))
        print('   SL-FEM Richardson      :', np.array2string(np.array(sllim), precision=5, floatmode='fixed'))
        rd = [abs(a - b) / abs(b) for a, b in zip(lim, sllim)]
        print('   relative difference scheme-limit vs SL-limit:', np.array2string(np.array(rd), precision=2))
        ok = max(rd[:3]) < 1e-3
        print('   PASS' if ok else '   FAIL', '(criterion: first three modes agree to < 1e-3)')
        # observed order of the scheme's fundamental eigenvalue
        d = np.diff(arr[:, 0]); pord = np.log2(d[:-1] / d[1:])
        print('   observed order of fundamental eigenvalue (N=50..400):', np.array2string(pord, precision=2))
        OUT[f'b1 {kind} {rc:g}'] = dict(scheme={N: res[N].tolist() for N in res}, scheme_lim=lim, sl={n: wsl[n].tolist() for n in wsl}, sl_lim=sllim,
                                         rel_diff=rd, order=pord.tolist(), pass_=bool(ok), M=Sx.tov['M'], R=Sx.tov['R'], L=Sx.L)
        # Newtonian-mode scheme vs report Galerkin
        Sn = Setup(kind, rho_c=rc, newt=True)
        resn = {}
        for N in (50, 100, 200, 400):
            star = make_star(kind, 'newton', N, newt=True, rho_c=rc, relax='db')
            om2, vec, imag = linear_modes(star, 4)
            resn[N] = om2 * s
        arrn = np.array([resn[N] for N in (50, 100, 200, 400)])
        limn = [rich(arrn[:, i]) for i in range(4)]
        gal = {n: RG.galerkin_modes(kind, rc, n, rmax=(None if rc >= 1e12 else 2.5e7))[0] for n in (500, 1000, 2000)}
        gl = np.array([gal[n] for n in (500, 1000, 2000)])
        gallim = [rich(gl[:, i]) for i in range(4)]
        print(f'   --- Newtonian limit (scheme in Newtonian mode on the report profile) vs report Galerkin')
        for N in (50, 100, 200, 400): print(f'   newt scheme N={N:4d}:', np.array2string(resn[N], precision=5, floatmode='fixed'))
        print('   newt scheme Richardson :', np.array2string(np.array(limn), precision=5, floatmode='fixed'))
        for n in (500, 1000, 2000): print(f'   report Galerkin N={n:4d}:', np.array2string(gal[n], precision=5, floatmode='fixed'))
        print('   report Galerkin Richardson:', np.array2string(np.array(gallim), precision=5, floatmode='fixed'))
        rdn = [abs(a - b) / abs(b) for a, b in zip(limn, gallim)]
        print('   relative difference:', np.array2string(np.array(rdn), precision=2), ' PASS' if max(rdn[:3]) < 2e-3 else ' FAIL', '(criterion < 2e-3, first three modes)')
        OUT[f'b1N {kind} {rc:g}'] = dict(scheme={N: resn[N].tolist() for N in resn}, scheme_lim=limn, galerkin={n: gal[n].tolist() for n in gal}, galerkin_lim=gallim, rel_diff=rdn)

if 'b2' in SEC:
    # ---------------- b2 polytrope
    print('\n# b2: Newtonian polytrope P = K rho^Gamma, Gamma = 4/3 + eps')
    K0 = 1.0e7   # SI, arbitrary scale
    rows = []
    for eps in (0.01, 0.02, 0.04, 0.08):
        Gam = 4 / 3 + eps
        rho_c = 1e10
        # choose K so that R ~ 1e7 m: scale-free, keep K0 * (adjust) -- compute R first
        K = K0 * (rho_c**(4 / 3 - Gam))
        wg, R, Mt = RG.galerkin_poly(Gam, K, rho_c, N=2000)
        # scheme in Newtonian mode on the same polytrope
        from scipy.integrate import solve_ivp
        G = M.G_SI
        def rhs(r, y):
            m, Hh = y
            rho = (max(Hh, 0) * (Gam - 1) / (Gam * K))**(1 / (Gam - 1))
            return [4 * np.pi * r**2 * rho, -G * m / r**2]
        Hc = Gam / (Gam - 1) * K * rho_c**(Gam - 1)
        ev = lambda r, y: y[1]; ev.terminal = True; ev.direction = -1
        sol = solve_ivp(rhs, [1.0, 1e10], [4 / 3 * np.pi * rho_c, Hc], rtol=1e-12, atol=1e-30, events=ev, dense_output=True, method='DOP853')
        Rs = sol.t_events[0][0]
        e = np.geomspace(1e-3, 1e16, 3001); p = K * e**Gam
        tab = dict(e=e, p=p, rhoN=e, cum=np.zeros_like(e))
        eos = M.EOS('poly', Rs, newt=True, table=tab)
        res = []
        for N in (100, 200, 400):
            r = np.linspace(0, 1, N + 1)
            m = sol.sol(np.maximum(r * Rs, 1.0))[0] * G / M.C_SI**2 / Rs; m[0] = 0
            star = Lag(eos, np.diff(m), r, newt=True)
            relax_equilibrium(star)
            om2, vec, im = linear_modes(star, 2)
            res.append(om2[0] * (M.C_SI / Rs)**2)
        lim = rich(res)
        rows.append((eps, wg[0], lim, res))
        print(f'  eps={eps:5.3f}  Galerkin(N=2000) omega0^2 = {wg[0]:.6e} s^-2   Newtonian scheme (N=100,200,400) = {res[0]:.6e} {res[1]:.6e} {res[2]:.6e}, Richardson {lim:.6e}   rel.diff {abs(lim-wg[0])/abs(wg[0]):.2e}   omega0^2 / (eps G rho_c) = {wg[0]/(eps*G*rho_c):.4f}')
    OUT['b2'] = [dict(eps=a, galerkin=b, scheme_lim=c, scheme=d) for a, b, c, d in rows]
    slope = [r[1] / (r[0] * M.G_SI * 1e10) for r in rows]
    print('  Pekeris limit: omega0^2 -> 0 linearly in eps (ratio above ~ const as eps -> 0):', np.array2string(np.array(slope), precision=4))
    print('  b2 PASS' if max(abs(r[2] - r[1]) / abs(r[1]) for r in rows) < 2e-3 else '  b2 FAIL', '(criterion: scheme-limit vs Galerkin(N=2000) < 2e-3)')

# ---------------- b3 time-domain frequency
def seed_run(kind, rc, N, mode, eps, T_in_periods=None, T=None, dt_scale=1.0, nout=400, unstable=False):
    star = make_star(kind, 'tov', N, rho_c=rc, relax='db')
    om2, vec, im = linear_modes(star, 3)
    v = vec[:, mode]; v = v / np.max(np.abs(v)); v *= np.sign(v[np.argmax(np.abs(v))])
    Rsurf = star.R[-1]
    star.R[1:] += eps * Rsurf * v
    F, a = star.force_R(star.R[1:] - eps * Rsurf * v)   # a at equilibrium
    if unstable:
        sig = np.sqrt(-om2[mode])
        star.U[1:] = sig * eps * Rsurf * v / a
    star.start()
    return star, om2[mode], a

if 'b3' in SEC:
    print('\n# b3: stable star rho_c = 1e11 kg/m^3: oscillation frequency of the evolution vs linear eigenvalue (small amplitude, two periods, fit with 2nd harmonic)')
    for kind in ('pure', 'neu'):
        Sx = Setup(kind, rho_c=1e11); s = (C_SI / Sx.L)**2
        rows = []
        for N, dts in ((100, 1.0), (100, 0.5), (200, 1.0), (200, 0.5), (400, 1.0)):
            star, w2, a = seed_run(kind, 1e11, N, 0, 1e-4, dt_scale=dts)
            per = 2 * np.pi / np.sqrt(w2)
            T = 2.0 * per
            star.start()
            st0 = dict(R=star.R.copy())
            dt0 = star.dt_cfl(star.geo, 0.4) * dts
            nst = 400 * int(np.ceil(T / dt0 / 400)); dt = T / nst
            ts, rc_ = [], []
            for n in range(nst):
                star.step(dt)
                if (n + 1) % (nst // 400) == 0: ts.append(star.t); rc_.append(star.geo['e'][0])
            ts = np.array(ts); y = np.array(rc_) / star.geo['e'][0] * 0 + np.array(rc_)
            y = y / y[0] - 1 if False else (np.array(rc_) - np.mean(rc_)) / np.mean(rc_)
            f = lambda t, A, w, ph, c, B, ps: A * np.cos(w * t + ph) + c + B * np.cos(2 * w * t + ps)
            p0 = [y.std() * 1.4, np.sqrt(w2), 0.0, 0.0, 0.0, 0.0]
            popt, _ = curve_fit(f, ts, y, p0=p0, maxfev=50000, xtol=1e-14, ftol=1e-14)
            w2m = popt[1]**2
            rows.append(dict(N=N, dt_scale=dts, w2_eig=w2 * s, w2_evol=w2m * s, rel=abs(w2m - w2) / w2, amp=abs(popt[0])))
            print(f'  {kind} N={N:4d} dt/dt_CFL*{dts:3.1f}: eigenvalue omega^2 = {w2*s:.5f} s^-2, evolution fit = {w2m*s:.5f} s^-2   rel.diff = {abs(w2m-w2)/w2:.2e}  (density amplitude {abs(popt[0]):.2e}, nsteps {nst})')
        OUT[f'b3 {kind}'] = rows
        # reference to SL
        wsl = LG.modes(Sx.eos, Sx.tov, n=2000, nmodes=1)[0][0] * s
        wsl1 = LG.modes(Sx.eos, Sx.tov, n=1000, nmodes=1)[0][0] * s
        print(f'  independent Sturm-Liouville (n=1000, 2000): {wsl1:.5f}, {wsl:.5f} -> Richardson {rich([wsl1, wsl]):.5f} s^-2;  period {2*np.pi/np.sqrt(rich([wsl1,wsl])):.4f} s')
        last = [r for r in rows if r['N'] == 400][0]
        print('  b3', kind, 'PASS' if abs(last['w2_evol'] - rich([wsl1, wsl])) / rich([wsl1, wsl]) < 3e-3 else 'FAIL', '(criterion: evolution frequency^2 at N=400 within 3e-3 of the SL limit)')

if 'b4' in SEC:
    # ---------------- b4 growth rate
    print('\n# b4: seeded instability of rho_c = 5e12 neutronised star (TOV, discretely balanced): growth of the central density vs e-folding rate')
    Sx = Setup('neu'); s = (C_SI / Sx.L)**2
    wsl = LG.modes(Sx.eos, Sx.tov, n=2000, nmodes=1)[0][0]; wsl1 = LG.modes(Sx.eos, Sx.tov, n=1000, nmodes=1)[0][0]
    sig_ref = np.sqrt(-rich([wsl1, wsl]))               # dimensionless (units c/L)
    print(f'  reference growth rate (SL limit): omega^2 = {rich([wsl1,wsl])*s:.4f} s^-2 -> sigma = {sig_ref*C_SI/Sx.L:.5f} s^-1 = {sig_ref:.6f} c/L, e-folding {1/sig_ref:.3f} L/c = {1/sig_ref*Sx.L/C_SI:.5f} s')
    rows = []
    for N, dts, eps in ((100, 1.0, 1e-4), (200, 1.0, 1e-4), (400, 1.0, 1e-4), (400, 0.5, 1e-4), (200, 1.0, 1e-5), (200, 1.0, 1e-3)):
        star, w2, a = seed_run('neu', 5e12, N, 0, eps, unstable=True)
        sig = np.sqrt(-w2)
        T = 60.0
        star.start()
        dt0 = star.dt_cfl(star.geo, 0.4) * dts
        nst = 300 * int(np.ceil(T / dt0 / 300)); dt = T / nst
        ts, y = [], []
        e0 = star.geo['e'][0]
        # unperturbed central density from the equilibrium state
        st_eq = make_star('neu', 'tov', N, relax='db'); st_eq.start(); e_eq = st_eq.geo['e'][0]
        for n in range(nst):
            star.step(dt)
            if (n + 1) % (nst // 300) == 0: ts.append(star.t); y.append(star.geo['e'][0] / e_eq - 1)
        ts = np.array(ts); y = np.array(y)
        sel = (ts > 10) & (ts < 50) & (np.abs(y) > 0)
        slope = np.polyfit(ts[sel], np.log(np.abs(y[sel])), 1)[0]
        # amplitude check: y(t) vs y0 cosh-like pure growing
        rows.append(dict(N=N, dt_scale=dts, eps=eps, sigma_eig=sig, sigma_fit=slope, rel=abs(slope - sig) / sig, ymax=float(np.max(np.abs(y)))))
        print(f'  N={N:4d} dtx{dts:3.1f} eps={eps:.0e}: discrete sigma = {sig:.6f}, fit of ln|drho_c| over t in (10,50) = {slope:.6f}  rel.diff {abs(slope-sig)/sig:.2e}  (max |drho_c/rho_c| = {np.max(np.abs(y)):.2e})  vs SL sigma {sig_ref:.6f}: {abs(slope-sig_ref)/sig_ref:.2e}')
    OUT['b4'] = dict(rows=rows, sigma_ref=sig_ref, omega2_ref=rich([wsl1, wsl]) * s)
    good = [r for r in rows if r['N'] == 400 and r['dt_scale'] == 1.0][0]
    print('  b4', 'PASS' if abs(good['sigma_fit'] - sig_ref) / sig_ref < 5e-3 else 'FAIL', '(criterion: fitted growth rate at N=400 within 5e-3 of the SL growth rate)')
json.dump(OUT, open('../results/test_b_linear.json' if SEC == 'b1b2b3b4' else '../results/test_b_%s.json' % SEC, 'w'), indent=1, default=float)
