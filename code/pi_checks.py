"""Numerical checks on the Painleve-I fold passage of Paper C, Sec. 5.6 (complements pi_tritronquee.py, which gives s_* = 3.41167).
Usage: python3 pi_checks.py      (~30 s)
 (1) Literature cross-check: Boutroux tritronquee data y(0) = -0.1875543083404949, y'(0) = 0.3049055602612289 of y'' = 6y^2 + x
     (Joshi & Kitaev 2001) integrated forward to the first real pole x_*; s_* = 6^(1/5) x_* (paper: x_* = 2.384169, s_* = 3.41167);
     backward integration shows the solution follows y = -sqrt(-x/6) (tritronquee) for x -> -oo.
 (2) Universal solution u'' = u^2 + s from the slow-manifold start (u = -sqrt(S) - 1/(8 S^2), u' = 1/(2 sqrt S), S = 100): the value of s at which
     u reaches u_h = y/alpha for y = 1e-2 and Mdot = 1e-7 (alpha = 9.2e-6, u_h = 1086; paper: 3.337 against 3.337 of the rescaled model), and the
     normal-form time to lambda = 0.5 (y = 7), paper: 367 s against 366 s of the full model.
 (3) Oscillation about the slow manifold: an offset A in u at s = -100 moves the pole to 3.30 (A=+0.1), 3.46 (A=-0.1), 2.56 (A=1);
     the oscillation amplitude decays as (-s)^(-1/8) (WKB: frequency 2^(1/2) (-s)^(1/4), amplitude ~ frequency^(-1/2)).
"""
import numpy as np
from scipy.integrate import solve_ivp
f = lambda s, u: [u[1], u[0]**2 + s]

# (1) literature data
B = 6.0**-0.2
y0, yp0 = -0.1875543083404949, 0.3049055602612289
g = lambda x, Y: [Y[1], 6*Y[0]**2 + x]
ev = lambda x, Y: Y[0] - 1e8; ev.terminal = True
sol = solve_ivp(g, (0, 6), [y0, yp0], method='DOP853', events=ev, rtol=1e-13, atol=1e-14, max_step=0.005)
xs = sol.t[-1] + sol.y[0, -1]**-0.5
print('Boutroux tritronquee data: first real pole x_* = %.6f  (paper 2.384169)   s_* = 6^(1/5) x_* = %.5f  (paper 3.41167)' % (xs, 6.0**0.2*xs), flush=True)
back = solve_ivp(g, (0, -40), [y0, yp0], method='DOP853', rtol=1e-12, atol=1e-13, dense_output=True)
for xx in (-10, -20, -40):
    Y = back.sol(xx)[0]; print('   backward: x=%4d  y=%.6f  -sqrt(-x/6)=%.6f  (difference %.1e)' % (xx, Y, -np.sqrt(-xx/6), Y + np.sqrt(-xx/6)))

# (2) universal solution
S = 100.0; x = np.sqrt(S)
sol = solve_ivp(f, (-S, 3.5), [-x - 1/(8*S**2), 1/(2*x)], method='DOP853', rtol=1e-12, atol=1e-13, max_step=0.01, dense_output=True)
kap = 18.733; a = 0.004027; yr = 3.15576e7; Md = 1e-7/yr
be = (4*a/(kap**2*Md))**0.2; al = 2/(kap*be**2); uh = 0.01/al
ss = np.linspace(2.5, 3.41, 910001); uu = sol.sol(ss)[0]
i = np.argmax(uu > uh); print('Mdot=1e-7: beta = %.2f s, alpha = %.3e, u_h = %.1f ; universal solution reaches y = 1e-2 at s = %.4f  (paper 3.337)' % (be, al, uh, ss[i]))
u7 = 7.0/al; i = np.argmax(uu > u7)
print('normal form reaches y = 7 (lambda = 0.5) at s = %.4f -> t = %.1f s  (paper: 366 s full model, 367 s = 3.412 beta pole)  ;  3.41167 beta = %.1f s' % (ss[i], ss[i]*be, 3.41167*be))

# (3) offsets about the slow manifold
for A in (0.1, -0.1, 1.0):
    ev = lambda s, u: u[0] - 1e7; ev.terminal = True
    sl = solve_ivp(f, (-S, 60), [-x - 1/(8*S**2) + A, 1/(2*x)], method='DOP853', events=ev, rtol=1e-11, atol=1e-12, max_step=0.05)
    print('offset A = %+.1f at s=-100: pole s_* = %.3f   (paper: +0.1 -> 3.30, -0.1 -> 3.46, +1 -> 2.56)' % (A, sl.t[-1] + np.sqrt(6/sl.y[0, -1])))
    if A == 0.1:
        # decay of the oscillation amplitude about the slow manifold, measured as the envelope of u + sqrt(-s) + 1/(8 s^2)
        sl2 = solve_ivp(f, (-S, -2), [-x - 1/(8*S**2) + A, 1/(2*x)], method='DOP853', rtol=1e-12, atol=1e-13, max_step=0.01, dense_output=True)
        tt = np.linspace(-100, -3, 970001); dev = sl2.sol(tt)[0] + np.sqrt(-tt) + 1/(8*tt**2)
        # envelope: maxima of |dev| in windows
        wins = [(-100, -85), (-60, -50), (-35, -28), (-20, -16), (-12, -9)]
        env = []
        for lo, hi in wins:
            m = (tt >= lo) & (tt <= hi); env.append((0.5*(lo + hi), np.max(np.abs(dev[m]))))
        env = np.array(env); sl_ = np.polyfit(np.log(-env[:, 0]), np.log(env[:, 1]), 1)[0]
        print('   oscillation envelope vs -s:', ' '.join('(%.1f, %.4f)' % tuple(e) for e in env), ' fitted exponent %.3f (WKB: -1/8 = -0.125)' % sl_)
