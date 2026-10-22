"""Fig. (friedman-correspondence) of Paper C: Friedman's Newtonian canonical energy E and the classical virial pulsation constant
I w^2 = (3 Gamma1 - 4)|W| evaluated on self-consistent uniform-Gamma1 Lane-Emden polytropes (Gamma1 = 1 + 1/n).
Lengths are in Lane-Emden units (xi); |W| and I are quoted without the common factor 4 pi, which cancels in w^2.
The ratio w^2_virial/(E/I) = 27 int P r^2 dr / int P r^4 dr carries units of length^-2, so its numerical value
(5.80 at Gamma1=1.25 ... 9.19 at 5/3) refers to these units.  Each curve in the left panel is normalized to its own maximum magnitude over the Gamma1 values shown. Writes ../fig_friedman_correspondence.png.  Runtime: seconds."""
import os
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, quad
HERE = os.path.dirname(os.path.abspath(__file__))

def lane_emden(n):
    f = lambda x, y: [y[1], -max(y[0], 0)**n - 2*y[1]/x]
    ev = lambda x, y: y[0]; ev.terminal = True
    x0 = 1e-6
    s = solve_ivp(f, [x0, 60], [1 - x0**2/6, -x0/3], events=ev, rtol=1e-12, atol=1e-14, dense_output=True)
    return s, s.t_events[0][0]

def quantities(G1):
    n = 1/(G1 - 1); s, x1 = lane_emden(n)
    th = lambda x: max(s.sol(x)[0], 0.0)
    P = lambda x: th(x)**(n + 1)
    I = quad(lambda x: th(x)**n*x**4, 1e-6, x1, limit=200)[0]
    W = 3*quad(lambda x: P(x)*x**2, 1e-6, x1, limit=200)[0]
    Eint = quad(lambda x: P(x)*x**4, 1e-6, x1, limit=200)[0]
    om2 = (3*G1 - 4)*W/I; EI = (3*G1 - 4)*Eint/9/I
    return om2, EI

gam = np.array([1.25, 1.28, 1.30, 1.32, 4/3, 1.34, 1.36, 1.40, 1.50, 5/3])
vals = np.array([quantities(g) for g in gam]); om2, EI = vals[:, 0], vals[:, 1]

for g, a, b in zip(gam, om2, EI): print(f"Gamma1={g:.4f}  omega2_virial={a:+.5e}  E/I={b:+.5e}  ratio={a/b if abs(b)>0 else float('nan'):.4f}")
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].plot(gam, om2/np.max(np.abs(om2)), 'o-', label=r'virial $\omega^2=(3\Gamma_1-4)|W|/I$ (normalized)')
ax[0].plot(gam, EI/np.max(np.abs(EI)), 's--', label=r'Friedman canonical energy $E/I$ (normalized)')
ax[0].axvline(4/3, ls=':', c='gray', lw=0.8); ax[0].axhline(0, c='gray', lw=0.8)
ax[0].set_xlabel(r'$\Gamma_1$ (uniform, self-consistent polytrope)'); ax[0].set_ylabel('normalized value'); ax[0].legend(fontsize=8, loc='upper left')
ax[0].set_title('Both vanish at the Pekeris point')
mk = np.abs(gam - 4/3) > 1e-9
ax[1].plot(gam[mk], (om2/EI)[mk], 'd-', c='#2a9d2a'); ax[1].axvline(4/3, ls=':', c='gray', lw=0.8)
ax[1].set_xlabel(r'$\Gamma_1$'); ax[1].set_ylabel(r'ratio: virial $\omega^2$ / Friedman $E/I$  (Lane-Emden units)')
ax[1].set_title('Same sign throughout; smoothly varying ratio')
plt.tight_layout(); plt.savefig(os.path.join(HERE, '..', 'fig_friedman_correspondence.png'), dpi=200)
