"""Fig. (collapse-timescale) of Paper C: calibrated homologous model, delta_0 = 1e-6, integrated to 1.02*lambda_crit.
Also prints the fitted offset c in t(lambda=0.5) = [ln(1/delta_0) - c]/sigma from final_single.json (c = 0.63).
Writes ../fig8_collapse_timescale.png.  Runtime: ~1-3 min."""
import os, json
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import background as bg
from homologous_calibrated import CalibratedHomologousModel
HERE = os.path.dirname(os.path.abspath(__file__)); G, c = bg.G, bg.c
bg.mu_e = 2.0
m = CalibratedHomologousModel(5e13); f, kN, kF = m.calibrate_to_RR(); m.f_cal = f
lcrit = np.max(2*G*m.Mr_bg[1:]/(m.r[1:]*c**2))
ev = lambda t, y: y[0] - 1.02*lcrit; ev.terminal = True
s = solve_ivp(lambda t, y: [y[1], m.lambda_ddot(y[0], f)], (0, 6), [1 - 1e-6, 0], method='Radau', events=ev, max_step=0.02, rtol=1e-9, atol=1e-13, dense_output=True)
print("f_cal", f, "t_end", s.t[-1], "lambda_end", s.y[0, -1], "rho amplification", 1/s.y[0, -1]**3)
fs = json.load(open(os.path.join(HERE, 'final_single.json'))); sig = np.sqrt(fs['target'])
cs = [np.log(1/float(d0)) - t05*sig for d0, (t05, dt) in fs['traj'].items()]
print("offset c in t(0.5)=[ln(1/delta0)-c]/sigma:", cs)
t = np.linspace(0, s.t[-1], 4000); lam = s.sol(t)[0]
fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
ax[0].plot(t, lam, 'b'); ax[0].axhline(lcrit, c='r', ls='--', lw=1, label=r'$\lambda_{\rm crit}=\max_r 2GM_r/(rc^2)=%.1f\times10^{-3}$' % (lcrit*1e3))
ax[0].set_xlabel('t (s)'); ax[0].set_ylabel(r'$\lambda(t)$'); ax[0].set_title('Homologous scale factor'); ax[0].legend(fontsize=8)
ax[1].semilogy(t, lam**-3, 'g'); ax[1].set_xlabel('t (s)'); ax[1].set_ylabel(r'$\rho_c(t)/\rho_c(0)$'); ax[1].set_title('Central density amplification')
plt.tight_layout(); plt.savefig(os.path.join(HERE, '..', 'fig8_collapse_timescale.png'), dpi=150)
