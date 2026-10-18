import numpy as np, sys
from scipy.integrate import solve_ivp
from twomode import TwoMode
import background as bg
G, c = bg.G, bg.c
T = TwoMode(); m = T.m
K = T.jac(); w, v = np.linalg.eig(K); i = np.argmax(w.real); sig2 = w[i].real; vec = v[:, i].real; vec = vec/vec[0]
print("sigma^2", sig2, "sigma", np.sqrt(sig2), "1/sigma", 1/np.sqrt(sig2), "mu/lam eigvec", vec[1])
lam_h = lambda lam: 1-2*G*m.M/(lam*m.R*c**2)
def rhs(t, y):
    q = y[:2]
    if not T.valid(q): return [0, 0, 0, 0]
    a = T.acc(q); return [y[2], y[3], a[0], a[1]]
def ev_h(t, y): return lam_h(y[0])-0.02
ev_h.terminal = True; ev_h.direction = -1
def ev_valid(t, y): return min(y[0]+y[1], y[0]-2*y[1], y[0]+y[1]*(1-3)) - 1e-3*y[0]*0 - 1e-4
ev_valid.terminal = True; ev_valid.direction = -1
def ev_small(t, y): return y[0]-1e-3
ev_small.terminal = True; ev_small.direction = -1
out = {}
for d0 in (1e-3, 1e-6):
    y0 = np.concatenate([np.array([1.0, 0.0]) - d0*vec, [0, 0]])
    s = solve_ivp(rhs, (0, 6), y0, method='Radau', events=[ev_h, ev_valid, ev_small], max_step=0.02, rtol=1e-9, atol=1e-13, dense_output=True)
    lam, mu = s.y[0], s.y[1]
    th = s.t[np.argmax(lam < 0.5)] if (lam < 0.5).any() else np.nan
    t_end = s.t[-1]
    lamend = lam[-1]
    # plunge time: from lam=0.5 to end
    print(f"d0={d0:g}: status={s.status} t_half={th:.4f} t_end={t_end:.4f} tau_plunge={t_end-th:.4f} lam_end={lamend:.4e} mu/lam_end={mu[-1]/lam[-1]:.4f} events={[len(e) for e in s.t_events]}")
    for tc in (th-0.2, th, t_end-0.005):
        yy = s.sol(tc); print(f"   t={tc:.4f} lam={yy[0]:.4e} mu/lam={yy[1]/yy[0]:+.4f} rho_c/rho0={1/(yy[0]+yy[1])**3:.3e}")
