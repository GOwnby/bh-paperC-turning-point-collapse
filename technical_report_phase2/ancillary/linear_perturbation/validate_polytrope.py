"""Polytrope equilibrium builder with the surface-finding fix (bug #2 of Paper C, Sec. 3.2).
The older strong-form check that used to follow here (module lawe_solver, not part of the release) has been removed;
the Pekeris validation is validate_galerkin_sg.py (and code/fig1_pekeris.py for Fig. 1 of Paper C)."""
import numpy as np
from scipy.integrate import solve_ivp

G = 6.674e-11

def build_polytrope(K, Gamma, rho_c, rmax=1e8, n=3000, rho_floor_frac=1e-8):
    def rhs(r, y):
        m, rho = y
        if rho <= 1e-6: return [0.0, 0.0]
        P = K*rho**Gamma
        dPdrho = K*Gamma*rho**(Gamma-1)
        dPdr = -G*m*rho/r**2 if r>0 else 0.0
        drhodr = dPdr/dPdrho
        dmdr = 4*np.pi*r**2*rho
        return [dmdr, drhodr]
    surf_rho = max(rho_c*rho_floor_frac, 1e-6)
    def surf(r,y): return y[1]-surf_rho
    surf.terminal=True; surf.direction=-1
    r0=1.0
    # auto-expand rmax until the surface event actually triggers (don't silently truncate the star)
    for attempt in range(8):
        sol = solve_ivp(rhs, [r0,rmax], [4/3*np.pi*r0**3*rho_c, rho_c], max_step=rmax/n,
                         events=surf, rtol=1e-10, atol=1e-14, method='RK45')
        if len(sol.t_events[0]) > 0:
            return sol
        rmax *= 4   # surface not found yet -- integration hit rmax first; expand and retry
    raise RuntimeError(f"Surface not found after expanding rmax to {rmax:.2e}; rho at cutoff={sol.y[1][-1]:.3e}")
