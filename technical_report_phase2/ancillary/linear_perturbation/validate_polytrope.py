import numpy as np
from scipy.integrate import solve_ivp
from lawe_solver import solve_modes

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

print("Validating against the classical Pekeris (Gamma1=4/3) marginal-stability theorem")
print("for the FUNDAMENTAL radial mode of a self-gravitating polytrope.\n")
K = 1e10  # arbitrary polytropic constant, units chosen for a reasonably sized star
for Gamma in [1.25, 1.30, 4/3, 1.40, 1.50, 5/3]:
    sol = build_polytrope(K, Gamma, rho_c=1e3, rmax=5e8)
    r_raw, m_raw, rho_raw = sol.t, sol.y[0], np.maximum(sol.y[1], 1e-8)
    N = 400
    r = np.linspace(r_raw[1], r_raw[-1]*0.999, N)
    rho0 = np.interp(r, r_raw, rho_raw)
    m = np.interp(r, r_raw, m_raw)
    P0 = K*rho0**Gamma
    Gamma1 = np.full(N, Gamma)
    eigvals, eigvecs = solve_modes(r, rho0, P0, m, Gamma1, G, include_self_gravity=True, n_modes=6)
    # fundamental mode = smallest |omega^2| among the low-order, smooth (non-numerical-garbage) modes;
    # take the eigenvalue closest to zero as a proxy for the fundamental
    idx = np.argsort(np.abs(eigvals))[0]
    print(f"Gamma1={Gamma:.4f} (3*Gamma1-4={3*Gamma-4:+.4f}):  fundamental omega^2 ~ {eigvals[idx]:+.4e} 1/s^2   "
          f"{'UNSTABLE' if eigvals[idx]<0 else 'stable'}")
