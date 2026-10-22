"""
Phase A1 (cont.): TOV integration using the Chandrasekhar EOS of background.py.
Standard TOV:
  dP/dr  = -G (rho + P/c^2)(M_r + 4*pi*r^3*P/c^2) / [r^2 (1 - 2*G*M_r/(r*c^2))]
  dM_r/dr = 4*pi*r^2*rho
with rho = rho(P) obtained by inverting P_of_rho (monotonic).
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import background as bg

G, c = bg.G, bg.c

def rho_of_P(P, rho_guess):
    """Invert P_of_rho via bracketed root find around rho_guess."""
    if P <= 0:
        return 0.0
    lo, hi = rho_guess*1e-8, rho_guess*10
    # expand bracket until sign change
    def g(rho):
        return bg.P_of_rho(rho) - P
    f_lo, f_hi = g(lo), g(hi)
    tries = 0
    while f_lo*f_hi > 0 and tries < 60:
        lo /= 2.0
        hi *= 2.0
        f_lo, f_hi = g(lo), g(hi)
        tries += 1
    return brentq(g, lo, hi, xtol=1e-4, rtol=1e-12)

def tov_rhs(r, y, rho_c):
    P, M = y
    if P <= 0:
        return [0.0, 0.0]
    rho = rho_of_P(P, rho_c)
    if r < 1e-6:
        dPdr = 0.0
        dMdr = 4*np.pi*r**2*bg.eps_over_c2_of_rho(rho)
    else:
        eps = bg.eps_over_c2_of_rho(rho)
        num = G*(eps + P/c**2)*(M + 4*np.pi*r**3*P/c**2)
        den = r**2*(1 - 2*G*M/(r*c**2))
        dPdr = -num/den
        dMdr = 4*np.pi*r**2*eps
    return [dPdr, dMdr]

def surface_event(r, y, rho_c):
    return y[0] - 1e-10*y0_P  # relative pressure threshold, set below
surface_event.terminal = True
surface_event.direction = -1

def solve_tov(rho_c, r_max=2e7, n_eval=12000):
    global y0_P
    P_c = bg.P_of_rho(rho_c)
    y0_P = P_c
    r0 = 10.0  # start a few meters out with a tiny-sphere analytic approx
    M0 = (4.0/3.0)*np.pi*r0**3*rho_c
    P0 = P_c - 2*np.pi*G*(rho_c+P_c/c**2)*(rho_c+3*P_c/c**2)*r0**2/3.0  # tiny correction, negligible
    sol = solve_ivp(tov_rhs, [r0, r_max], [P0, M0], args=(rho_c,),
                     method='RK45', rtol=1e-10, atol=[P_c*1e-14, 1e-6],
                     events=surface_event, dense_output=True, max_step=r_max/2000)
    if len(sol.t_events[0]) == 0:
        raise RuntimeError("surface not found within r_max")
    R = sol.t_events[0][0]
    M = sol.y_events[0][0][1]
    r_arr = np.linspace(r0, R, n_eval)
    P_arr = sol.sol(r_arr)[0]
    M_arr = sol.sol(r_arr)[1]
    P_arr = np.clip(P_arr, 0, None)
    rho_arr = np.array([rho_of_P(P, rho_c) for P in P_arr])
    return r_arr, rho_arr, P_arr, M_arr, R, M

if __name__ == "__main__":
    rho_c = 5e13  # kg/m^3, Paper C's "corrected configuration"
    r_arr, rho_arr, P_arr, M_arr, R, M = solve_tov(rho_c)
    Msun = bg.M_sun
    compactness = G*M/(R*c**2)
    Gamma1_arr = np.array([bg.Gamma1_of_rho(rho) if rho>1 else 5/3 for rho in rho_arr])
    print(f"rho_c = {rho_c:.3e} kg/m^3")
    print(f"R     = {R:.6e} m")
    print(f"M     = {M:.6e} kg = {M/Msun:.6f} Msun")
    print(f"compactness GM/(Rc^2) = {compactness:.6f}")
    print(f"Gamma1(center) = {Gamma1_arr[0]:.6f}")
    print(f"Gamma1(min)    = {Gamma1_arr.min():.6f}")
    print(f"Gamma1(surface)= {Gamma1_arr[-5]:.6f}")
    print()
    print("Paper C: M=1.4143 Msun, R=815.8 km, 2GM/(Rc^2)=0.00512 (printed value above is GM/(Rc^2)), Gamma1_min=1.33409 (at center)")
