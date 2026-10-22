import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
import eos as E

# Precompute a fast, monotonic table of the Chandrasekhar EOS once (x from 1e-6 to 2000
# covers rho from far below white-dwarf densities to far past nuclear density).
xs = np.geomspace(1e-3, 2000, 4000)
rhos_tab = E.rho_of_x(xs)
Ps_tab = E.P_of_x(xs)
lr = np.log(rhos_tab); lP = np.log(Ps_tab)
logP_of_logrho = CubicSpline(lr, lP)
dlogPdlogrho = logP_of_logrho.derivative()

# Non-relativistic power-law slope (Gamma=5/3) anchored at the table's lower edge,
# used below the table so that P -> 0 smoothly and continuously as rho -> 0,
# instead of an artificial floor (which would otherwise create a spurious pressure
# discontinuity at any surface/vacuum boundary).
K_NR = Ps_tab[0] / rhos_tab[0]**(5/3)

def P_of_rho(rho):
    rho = np.asarray(rho, dtype=float)
    out = np.where(rho >= rhos_tab[0],
                   np.exp(logP_of_logrho(np.log(np.maximum(rho, rhos_tab[0])))),
                   K_NR*np.maximum(rho, 0.0)**(5/3))
    return out

def dPdrho(rho):
    rho = np.asarray(rho, dtype=float)
    Gamma_local = np.where(rho >= rhos_tab[0],
                            dlogPdlogrho(np.log(np.maximum(rho, rhos_tab[0]))),
                            5/3)
    P = P_of_rho(rho)
    return Gamma_local * P / np.maximum(rho, 1e-300)

def build_star(rho_c, rmax=2e7, n=6000):
    def rhs(r, y):
        m, rho = y
        if rho <= rhos_tab[0]*1.01:
            return [0.0, 0.0]
        dpdrho = dPdrho(rho)
        dPdr = -E.G*m*rho/r**2 if r > 0 else 0.0
        drhodr = dPdr/dpdrho if dpdrho > 0 else 0.0
        dmdr = 4*np.pi*r**2*rho
        return [dmdr, drhodr]
    def surface_event(r, y):
        return y[1] - rhos_tab[0]*2
    surface_event.terminal = True; surface_event.direction = -1
    r0 = 1.0
    sol = solve_ivp(rhs, [r0, rmax], [4/3*np.pi*r0**3*rho_c, rho_c], max_step=rmax/n,
                     events=surface_event, rtol=1e-8, atol=1e-3, method='RK45')
    return sol

if __name__ == "__main__":
    print("Scanning central density to find the mass-radius relation and the maximum mass (Chandrasekhar limit)...")
    rhocs = np.geomspace(1e8, 5e12, 16)
    masses = []
    for rc in rhocs:
        sol = build_star(rc)
        M = sol.y[0, -1]/E.Msun
        R = sol.t[-1]
        masses.append(M)
        print(f"  rho_c={rc:.3e} kg/m^3   M={M:.4f} Msun   R={R/1e3:.1f} km")
    masses = np.array(masses)
    i_max = np.argmax(masses)
    print(f"\nMaximum mass found in this scan: {masses[i_max]:.4f} Msun at rho_c={rhocs[i_max]:.3e} kg/m^3")
    print(f"(Chandrasekhar's classic analytic result for mu_e=2: M_Ch = 5.836/mu_e^2 = {5.836/2**2:.4f} Msun)")
