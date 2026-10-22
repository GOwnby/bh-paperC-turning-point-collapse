import numpy as np
from scipy.integrate import solve_ivp
import eos_neutronized as EN
import eos as E

def build_star(rho_c, rmax=2e7, n=6000):
    def rhs(r, y):
        m, rho = y
        if rho <= 1e-1: return [0.0, 0.0]
        dpdrho = EN.dPdrho(np.array([rho]))[0]
        dPdr = -E.G*m*rho/r**2 if r > 0 else 0.0
        drhodr = dPdr/dpdrho if dpdrho > 0 else 0.0
        dmdr = 4*np.pi*r**2*rho
        return [dmdr, drhodr]
    def surface_event(r, y):
        return y[1] - 1.0
    surface_event.terminal = True; surface_event.direction = -1
    r0 = 1.0
    sol = solve_ivp(rhs, [r0, rmax], [4/3*np.pi*r0**3*rho_c, rho_c], max_step=rmax/n,
                     events=surface_event, rtol=1e-8, atol=1e-3, method='RK45')
    return sol

if __name__ == "__main__":
    print("Mass-radius relation WITH neutronization softening:")
    for rc in np.geomspace(1e9, 5e12, 12):
        sol = build_star(rc)
        M = sol.y[0,-1]/E.Msun
        print(f"  rho_c={rc:.3e} kg/m^3   M={M:.4f} Msun   R={sol.t[-1]/1e3:.1f} km")
