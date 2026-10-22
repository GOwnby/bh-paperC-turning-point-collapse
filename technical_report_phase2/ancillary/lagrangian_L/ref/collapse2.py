import numpy as np, pickle
from scipy.interpolate import CubicSpline
import eos as E, initial_profile as IP

rho_c = 5e12
sol = IP.build_star(rho_c, rmax=2e6, n=8000)
r_raw, m_raw, rho_raw = sol.t, sol.y[0], np.maximum(sol.y[1],0)
r_surf = r_raw[-1]
print(f"Star: M={m_raw[-1]/E.Msun:.4f} Msun, R={r_surf/1e3:.1f} km")

N = 300
r_grid = np.linspace(0.01*r_surf, r_surf, N)   # start at 1% of R, not at the singular center
m_grid = np.interp(r_grid, r_raw, m_raw)
rho_grid_init = np.interp(r_grid, r_raw, rho_raw)

L = r_surf
r_lab = r_grid / L
m_units = (m_grid*E.G/E.c**2) / L
Gamma_shell = np.sqrt(np.maximum(1 - 2*m_units/r_lab, 1e-10))
print(f"compactness 2m/R at surface: {2*m_units[-1]/r_lab[-1]:.6e}")
print(f"rho_dimless at innermost grid point vs expected central value:")
Gscale = E.G*L**2/E.c**2
print(f"  computed from m-gradient at i=0 would be unreliable; using interpolated rho directly there instead")
print(f"  rho_c dimless expected = {Gscale*rho_c:.4e}")

pickle.dump(dict(r_lab=r_lab, m_units=m_units, Gamma_shell=Gamma_shell, L=L, rho_c=rho_c,
                  Mtot=m_raw[-1], rho_grid_init=rho_grid_init), open('initdata2.pkl','wb'))
print("saved initdata2.pkl")
