import sys; sys.path.insert(0,'.'); sys.path.insert(0,'..')
import numpy as np
from scipy.linalg import eigh
from galerkin_fem import assemble
import eos_neutronized as EN
import initial_profile_neu as IPN
import eos as E

def count_nodes(xi):
    s = np.sign(xi[2:-2]); s = s[s!=0]
    return np.sum(np.abs(np.diff(s))>0)

print("Cross-check: a robustly STABLE configuration (rho_c=1e11, same as the nonlinear-hydro control test)\n")
rho_c = 1e11
sol = IPN.build_star(rho_c, rmax=3e6, n=8000)
r_raw, m_raw, rho_raw = sol.t, sol.y[0], np.maximum(sol.y[1], 1.0)
R_star = r_raw[-1]; Mstar = m_raw[-1]
print(f"M={Mstar/E.Msun:.4f} Msun, R={R_star/1e3:.1f} km (classically stable branch, dM/drho_c>0)\n")

N = 500
r = np.concatenate([[0.0], np.linspace(r_raw[1], R_star*0.9999, N-1)])
rho0 = np.interp(r, r_raw, rho_raw); rho0[0] = rho_raw[0]
P0 = EN.P_of_rho(rho0)
Gamma1_arr = np.where(rho0 >= EN._rho_tab[0], EN._dlogPdlogrho(np.log(np.maximum(rho0,EN._rho_tab[0]))), 5/3)
print(f"Gamma1 range: [{Gamma1_arr.min():.5f}, {Gamma1_arr.max():.5f}]")

M, D = assemble(r, rho0, P0, Gamma1_arr, G=E.G, include_self_gravity=True)
M = M[1:,1:]; D = D[1:,1:]
eigvals, eigvecs = eigh(M, D)
order = np.argsort(eigvals)

fund = None
print("\nLowest eigenvalues with node counts:")
for idx in order[:8]:
    xi = eigvecs[:,idx]
    nn = count_nodes(xi)
    print(f"  omega^2={eigvals[idx]:+.6e}  nodes={nn}")
    if fund is None and nn==0: fund = eigvals[idx]

print(f"\nFundamental mode: omega^2 = {fund:+.6e} s^-2")
print("STABLE" if fund>0 else "UNSTABLE -- would contradict the expected classical stability!")
if fund>0: print(f"Oscillation period = {2*np.pi/np.sqrt(fund):.4e} s")
