import sys
sys.path.insert(0, 'lawe')
import numpy as np
from scipy.linalg import eigh
from galerkin_fem import assemble
import eos_neutronized as EN
import initial_profile_neu as IPN
import eos as E

def count_nodes(xi):
    s = np.sign(xi[2:-2]); s = s[s!=0]
    return np.sum(np.abs(np.diff(s))>0)

rho_c = 5e12
sol = IPN.build_star(rho_c, rmax=3e6, n=8000)
r_raw, m_raw, rho_raw = sol.t, sol.y[0], np.maximum(sol.y[1], 1.0)
R_star = r_raw[-1]
Mstar = m_raw[-1]
print(f"Target configuration: M={Mstar/E.Msun:.4f} Msun, R={R_star/1e3:.1f} km, rho_c={rho_c:.2e} kg/m^3\n")

N = 500
r = np.concatenate([[0.0], np.linspace(r_raw[1], R_star*0.9999, N-1)])
rho0 = np.interp(r, r_raw, rho_raw); rho0[0] = rho_raw[0]
P0 = EN.P_of_rho(rho0)

# LOCAL (real, varying) Gamma1 from the actual neutronized EOS
Gamma1_arr = np.where(rho0 >= EN._rho_tab[0],
                       EN._dlogPdlogrho(np.log(np.maximum(rho0, EN._rho_tab[0]))),
                       5/3)
print(f"Gamma1 range across the star: [{Gamma1_arr.min():.5f}, {Gamma1_arr.max():.5f}]")
print(f"Gamma1 at center: {Gamma1_arr[0]:.5f}, at r=R/2: {Gamma1_arr[N//2]:.5f}, near surface: {Gamma1_arr[-5]:.5f}")

M, D = assemble(r, rho0, P0, Gamma1_arr, G=E.G, include_self_gravity=True)
M = M[1:,1:]; D = D[1:,1:]
eigvals, eigvecs = eigh(M, D)
order = np.argsort(eigvals)

print("\nLowest eigenvalues with node counts:")
fund = None
for idx in order[:15]:
    xi = eigvecs[:,idx]
    nn = count_nodes(xi)
    print(f"  omega^2={eigvals[idx]:+.6e}  nodes={nn}")
    if fund is None and nn == 0:
        fund = eigvals[idx]
        fund_vec = xi

print(f"\nFundamental (n=0) mode: omega^2 = {fund:+.6e} s^-2")
if fund is not None and fund < 0:
    growth_rate = np.sqrt(-fund)
    print(f"UNSTABLE. e-folding timescale = 1/sqrt(|omega^2|) = {1/growth_rate:.4e} s")
elif fund is not None:
    print(f"STABLE. Oscillation period = {2*np.pi/np.sqrt(fund):.4e} s")
