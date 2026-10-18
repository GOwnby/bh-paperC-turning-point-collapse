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

rho_c = 5e12
sol = IPN.build_star(rho_c, rmax=3e6, n=8000)
r_raw, m_raw, rho_raw = sol.t, sol.y[0], np.maximum(sol.y[1], 1.0)
R_star = r_raw[-1]

N = 500
r = np.concatenate([[0.0], np.linspace(r_raw[1], R_star*0.9999, N-1)])
rho0 = np.interp(r, r_raw, rho_raw); rho0[0] = rho_raw[0]
P0 = EN.P_of_rho(rho0)
Gamma1_arr = np.where(rho0 >= EN._rho_tab[0], EN._dlogPdlogrho(np.log(np.maximum(rho0,EN._rho_tab[0]))), 5/3)

M, D = assemble(r, rho0, P0, Gamma1_arr, G=E.G, include_self_gravity=True)
M = M[1:,1:]; D = D[1:,1:]
eigvals, eigvecs = eigh(M, D)
order = np.argsort(eigvals)
for idx in order[:3]:
    xi = eigvecs[:,idx]
    if count_nodes(xi)==0:
        fund_idx = idx; break

xi = eigvecs[:,fund_idx]
xi = xi/np.max(np.abs(xi))
r_plot = r[1:]

print(f"Fundamental eigenfunction shape (xi/xi_max vs r/R), omega^2={eigvals[fund_idx]:.4f}:\n")
for frac in [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99,1.0]:
    i = min(int(frac*(N-2)), N-2)
    rho_here = np.interp(r_plot[i], r_raw, rho_raw)
    print(f"  r/R={r_plot[i]/R_star:.3f}  xi/xi_max={xi[i]:+.4f}  xi/r (relative compression)={xi[i]/(r_plot[i]/R_star):+.4f}  rho(r)/rho_c={rho_here/rho_c:.3e}")

# compare to the simple homologous trial xi=r for context
print(f"\nFor comparison, homologous (xi=r) would give xi/xi_max = r/R exactly, and xi/r = 1 everywhere.")
print(f"Deviation from homologous indicates how much the TRUE mode differs from uniform compression.")
