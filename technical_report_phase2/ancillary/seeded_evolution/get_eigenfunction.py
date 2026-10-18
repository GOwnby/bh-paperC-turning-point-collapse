import sys; sys.path.insert(0,'..'); sys.path.insert(0,'../linear_perturbation')
import numpy as np
from scipy.linalg import eigh
from scipy.interpolate import interp1d
from galerkin_fem import assemble
import eos_neutronized as EN
import initial_profile_neu as IPN
import eos as E
import pickle

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
for idx in order[:5]:
    xi = eigvecs[:,idx]
    if count_nodes(xi) == 0:
        fund_idx, fund_val = idx, eigvals[idx]; break

xi_fund = eigvecs[:,fund_idx]
r_fund = r[1:]
# normalize so max|xi| = 1
xi_fund = xi_fund / np.max(np.abs(xi_fund))
# prepend xi(0)=0 (the point we removed from the eigenvalue problem by the regularity BC)
r_full = np.concatenate([[0.0], r_fund])
xi_full = np.concatenate([[0.0], xi_fund])

print(f"Fundamental mode: omega^2={fund_val:.6e} s^-2, growth rate={np.sqrt(-fund_val):.6e} s^-1")
print(f"xi(r) extracted on {len(r_full)} points, r from 0 to {r_full[-1]:.4e} m (R_star={R_star:.4e} m)")

pickle.dump(dict(r=r_full, xi=xi_full, omega2=fund_val, rho0_grid=rho0, r_grid=r,
                  R_star=R_star, rho_c=rho_c, m_raw=m_raw, r_raw=r_raw, rho_raw=rho_raw),
            open('eigenfunction_data.pkl','wb'))
print("saved eigenfunction_data.pkl")
