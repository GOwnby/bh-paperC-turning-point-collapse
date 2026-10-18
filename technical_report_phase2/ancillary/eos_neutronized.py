import numpy as np
from scipy.interpolate import CubicSpline
import eos as E
import neutronization as NZ

# Precompute Y_e(rho) on a fast grid using the SEMF equilibrium (A=12), then
# build the modified pressure P(rho) = P_electron(n_e = Y_e(rho)*n_B).
rho_grid = np.geomspace(1e2, 5e14, 120)
nB_grid = rho_grid / NZ.u
Ye_grid = np.array([NZ.equilibrium_Z(NZ.A_rep, nB)/NZ.A_rep for nB in nB_grid])
Ye_of_lnrho = CubicSpline(np.log(rho_grid), Ye_grid)

def Ye_of_rho(rho):
    rho = np.clip(rho, rho_grid[0], rho_grid[-1])
    return Ye_of_lnrho(np.log(rho))

# electron pressure as a function of (rho, Ye): same relativistic electron gas,
# but now with the REDUCED electron density set by Ye(rho) instead of fixed mu_e=2 (Ye=0.5)
def x_of_ne(n_e):
    return (3*n_e/(8*np.pi))**(1/3) * E.h/(E.m_e*E.c)

def P_of_rho_neutronized(rho):
    rho = np.atleast_1d(np.asarray(rho, dtype=float))
    Ye = Ye_of_rho(rho)
    n_B = rho/NZ.u
    n_e = Ye*n_B
    x = x_of_ne(np.maximum(n_e, 1e-10))
    P = (np.pi*E.m_e**4*E.c**5)/(3*E.h**3) * E.f_cs(x)
    return P

# build a fast table + spline for dP/drho (needed by the hydro code), same interface as before
_rho_tab = np.geomspace(1e-2, 5e14, 3000)
_P_tab = P_of_rho_neutronized(_rho_tab)
_lr = np.log(np.maximum(_rho_tab,1e-300)); _lP = np.log(np.maximum(_P_tab,1e-300))
_logP_spline = CubicSpline(_lr, _lP)
_dlogPdlogrho = _logP_spline.derivative()
_K_NR = _P_tab[0]/_rho_tab[0]**(5/3)

def P_of_rho(rho):
    rho = np.asarray(rho, dtype=float)
    return np.where(rho >= _rho_tab[0],
                     np.exp(_logP_spline(np.log(np.maximum(rho,_rho_tab[0])))),
                     _K_NR*np.maximum(rho,0.0)**(5/3))

def dPdrho(rho):
    rho = np.asarray(rho, dtype=float)
    Gamma_local = np.where(rho >= _rho_tab[0],
                            _dlogPdlogrho(np.log(np.maximum(rho,_rho_tab[0]))), 5/3)
    P = P_of_rho(rho)
    return Gamma_local*P/np.maximum(rho,1e-300)

if __name__=="__main__":
    print("Comparing pure-Chandrasekhar vs neutronized pressure at several densities:")
    import initial_profile as IP
    for rho in [1e10, 1e11, 1e12, 3e12, 1e13, 3e13, 1e14]:
        P_pure = IP.P_of_rho(rho)
        P_neu = P_of_rho(rho)
        print(f"  rho={rho:.1e} kg/m^3   P_pure={P_pure:.4e} Pa   P_neutronized={P_neu:.4e} Pa   ratio={P_neu/P_pure:.4f}   Ye={Ye_of_rho(rho):.4f}")
