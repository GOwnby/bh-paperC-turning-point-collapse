"""The report's own Galerkin (Newtonian, l=0) solver, called exactly as in ancillary/linear_perturbation/apply_to_real_star.py and
stable_crosscheck.py, generalised only in (EOS kind, rho_c, number of nodes N).  Used as the reference normal-mode frequency."""
import sys, os, io, contextlib
import numpy as np
from scipy.linalg import eigh
HERE = os.path.dirname(os.path.abspath(__file__)); REF = os.path.join(HERE, '..', 'ref')
sys.path.insert(0, REF)
with contextlib.redirect_stdout(io.StringIO()):
    import eos as E, eos_neutronized as EN, initial_profile as IP, initial_profile_neu as IPN
    from galerkin_fem import assemble


def galerkin_modes(kind, rho_c, N=500, nmodes=4, rmax=None):
    if rmax is None:
        rmax, nn = 3e6, 8000          # as in the report's scripts
    else:
        nn = int(rmax / 250)
    if kind == 'neu':
        sol = IPN.build_star(rho_c, rmax=rmax, n=nn)
        Pf = EN.P_of_rho
        G1 = lambda rho: np.where(rho >= EN._rho_tab[0], EN._dlogPdlogrho(np.log(np.maximum(rho, EN._rho_tab[0]))), 5 / 3)
        floor = 1.0
    else:
        sol = IP.build_star(rho_c, rmax=rmax, n=nn)
        Pf = IP.P_of_rho
        G1 = lambda rho: np.where(rho >= IP.rhos_tab[0], IP.dlogPdlogrho(np.log(np.maximum(rho, IP.rhos_tab[0]))), 5 / 3)
        floor = 1.0
    r_raw, m_raw, rho_raw = sol.t, sol.y[0], np.maximum(sol.y[1], floor)
    R_star = r_raw[-1]
    r = np.concatenate([[0.0], np.linspace(r_raw[1], R_star * 0.9999, N - 1)])
    rho0 = np.interp(r, r_raw, rho_raw); rho0[0] = rho_raw[0]
    P0 = Pf(rho0)
    Gamma1 = G1(rho0)
    M, D = assemble(r, rho0, P0, Gamma1, G=E.G, include_self_gravity=True)
    M = M[1:, 1:]; D = D[1:, 1:]
    w = eigh(M, D, eigvals_only=True, subset_by_index=[0, nmodes - 1])
    return np.sort(w)[:nmodes], R_star, m_raw[-1]


def galerkin_poly(Gam, K, rho_c, N=500, nmodes=4):
    """Newtonian polytrope P = K rho^Gam, same assembler (for the n ~ 3 / Gamma = 4/3 + eps check)."""
    from scipy.integrate import solve_ivp
    G = E.G
    def rhs(r, y):
        m, Hh = y
        rho = (max(Hh, 0) * (Gam - 1) / (Gam * K))**(1 / (Gam - 1))
        return [4 * np.pi * r**2 * rho, -G * m / r**2]
    Hc = Gam / (Gam - 1) * K * rho_c**(Gam - 1)
    r0 = 1.0
    ev = lambda r, y: y[1]; ev.terminal = True; ev.direction = -1
    sol = solve_ivp(rhs, [r0, 1e9], [4 / 3 * np.pi * r0**3 * rho_c, Hc], rtol=1e-11, atol=1e-30, events=ev, dense_output=True, method='DOP853')
    R = sol.t_events[0][0]
    r = np.concatenate([[0.0], np.linspace(r0, R * 0.9999, N - 1)])
    y = sol.sol(np.maximum(r, r0))
    rho0 = (np.maximum(y[1], 1e-300) * (Gam - 1) / (Gam * K))**(1 / (Gam - 1)); rho0[0] = rho_c
    P0 = K * rho0**Gam
    M, D = assemble(r, rho0, P0, np.full_like(r, Gam), G=G, include_self_gravity=True)
    w = eigh(M[1:, 1:], D[1:, 1:], eigvals_only=True, subset_by_index=[0, nmodes - 1])
    return np.sort(w)[:nmodes], R, sol.sol(R)[0]


if __name__ == '__main__':
    for N in (500, 1000, 2000):
        print(N, galerkin_modes('neu', 5e12, N), )
