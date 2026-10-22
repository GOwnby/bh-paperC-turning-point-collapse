"""Equilibrium stellar structure for the barotropic EOS of eosmodel.py (dimensionless units G=c=1, length L).

tov(eos, e_c)       exact general-relativistic hydrostatic solution (Tolman-Oppenheimer-Volkoff) in the enthalpy
                    variable H = ln h = int dp/(e+p), for which the surface is the exact zero H = 0:
                        dm/dr = 4 pi r^2 e,   dH/dr = -(m + 4 pi r^3 p) / (r^2 (1 - 2m/r)),
                        dM_b/dr = 4 pi r^2 rho_N / sqrt(1 - 2m/r)           (rest-mass/number content)
                    (Tolman 1939; Oppenheimer & Volkoff 1939).  Lapse: a h = const  (a = e^{nu/2}).
newtonian_profile   the report's own Newtonian profile (initial_profile*.build_star, RK45, rtol 1e-8), sampled by
                    linear interpolation.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
import eosmodel as M


def central_state(eos, e_c_SI):
    e_c = eos.Ge * e_c_SI
    # invert e(rho_N) at the centre by monotone interpolation of the table
    rN = np.exp(np.interp(np.log(e_c), np.log(eos.e_tab), np.log(eos.rN_tab)))
    for _ in range(3):  # polish with Newton on thermo (smooth)
        e, p, g1 = eos.thermo(np.array([rN]))
        rN = rN * (e_c / e[0])**(1.0 / (1.0 + p[0] / e[0]))   # d ln e/d ln rN = (e+p)/e
    e, p, g1 = eos.thermo(np.array([rN]))
    return float(rN), float(e[0]), float(p[0])


def tov(eos, e_c_SI, rtol=1e-12):
    rN_c, e_c, p_c = central_state(eos, e_c_SI)
    H_c = np.log((e_c + p_c) / rN_c)

    def rhs(r, y):
        m, H, Mb = y
        rN, e, p = eos.state_from_H(H)
        f = 1.0 - 2.0 * m / r
        return [4 * np.pi * r**2 * float(e), -(m + 4 * np.pi * r**3 * float(p)) / (r**2 * f), 4 * np.pi * r**2 * float(rN) / np.sqrt(f)]

    def surf(r, y): return y[1]
    surf.terminal = True; surf.direction = -1
    r0 = 1e-6
    y0 = [4 * np.pi / 3 * e_c * r0**3, H_c - 2 * np.pi / 3 * (e_c + 3 * p_c) * r0**2, 4 * np.pi / 3 * rN_c * r0**3]
    sol = solve_ivp(rhs, [r0, 5.0], y0, method='DOP853', rtol=rtol, atol=1e-22, events=surf, dense_output=True)
    R = sol.t_events[0][0]
    yR = sol.sol(R)
    out = dict(R=R, M=yR[0], Mb=yR[2], sol=sol, e_c=e_c, p_c=p_c, rN_c=rN_c, H_c=H_c)
    return out


def tov_at(st, eos, r):
    """m, H, Mb, rho_N, e, p at radii r (array, 0<=r<=R)"""
    r = np.asarray(r, float)
    rr = np.maximum(r, 1e-6)
    y = st['sol'].sol(np.minimum(rr, st['R']))
    m, H, Mb = y
    # center: use series
    c = r < 1e-6
    m = np.where(c, 0.0, m); Mb = np.where(c, 0.0, Mb); H = np.where(c, st['H_c'], H)
    H = np.where(r >= st['R'], 0.0, H)
    m = np.where(r >= st['R'], st['M'], m); Mb = np.where(r >= st['R'], st['Mb'], Mb)
    rN, e, p = eos.state_from_H(np.maximum(H, 1e-300))
    return m, H, Mb, rN, e, p


def newtonian_profile(kind, rho_c=5e12):
    """The report's profile: returns L (=R_star, m), M (kg), and dimensionless callable m_d(r_d)"""
    L, sol = M.kind_L(kind, rho_c)
    r_raw, m_raw, rho_raw = sol.t, sol.y[0], np.maximum(sol.y[1], 0)
    Gm = M.G_SI / M.C_SI**2 / L
    # cubic-spline interpolation of the report's RK45 solution (it is sampled every ~250 m; linear interpolation
    # would inject O(1e-3) noise into m(r) at small r)
    m_spl = CubicSpline(np.concatenate([[0.0], r_raw]), np.concatenate([[0.0], m_raw]))
    def m_of_r(r_d, linear=False):
        if linear:
            return np.interp(np.asarray(r_d) * L, r_raw, m_raw) * Gm
        return m_spl(np.asarray(r_d) * L) * Gm
    def e_of_r(r_d):   # dimensionless energy density from the Newtonian rho
        return np.interp(np.asarray(r_d) * L, r_raw, rho_raw) * M.G_SI * L**2 / M.C_SI**2
    return dict(L=L, M_kg=m_raw[-1], m_of_r=m_of_r, e_of_r=e_of_r, Rd=1.0, Md=m_raw[-1] * Gm)


if __name__ == '__main__':
    for kind in ('pure', 'neu'):
        Lr, sol = M.kind_L(kind)
        eos = M.EOS(kind, Lr)
        st = tov(eos, 5e12)
        print(kind, f"Newtonian (report): R={Lr/1e3:.1f} km M={sol.y[0,-1]/M.MSUN:.4f} Msun")
        print(kind, f"TOV              : R={st['R']*Lr/1e3:.1f} km M={st['M']*Lr*M.C_SI**2/M.G_SI/M.MSUN:.4f} Msun  Mb={st['Mb']*Lr*M.C_SI**2/M.G_SI/M.MSUN:.4f}  2M/R={2*st['M']/st['R']:.3e}  e_c={st['e_c']:.4e}")
