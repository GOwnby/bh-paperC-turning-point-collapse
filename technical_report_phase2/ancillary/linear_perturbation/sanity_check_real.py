import sys; sys.path.insert(0,'.'); sys.path.insert(0,'..')
import numpy as np
import eos_neutronized as EN
import initial_profile_neu as IPN
import eos as E

rho_c = 5e12
sol = IPN.build_star(rho_c, rmax=3e6, n=8000)
r_raw, m_raw, rho_raw = sol.t, sol.y[0], np.maximum(sol.y[1], 1.0)

N = 2000
r = np.linspace(r_raw[1], r_raw[-1]*0.999, N)
rho0 = np.interp(r, r_raw, rho_raw)
m = np.interp(r, r_raw, m_raw)
P0 = EN.P_of_rho(rho0)
Gamma1 = np.where(rho0 >= EN._rho_tab[0], EN._dlogPdlogrho(np.log(np.maximum(rho0,EN._rho_tab[0]))), 5/3)

xi = r.copy()
divxi = np.ones(N)*3
Pp = -Gamma1*P0*divxi - xi*np.gradient(P0,r)
rho_prime = -rho0*divxi - xi*np.gradient(rho0,r)
phip = -4*np.pi*E.G*rho0*xi   # exact l=0 self-gravity closure, independently verified earlier
dPhi0dr = E.G*m/r**2
rhs = np.gradient(Pp,r) + rho0*phip + rho_prime*dPhi0dr
omega2_pt = rhs/(rho0*xi)

for frac in [0.1,0.3,0.5,0.7,0.9]:
    i = int(frac*(N-1))
    print(f"r/R={frac:.1f}: local omega^2 = {omega2_pt[i]:+.4e} s^-2   (Galerkin fundamental: -18.68 s^-2)")

rq = np.trapezoid(r**2*xi*rhs, r) / np.trapezoid(r**2*rho0*xi**2, r)
print(f"\nIntegrated Rayleigh quotient (independent, non-matrix): omega^2 = {rq:+.4e} s^-2")
print(f"Galerkin fundamental (for comparison):                   omega^2 = -1.8677e+01 s^-2")
