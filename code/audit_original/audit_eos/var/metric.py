"""
Phase A2 part 1: metric functions a(r), b(r) on the TOV background.
  e^{-2b} = 1 - 2*G*M_r/(r*c^2)
  da/dr   = -dP/dr / (rho*c^2 + P)      (standard TOV hydrostatic-equilibrium form)
with a(R) fixed by matching the exterior Schwarzschild solution,
  a(R) = 0.5 * ln(1 - 2*G*M/(R*c^2)).
"""
import numpy as np
import background as bg
from tov_solve import solve_tov, tov_rhs

G, c = bg.G, bg.c

def build_metric(rho_c):
    r, rho, P, Mr, R, M = solve_tov(rho_c)
    n = len(r)

    # b(r) directly algebraic
    b = -0.5*np.log(1 - 2*G*Mr/(r*c**2))

    # dP/dr from the RHS function directly (consistent with the TOV integration)
    dPdr = np.array([tov_rhs(ri, [Pi, Mi], rho_c)[0] for ri, Pi, Mi in zip(r, P, Mr)])

    eps = bg.eps_over_c2_of_rho(rho)
    integrand_a = -dPdr/(eps*c**2 + P)
    # integrate cumulatively from center (a(0)=0 arbitrary) then shift
    a_rel = np.concatenate([[0.0], np.cumsum(0.5*(integrand_a[1:]+integrand_a[:-1])*np.diff(r))])
    a_R_target = 0.5*np.log(1 - 2*G*M/(R*c**2))
    a = a_rel - a_rel[-1] + a_R_target

    return dict(eps=eps, r=r, rho=rho, P=P, Mr=Mr, R=R, M=M, a=a, b=b, dPdr=dPdr)

if __name__ == "__main__":
    bgd = build_metric(5e13)
    print("Surface check: a(R) =", bgd['a'][-1], " vs exterior value",
          0.5*np.log(1-2*G*bgd['M']/(bgd['R']*c**2)))
    print("a(0) =", bgd['a'][0], " b(0) =", bgd['b'][0], " (should be ~0 at center for weak field, small)")
    print("a(R)-a(0) =", bgd['a'][-1]-bgd['a'][0])
    print("compactness 2GM/(Rc^2) =", 2*G*bgd['M']/(bgd['R']*c**2))
