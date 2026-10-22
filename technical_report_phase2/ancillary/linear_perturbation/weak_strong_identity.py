"""Numerical check of the weak form of the Newtonian radial pulsation problem against the strong form, on self-consistent polytropes
(Paper C, Sec. 3.2: 'checking the weak form against the strong form as an exact algebraic identity ... holds to 5-6 significant figures for every
Gamma1 tested'; the symbolic derivation of the same identity is code/symbolic/sympy_weak_form.py).

For a trial displacement xi(r) (the identity holds for ANY trial function, so a mismatch cannot be blamed on a poor trial function):
   strong:  S[xi] = int_0^R r^2 xi { (P')' + rho' g - 4 pi G rho0^2 xi } dr,      P' = -Gamma1 P0 (xi' + 2 xi/r) - xi P0',  rho' = -rho0 (xi'+2xi/r) - xi rho0'
   weak:    B[xi]  = int_0^R [A xi^2 + 2 C xi xi' + D xi'^2] dr  - R^2 P0'(R) xi(R)^2          (A, C, D of Eq. (ACD); boundary term of Eq. (boundary-term))
The two must agree (S = B) when the equilibrium reaches its true surface P0(R) = 0.  All background derivatives are taken analytically from the
equilibrium ODEs (no numerical differentiation).  Also shown: (i) B without the boundary term, (ii) the SAME check on a 'truncated' polytrope whose integration
is stopped before the true surface (error #2 of the paper: builder stopping at a fixed outer radius), which fails by a large resolution-independent factor.
Usage (from this directory):  python3 weak_strong_identity.py         (~5 s)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.integrate import solve_ivp
from validate_polytrope import build_polytrope, G

K, rho_c = 1e10, 1e3

def equilibrium(Gam, rmax=None):
    """dense equilibrium (m, rho) of the polytrope P = K rho^Gam, integrated to the true surface (rho = 1e-8 rho_c, as in build_polytrope)
    or, if rmax is given, stopped at rmax (a 'truncated' star)"""
    def rhs(r, y):
        m, rho = y
        dPdr = -G*m*rho/r**2
        return [4*np.pi*r**2*rho, dPdr/(K*Gam*rho**(Gam - 1))]
    surf = lambda r, y: y[1] - 1e-8*rho_c
    surf.terminal = True; surf.direction = -1
    r0 = 1.0
    s = solve_ivp(rhs, [r0, 1e11 if rmax is None else rmax], [4/3*np.pi*r0**3*rho_c, rho_c], events=None if rmax else surf,
                  rtol=1e-12, atol=[1e-6, 1e-14], method='DOP853', dense_output=True)
    return s

def identity(Gam, trial_c, rmax=None):
    s = equilibrium(Gam, rmax); R = s.t[-1]
    r = np.linspace(1e-3*R, R, 80001)
    m, rho = s.sol(r)
    P0 = K*rho**Gam
    P0p = -G*m*rho/r**2
    rho_p = P0p/(K*Gam*rho**(Gam - 1))
    P0pp = -G*(4*np.pi*r**2*rho*rho + m*rho_p)/r**2 + 2*G*m*rho/r**3
    g = G*m/r**2
    xi = r + trial_c*r**3/R**2; xip = 1 + 3*trial_c*r**2/R**2; div = 3 + 5*trial_c*r**2/R**2; divp = 10*trial_c*r/R**2
    Pe_p = -Gam*P0p*div - Gam*P0*divp - xip*P0p - xi*P0pp          # d/dr P'
    rho_e = -rho*div - xi*rho_p
    strong = np.trapezoid(r**2*xi*(Pe_p + rho_e*g - 4*np.pi*G*rho**2*xi), r)
    A = 4*Gam*P0 + r**2*P0p*rho_p/rho + 4*r*P0p - 4*np.pi*G*rho**2*r**2
    C = r*(2*Gam*P0 + r*P0p); D = Gam*P0*r**2
    weak_nb = np.trapezoid(A*xi**2 + 2*C*xi*xip + D*xip**2, r)
    weak = weak_nb - R**2*P0p[-1]*xi[-1]**2
    scale = np.trapezoid(np.abs(A*xi**2) + np.abs(2*C*xi*xip) + np.abs(D*xip**2), r)
    return strong, weak, weak_nb, scale, R

print('Gamma1  trial  R[m]       strong S            weak B (with bdry)   (S-B)/scale   (S-B_nobdry)/scale')
for Gam in (1.25, 4/3, 1.40, 1.50, 5/3):
    for tc in (0.0, 0.3):
        S, B, Bnb, sc, R = identity(Gam, tc)
        print('%.4f  %.1f  %.4e  %+.8e  %+.8e   %.2e      %.2e' % (Gam, tc, R, S, B, (S - B)/sc, (S - Bnb)/sc))
# the builder of the repository stops at the same surface
for Gam in (1.25, 5/3):
    sol = build_polytrope(K, Gam, rho_c=rho_c, rmax=5e8)
    print('validate_polytrope.build_polytrope R = %.6e m  (this script: %.6e m)' % (sol.t[-1], equilibrium(Gam).t[-1]))
print('\nTruncated star (integration stopped at 0.7 R): the identity fails, resolution-independently (error #2 of the paper)')
for Gam in (1.25, 4/3, 5/3):
    R = equilibrium(Gam).t[-1]
    S, B, Bnb, sc, Rt = identity(Gam, 0.0, rmax=0.7*R)
    print('Gamma1 = %.4f: R_trunc/R = %.2f  (S-B)/scale = %.2e' % (Gam, Rt/R, (S - B)/sc))
