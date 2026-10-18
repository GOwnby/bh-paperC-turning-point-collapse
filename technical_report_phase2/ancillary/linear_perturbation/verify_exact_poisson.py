"""
For purely radial (l=0) adiabatic perturbations, claim: phi'(r) = -4*pi*G*rho0(r)*xi(r)
satisfies the perturbed Poisson equation EXACTLY, with no need to solve a separate
coupled field equation. This follows from Delta_rho being exactly a total derivative
(a consequence of radial mass conservation), not from any additional approximation.
Verify this symbolically before relying on it.
"""
import sympy as sp

r = sp.symbols('r', positive=True)
rho0 = sp.Function('rho0')(r)
xi = sp.Function('xi')(r)
G = sp.symbols('G', positive=True)

divxi = sp.diff(xi,r) + 2*xi/r
rho_prime = -rho0*divxi - xi*sp.diff(rho0,r)   # Eulerian density perturbation, as derived/used throughout

phi_prime_claimed = -4*sp.pi*G*rho0*xi
# Poisson LHS using the claimed phi':  (1/r^2) d(r^2 * phi')/dr
poisson_lhs = sp.simplify(sp.diff(r**2*phi_prime_claimed, r)/r**2)
poisson_rhs = sp.simplify(4*sp.pi*G*rho_prime)

residual = sp.simplify(poisson_lhs - poisson_rhs)
print("Poisson equation residual using phi' = -4 pi G rho0 xi  (should be exactly 0):")
print(residual)
