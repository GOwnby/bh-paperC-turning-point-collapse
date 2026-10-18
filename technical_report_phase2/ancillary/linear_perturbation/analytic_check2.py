"""
Second analytic test case, with rho0'  != 0 (the constant-density test had
rho0'=0 identically, so it never exercised that term in the A coefficient).
Choose an arbitrary smooth rho0(r) = rho_c*(1-(r/R)^2), then build m(r) and
P0(r) self-consistently from the defining ODEs via symbolic integration --
no np.gradient anywhere, no numerical ODE solve, fully closed-form.
"""
import sympy as sp

r, R, rho_c, G, Gamma1 = sp.symbols('r R rho_c G Gamma1', positive=True)
rp = sp.symbols('rp', positive=True)

rho0 = rho_c*(1 - (r/R)**2)
rho0_rp = rho0.subs(r, rp)
m = sp.integrate(4*sp.pi*rp**2*rho0_rp, (rp, 0, r))
m = sp.simplify(m)
print("m(r) =", m)

# P0(r) = integral_r^R  G*m(r')*rho0(r')/r'^2  dr'   (so that P0(R)=0)
m_rp = m.subs(r, rp)
integrand_P0 = G*m_rp*rho0_rp/rp**2
P0 = sp.integrate(integrand_P0, (rp, r, R))
P0 = sp.simplify(P0)
print("\nP0(r) =", P0)
print("\nP0(0) [central pressure, sanity check -- should be positive]:", sp.simplify(P0.subs(r,0)))

xi = r
divxi = sp.diff(xi,r) + 2*xi/r
Pp = -Gamma1*P0*divxi - xi*sp.diff(P0,r)
rho_prime = -rho0*divxi - xi*sp.diff(rho0,r)
dPhi0dr = G*m/r**2

strong_rhs = sp.diff(Pp, r) + rho_prime*dPhi0dr
LHS_exact = sp.integrate(sp.expand(r**2*xi*strong_rhs), (r, 0, R))
LHS_exact = sp.simplify(LHS_exact)
print("\nLHS_exact =", LHS_exact)

A = 4*Gamma1*P0 + sp.diff(P0,r)*r**2*sp.diff(rho0,r)/rho0 + 4*sp.diff(P0,r)*r
C = r*(2*Gamma1*P0 + sp.diff(P0,r)*r)
Dc = Gamma1*P0*r**2
xip = sp.diff(xi, r)
weak_integrand = sp.expand(A*xi*xi + 2*C*xi*xip + Dc*xip*xip)
weak_no_boundary = sp.integrate(weak_integrand, (r, 0, R))
weak_no_boundary = sp.simplify(weak_no_boundary)
print("\nweak_no_boundary =", weak_no_boundary)

boundary_term = sp.simplify(-R**2*sp.diff(P0,r).subs(r,R)*xi.subs(r,R)**2)
print("\nboundary_term =", boundary_term)

residual = sp.simplify(LHS_exact - (weak_no_boundary + boundary_term))
print("\nLHS_exact - (weak_no_boundary + boundary_term) =", residual)
