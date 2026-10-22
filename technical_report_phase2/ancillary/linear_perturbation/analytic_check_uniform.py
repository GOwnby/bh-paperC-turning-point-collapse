"""First analytic test case of the weak form (Paper C, Sec. 3.2; companion of analytic_check2.py, whose density rho0 = rho_c (1 - (r/R)^2) vanishes at
the surface so that its boundary term is zero).  Here the star has UNIFORM density, rho0 = rho_c = const, i.e. a finite surface density: m(r) = (4 pi/3) rho_c r^3,
P0(r) = (2 pi/3) G rho_c^2 (R^2 - r^2), P0(R) = 0 but P0'(R) = -(4 pi/3) G rho_c^2 R != 0.  Fully closed form, no numerical differentiation.
Checks, for the trial functions xi = r and xi = r + k r^3/R^2 (symbolic k) and symbolic Gamma1 (constant):
   S_exact := int_0^R r^2 xi { (P')' + rho' g - 4 pi G rho0^2 xi } dr        (strong form, including the exact l=0 self-gravity -4 pi G rho0^2 xi)
   B_nobdry := int_0^R [A xi^2 + 2 C xi xi' + D xi'^2] dr                    (weak form without surface term)
   B_bdry   := B_nobdry - R^2 P0'(R) xi(R)^2                                  (weak form with the surface term of Eq. (boundary-term))
Expected: S_exact - B_bdry = 0 identically, and S_exact - B_nobdry = -R^2 P0'(R) xi(R)^2 != 0 (the first implementation error of the paper).
Run: python3 analytic_check_uniform.py        (needs sympy; exit status 1 on failure)
"""
import sys
import sympy as sp

r, R, rho_c, G, Gam, k = sp.symbols('r R rho_c G Gamma1 k', positive=True)
rho0 = rho_c
m = sp.Rational(4, 3)*sp.pi*rho_c*r**3
P0 = sp.Rational(2, 3)*sp.pi*G*rho_c**2*(R**2 - r**2)
assert sp.simplify(sp.diff(P0, r) + G*m*rho0/r**2) == 0       # hydrostatic equilibrium
ok = True
for name, xi in (('xi = r', r), ('xi = r + k r^3/R^2', r + k*r**3/R**2)):
    div = sp.diff(xi, r) + 2*xi/r
    Pp = -Gam*P0*div - xi*sp.diff(P0, r)
    rho_p = -rho0*div - xi*sp.diff(rho0, r)
    g = G*m/r**2
    S_exact = sp.simplify(sp.integrate(sp.expand(r**2*xi*(sp.diff(Pp, r) + rho_p*g - 4*sp.pi*G*rho0**2*xi)), (r, 0, R)))
    A = 4*Gam*P0 + sp.diff(P0, r)*r**2*sp.diff(rho0, r)/rho0 + 4*sp.diff(P0, r)*r - 4*sp.pi*G*rho0**2*r**2
    C = r*(2*Gam*P0 + sp.diff(P0, r)*r); D = Gam*P0*r**2
    xip = sp.diff(xi, r)
    B_nb = sp.simplify(sp.integrate(sp.expand(A*xi*xi + 2*C*xi*xip + D*xip*xip), (r, 0, R)))
    bnd = sp.simplify(-R**2*sp.diff(P0, r).subs(r, R)*xi.subs(r, R)**2)
    res_with = sp.simplify(S_exact - (B_nb + bnd)); res_without = sp.simplify(S_exact - B_nb)
    print(f"trial {name}:  S_exact = {S_exact}")
    print(f"   boundary term = {bnd}")
    print(f"   S_exact - B(with boundary term)    = {res_with}")
    print(f"   S_exact - B(without boundary term) = {res_without}   (should equal the boundary term {bnd})")
    good = (res_with == 0) and (sp.simplify(res_without - bnd) == 0) and (bnd != 0)
    print("   [%s]" % ("PASS" if good else "FAIL"))
    ok &= good
sys.exit(0 if ok else 1)
