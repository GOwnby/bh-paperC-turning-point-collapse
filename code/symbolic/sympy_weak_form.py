"""Symbolic verification of the weak form of the radial Newtonian pulsation problem of Paper C, Sec. 3.1-3.2.

Claims verified (each prints PASS/FAIL; exit status 1 if any fails):
 (1) For purely radial perturbations, with xi the Lagrangian displacement, the Eulerian perturbations are
       rho' = -rho0 (xi' + 2 xi/r) - xi rho0',   P' = -Gamma1 P0 (xi' + 2 xi/r) - xi P0',   phi' = -4 pi G rho0 xi
     (the last one is the exact solution of the perturbed Poisson equation, see verify_exact_poisson.py).
     The linearised momentum equation is   omega^2 rho0 xi = (P')' + rho0 phi'' ... written here as
       omega^2 rho0 xi = d/dr P' + rho' g - 4 pi G rho0^2 xi,     g = G M/r^2 = -P0'/rho0.
 (2) Multiplying by eta r^2 and integrating by parts once, the integrand of the bilinear form is, using hydrostatic
     equilibrium (P0' = -G M rho0/r^2, M' = 4 pi r^2 rho0) and NOTHING else,
       W = A xi eta + C (xi eta' + xi' eta) + D xi' eta',
       A = 4 Gamma1 P0 + r^2 P0' rho0'/rho0 + 4 r P0' - 4 pi G rho0^2 r^2,   C = r (2 Gamma1 P0 + r P0'),   D = Gamma1 P0 r^2
     (Eq. (ACD) of the paper), which is symmetric under xi <-> eta.
 (3) The boundary term of the integration by parts is  [r^2 eta P']_R = -R^2 P0'(R) xi(R) eta(R)  when Delta P(R) = 0
     (Eq. (boundary-term) of the paper); omitting it (taking P'(R)=0) is the first of the two implementation errors.
Run:  python3 sympy_weak_form.py
"""
import sys
import sympy as sp
from _zero import is_zero, report

r, G, R = sp.symbols('r G R', positive=True)
rho0 = sp.Function('rho0')(r); M = sp.Function('M')(r); P0 = sp.Function('P0')(r); Gam = sp.Function('Gamma1')(r)
xi = sp.Function('xi')(r); eta = sp.Function('eta')(r)
ok = True

div = sp.diff(xi, r) + 2*xi/r
P_E = -Gam*P0*div - xi*sp.diff(P0, r)                 # Eulerian pressure perturbation
rho_E = -rho0*div - xi*sp.diff(rho0, r)               # Eulerian density perturbation
g = G*M/r**2
S = sp.diff(P_E, r) + rho_E*g - 4*sp.pi*G*rho0**2*xi   # omega^2 rho0 xi = S   (with phi' = -4 pi G rho0 xi)

def hydro(e):
    """impose hydrostatic equilibrium and M' = 4 pi r^2 rho0"""
    Pp = -G*M*rho0/r**2
    Ppp = sp.diff(Pp, r)
    e = e.subs(sp.Derivative(P0, (r, 2)), Ppp).subs(sp.Derivative(P0, r), Pp)
    return e.subs(sp.Derivative(M, r), 4*sp.pi*r**2*rho0)

# (2) integrand after one integration by parts:  r^2 eta S = d/dr(r^2 eta P') + W
W_raw = -sp.diff(r**2*eta, r)*P_E + r**2*eta*(rho_E*g - 4*sp.pi*G*rho0**2*xi)
ok &= report("(2a) r^2 eta S - d/dr(r^2 eta P') - W_raw = 0 (integration by parts, no assumption)",
             *is_zero(sp.expand(r**2*eta*S - sp.diff(r**2*eta*P_E, r) - W_raw), (r, G)))

A = 4*Gam*P0 + r**2*sp.diff(P0, r)*sp.diff(rho0, r)/rho0 + 4*r*sp.diff(P0, r) - 4*sp.pi*G*rho0**2*r**2
C = r*(2*Gam*P0 + r*sp.diff(P0, r))
D = Gam*P0*r**2
W_sym = A*xi*eta + C*(xi*eta.diff(r) + xi.diff(r)*eta) + D*xi.diff(r)*eta.diff(r)
diffW = hydro(W_raw - W_sym)
# the difference contains only xi, xi', eta, eta', rho0, rho0', M, Gamma1, P0: treat them as independent symbols
x0, x1, e0, e1, rr0, rr1, MM, GG, PP = sp.symbols('x0 x1 e0 e1 rho0_ rho0p_ M_ Gam_ P0_')
repl = {sp.Derivative(xi, r): x1, sp.Derivative(eta, r): e1, sp.Derivative(rho0, r): rr1}
d2 = sp.expand(diffW).subs(repl).subs({xi: x0, eta: e0, rho0: rr0, M: MM, Gam: GG, P0: PP})
ok &= report("(2b) W_raw = A xi eta + C(xi eta' + xi' eta) + D xi' eta'  given only P0' = -G M rho0/r^2",
             sp.simplify(d2) == 0, f"residual = {sp.simplify(d2)}")

# symmetry of the bilinear form: swap (xi, xi') <-> (eta, eta')
Ws = sp.expand(W_sym).subs(repl).subs({xi: x0, eta: e0, rho0: rr0, M: MM, Gam: GG, P0: PP})
Ws_sw = Ws.subs({x0: e0, e0: x0, x1: e1, e1: x1}, simultaneous=True)
ok &= report("(2c) W is symmetric under xi <-> eta", sp.simplify(Ws - Ws_sw) == 0, f"difference = {sp.simplify(Ws - Ws_sw)}")

# (3) boundary term: with Delta P(R) = -Gamma1 P0 div = 0 at the surface the boundary value of r^2 eta P' is -R^2 eta xi P0'
bnd = (r**2*eta*(P_E + Gam*P0*div)).subs(r, R)
target = (-R**2*eta*xi*sp.diff(P0, r)).subs(r, R)
ok &= report("(3) [r^2 eta P']_R = -R^2 P0'(R) xi(R) eta(R) when Delta P(R) = 0", sp.simplify(bnd - target) == 0,
             f"difference = {sp.simplify(bnd - target)}")
# (3b) what the erroneous assumption P'(R)=0 would drop: the term xi(R) P0'(R) R^2 eta(R) is nonzero unless P0'(R)=0
print("note: the omitted term is -R^2 P0'(R) xi(R) eta(R); for a star with P0'(R) != 0 (finite-density surface) it is nonzero.")
sys.exit(0 if ok else 1)
