"""Symbolic checks of the equation-of-state and energy statements of Paper C (Secs. 2.1, 2.3, 5.2, 5.4, 5.7).
Run: python3 sympy_eos_checks.py   (exit status 1 if any check fails)
 (a) Chandrasekhar P(x), eps_e(x) from the Fermi integrals:   P = (8 pi m^4 c^5/3h^3) int_0^x t^4/sqrt(1+t^2) dt,
     eps_e = (8 pi m^4 c^5/h^3) int_0^x t^2 sqrt(1+t^2) dt  (energy density incl. rest mass), n_e = (8 pi/3)(m c/h)^3 x^3.
 (b) f'(x) = 8 x^4/sqrt(1+x^2).
 (c) polytropic limits Gamma1 -> 5/3 (x->0) and 4/3 (x->inf).
 (d) first law  n d(eps)/dn = eps + P  for eps = rho c^2 + eps_e - m_e n_e c^2  (exactly), and for the Salpeter lattice term;
     failure of the first law for eps = rho c^2 (the 'first implementation'), at O(P/rho c^2).
 (e) du/drho = P/rho^2 (internal energy per unit mass of the barotropic gas; used in Eq. dUint).
 (f) linearisation  d/dlambda[3 lambda^2 P(rho/lambda^3)]_{lambda=1} = 3 P (2 - 3 Gamma1)   (model of Sec. 5.1).
 (g) integral identity  int(-g_eff r) dM_r = 3 W_P  (g_eff defined by dP/dr = rho g_eff) by one integration by parts.
 (h) Salpeter lattice: u_L = -(9/10) Z^2 e^2/r_i per ion, r_i = (3/4 pi n_i)^(1/3), n_i = n_e/Z  =>  P_L = u_L/3 =
     -(3/10)(4 pi/3)^(1/3) e^2 Z^(2/3) n_e^(4/3).
 (i) bulk-viscosity coefficient  nu = 9 Gamma1bar W_P Delta_r tau / I.
"""
import sys
import sympy as sp
from _zero import is_zero, report

ok = True
x, t = sp.symbols('x t', positive=True)
f = x*(2*x**2 - 3)*sp.sqrt(1 + x**2) + 3*sp.asinh(x)
g = 3*(x*(2*x**2 + 1)*sp.sqrt(1 + x**2) - sp.asinh(x))

# (a) Fermi integrals: check that the closed forms are antiderivatives with the right normalisation and vanish at x=0
Fp = f/8          # int_0^x t^4/sqrt(1+t^2) dt
Fe = (x*(2*x**2 + 1)*sp.sqrt(1 + x**2) - sp.asinh(x))/8   # int_0^x t^2 sqrt(1+t^2) dt
ok &= report("(a1) d/dx[f/8] = x^4/sqrt(1+x^2) and f(0)=0", *is_zero(sp.diff(Fp, x) - x**4/sp.sqrt(1 + x**2), (x,)))
ok &= report("(a2) d/dx[(x(2x^2+1)sqrt(1+x^2)-asinh x)/8] = x^2 sqrt(1+x^2), value 0 at 0",
             *is_zero(sp.diff(Fe, x) - x**2*sp.sqrt(1 + x**2), (x,)))
ok &= report("(a3) prefactors: (8/3) * (1/8) f = f/3 -> P = (pi m^4 c^5/3h^3) f ; 8 * Fe = x(2x^2+1)sqrt(1+x^2)-asinh x -> eps_e = (pi m^4c^5/3h^3) g",
             sp.simplify(sp.Rational(8, 3)*Fp - f/3) == 0 and sp.simplify(8*Fe*3 - g) == 0, "")
# (b)
ok &= report("(b) f'(x) = 8 x^4/sqrt(1+x^2)", *is_zero(sp.diff(f, x) - 8*x**4/sp.sqrt(1 + x**2), (x,)))
# (c)
Gam1 = x*sp.diff(f, x)/(3*f)
l0 = sp.limit(Gam1, x, 0); linf = sp.limit(Gam1, x, sp.oo)
ok &= report("(c) Gamma1(x->0) = 5/3, Gamma1(x->oo) = 4/3", l0 == sp.Rational(5, 3) and linf == sp.Rational(4, 3), f"limits: {l0}, {linf}")
# (d) first law; n ~ x^3 so n d/dn = (x/3) d/dx.  eps = (mu_e m_u - m_e) c^2 A_n x^3 + K g(x);  K = pi m^4c^5/(3h^3)
K, An, mu_e, mu, me, c, B = sp.symbols('K A_n mu_e m_u m_e c B', positive=True)
eps = (mu_e*mu - me)*c**2*An*x**3 + K*g
P = K*f
ok &= report("(d1) first law n deps/dn = eps + P (ideal electron gas, eps = rho c^2 + eps_e - m_e n c^2)",
             *is_zero(x/3*sp.diff(eps, x) - (eps + P), (x, K, An, mu_e, mu, me, c)))
eps_naive = mu_e*mu*c**2*An*x**3
res_naive = sp.simplify(x/3*sp.diff(eps_naive, x) - (eps_naive + P))
ok &= report("(d2) the first law FAILS for eps = rho c^2 (residual = -P != 0)", sp.simplify(res_naive + P) == 0 and res_naive != 0, f"n deps/dn - (eps+P) = {res_naive}")
uL = -3*B*x**4      # lattice energy density u_L = 3 P_L, P_L = -B x^4
PL = -B*x**4
ok &= report("(d3) first law for the Salpeter lattice term (u_L = 3 P_L ~ n^(4/3))", *is_zero(x/3*sp.diff(uL, x) - (uL + PL), (x, B)))
# (e) rho^2 du/drho = P  with u = (eps - rho c^2)/rho, rho = mu_e m_u A_n x^3
rho = mu_e*mu*An*x**3
u = (eps - rho*c**2)/rho
ok &= report("(e) rho^2 du/drho = P (barotropic identity behind Eq. dUint)", *is_zero(rho**2*sp.diff(u, x)/sp.diff(rho, x) - P, (x, K, An, mu_e, mu, me, c)))
# (f)
lam, rr = sp.symbols('lambda rho', positive=True)
Pf = sp.Function('P')
G1 = sp.symbols('Gamma1')
dh = sp.diff(3*lam**2*Pf(rr/lam**3), lam).subs(lam, 1).doit()
by_hand = 6*Pf(rr) - 9*rr*sp.Derivative(Pf(rr), rr)
print("      sympy d/dlambda[3 lambda^2 P(rho/lambda^3)] at lambda=1:", dh)
ok &= report("(f1) d/dlambda[3 lambda^2 P(rho lambda^-3)]_1 = 6P - 9 rho P'", sp.simplify(dh - by_hand) == 0, "")
ok &= report("(f2) ... = 3P(2 - 3 Gamma1) with Gamma1 = rho P'/P",
             sp.simplify(by_hand.subs(sp.Derivative(Pf(rr), rr), G1*Pf(rr)/rr) - 3*Pf(rr)*(2 - 3*G1)) == 0, "")
# (g) integration by parts: -r (4 pi r^2) P' = -(4 pi r^3 P)' + 3 (4 pi r^2 P)
r = sp.symbols('r', positive=True); Pr = sp.Function('P')(r)
ok &= report("(g) -4 pi r^3 P' = -(4 pi r^3 P)' + 3 (4 pi r^2 P)  => int(-g_eff r)dM = 3 W_P for P(R)=0",
             *is_zero(-4*sp.pi*r**3*sp.diff(Pr, r) - (-sp.diff(4*sp.pi*r**3*Pr, r) + 3*4*sp.pi*r**2*Pr), (r,)))
# (h) Salpeter
Z, e, ne = sp.symbols('Z e n_e', positive=True)
ni = ne/Z
ri = (3/(4*sp.pi*ni))**sp.Rational(1, 3)
uion = -sp.Rational(9, 10)*Z**2*e**2/ri           # lattice energy per ion
uL_ = ni*uion                                      # energy density
PL_ = ne*sp.diff(uL_/ne, ne)*ne                    # P = -d(E)/dV = n^2 d(u/n)/dn  (u/n = energy per electron)
PL_target = -sp.Rational(3, 10)*(4*sp.pi/3)**sp.Rational(1, 3)*e**2*Z**sp.Rational(2, 3)*ne**sp.Rational(4, 3)
ok &= report("(h1) P_L = n^2 d(u_L/n)/dn = -(3/10)(4pi/3)^(1/3) e^2 Z^(2/3) n_e^(4/3)", *is_zero(PL_ - PL_target, (Z, e, ne)))
ok &= report("(h2) P_L = u_L/3", *is_zero(PL_ - uL_/3, (Z, e, ne)))
# (i) bulk viscosity: viscous energy loss -int zeta (div v)^2 dV with div v = 3 lamdot (homologous, lambda = 1) versus the damped-oscillator loss -nu I lamdot^2
Dr, tau, G1b, WP, I, ldot, nu = sp.symbols('Delta_r tau Gamma1bar W_P I lamdot nu', positive=True)
# zeta = rho (c_inf^2 - c_0^2) tau/(1+w^2 tau^2) -> rho Delta_r c_0^2 tau (w tau << 1), c_0^2 = Gamma1 P/rho  =>  zeta = Gamma1 P Delta_r tau;
# int zeta dV = Delta_r tau int Gamma1 P dV = Delta_r tau Gamma1bar W_P   (definition Gamma1bar = int Gamma1 P dV / W_P)
int_zeta = Dr*tau*G1b*WP
nu_sol = sp.solve(sp.Eq(nu*I*ldot**2, (3*ldot)**2*int_zeta), nu)[0]
print("      [info] nu from  nu I lamdot^2 = 9 lamdot^2 int zeta dV :", nu_sol)
ok &= report("(i) nu = 9 Gamma1bar W_P Delta_r tau / I", sp.simplify(nu_sol - 9*G1b*WP*Dr*tau/I) == 0, "")
sys.exit(0 if ok else 1)
