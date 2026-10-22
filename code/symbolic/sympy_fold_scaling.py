"""Symbolic checks of the fold-passage (Painleve I) scaling and the dissipation (Airy) law of Paper C, Secs. 5.6-5.7.
Run: python3 sympy_fold_scaling.py   (exit status 1 if any check fails; the numerical pole values are in pi_tritronquee.py,
pi_boutroux.py and damped_pi.py)
 (1) normal form  y'' = (kappa/2) y^2 + kappa/(2a) (M - Mmax), M - Mmax = Mdot t: equilibria and sigma^2 = -kappa y on both branches.
 (2) scaling y = alpha u, t = beta s  =>  u'' = u^2 + s  with  beta = (4a/(kappa^2 Mdot))^(1/5), alpha = 2/(kappa beta^2).
 (3) slow manifold of u'' = u^2 + s: u = -sqrt(-s) - 1/(8 s^2) + ...   (the 'first adiabatic correction').
 (4) Boutroux map: u = 6 B^2 Y, s = x/B, B = 6^(-1/5) turns u'' = u^2 + s into Y'' = 6 Y^2 + x  (so s_* = 6^(1/5) x_*).
 (5) damped form u'' + nu~ u' = u^2 + s, nu~ = nu beta;  overdamped scaling beta' = beta nu~^(2/3) = (4 a nu^2/(kappa^2 Mdot))^(1/3) -> w' = w^2 + tau.
 (6) w = -phi'/phi turns w' = w^2 + tau into phi'' = -tau phi, solved by phi = Ai(-tau); first pole at the first zero of Ai(-tau), tau_* = -a_1 = 2.33811.
"""
import sys
import sympy as sp
import mpmath as mp
from _zero import is_zero, report

ok = True
y, t, kap, a, Md, al, be, s = sp.symbols('y t kappa a Mdot alpha beta s', positive=True)
# (1)
Mm = sp.symbols('DM')    # M - Mmax
yeq = sp.symbols('y_e', real=True)
F = kap/2*y**2 + kap/(2*a)*Mm
Meq = sp.solve(sp.Eq(F.subs(y, yeq), 0), Mm)[0]
ok &= report("(1a) equilibria of y''=0: M - Mmax = -a y^2", sp.simplify(Meq + a*yeq**2) == 0, f"M-Mmax = {Meq}")
dF = sp.diff(F, y).subs(y, yeq)
ok &= report("(1b) d(y'')/dy at equilibrium = kappa y = -sigma^2 with sigma^2 = -kappa y (both branches)", sp.simplify(dF - kap*yeq) == 0, f"= {dF}")
# (2)
u = sp.Function('u')
# y(t) = alpha u(t/beta)  =>  y'' = (alpha/beta^2) u''(s)
lhs_coeff = al/be**2
c_u2 = (kap/2*al**2)/lhs_coeff
c_s = (kap/(2*a)*Md*be)/lhs_coeff            # M - Mmax = Mdot t = Mdot beta s
sol = sp.solve([sp.Eq(c_u2, 1), sp.Eq(c_s, 1)], [al, be], dict=True)
print("      solutions:", sol)
good = False
for so in sol:
    if sp.simplify(so[be]**5 - 4*a/(kap**2*Md)) == 0 and sp.simplify(so[al] - 2/(kap*so[be]**2)) == 0:
        good = True
ok &= report("(2) y = alpha u, t = beta s gives u'' = u^2 + s iff beta^5 = 4a/(kappa^2 Mdot), alpha = 2/(kappa beta^2)", good, "")
# (3) slow manifold: s = -X^2, d/ds = -(1/(2X)) d/dX,  U(X) = -X + c/X^4  (= -sqrt(-s) + c/s^2)
X = sp.Symbol('X', positive=True)
c = sp.Symbol('c')
D = lambda F_: -sp.diff(F_, X)/(2*X)
U = -X + c/X**4
resid = sp.expand(D(D(U)) - U**2 - (-X**2))
lead = sp.expand(resid*X**3).subs(X, sp.oo) if False else sp.limit(sp.expand(resid*X**3), X, sp.oo)
csol = sp.solve(sp.Eq(lead, 0), c)
print("      residual u''-u^2-s of u = -sqrt(-s) + c/s^2:", resid, " ; leading X^-3 coefficient:", lead)
ok &= report("(3) first adiabatic correction of u''=u^2+s on the slow manifold is c = -1/8 (u = -sqrt(-s) - 1/(8 s^2))",
             csol == [-sp.Rational(1, 8)], f"c = {csol}")
# (4) Boutroux
Bc = sp.Rational(6)**sp.Rational(-1, 5)
x_, Y = sp.symbols('x', positive=True), sp.Function('Y')
Yx = Y(x_)
uu = 6*Bc**2*Yx            # u as function of x, with s = x/B => d/ds = B d/dx
upp = Bc**2*sp.diff(uu, x_, 2)
resid = sp.simplify(upp - (uu**2 + x_/Bc))      # = 6 B^4 [Y'' - 6 Y^2 - x/(6 B^5)]
target = 6*Bc**4*(sp.diff(Yx, x_, 2) - 6*Yx**2 - x_)
ok &= report("(4) u''=u^2+s with u = 6B^2 Y, s = x/B, B=6^(-1/5) is equivalent to Y'' = 6Y^2 + x", sp.simplify(resid - target) == 0, "")
ok &= report("(4b) s_* = x_*/B = 6^(1/5) x_*", sp.simplify(1/Bc - sp.Rational(6)**sp.Rational(1, 5)) == 0, f"6^(1/5) = {sp.N(sp.Rational(6)**sp.Rational(1,5), 10)}")
# (5)
nu, nut = sp.symbols('nu nu_t', positive=True)
A_, B_ = sp.symbols('A_ B_', positive=True)
# u = A w(tau), s = B tau :  nu~ (A/B) w' = A^2 w^2 + B tau   -> need nu~ A/B = A^2 = B
sol5 = sp.solve([sp.Eq(nut*A_/B_, A_**2), sp.Eq(A_**2, B_)], [A_, B_], dict=True)
print("      overdamped rescaling solutions:", sol5)
good5 = any(sp.simplify(so[B_] - nut**sp.Rational(2, 3)) == 0 and sp.simplify(so[A_] - nut**sp.Rational(1, 3)) == 0 for so in sol5)
ok &= report("(5a) overdamped: u = nu~^(1/3) w, s = nu~^(2/3) tau gives w' = w^2 + tau", good5, "")
bet = (4*a/(kap**2*Md))**sp.Rational(1, 5)
bprime = bet*(nu*bet)**sp.Rational(2, 3)
ok &= report("(5b) beta' = beta nu~^(2/3) = (4 a nu^2/(kappa^2 Mdot))^(1/3)",
             *is_zero(bprime**3 - 4*a*nu**2/(kap**2*Md), (a, kap, Md, nu)))
# (6)
tau = sp.symbols('tau', real=True)
phi = sp.airyai(-tau)
w = -sp.diff(phi, tau)/phi
ok &= report("(6a) w = -phi'/phi with phi = Ai(-tau) satisfies w' = w^2 + tau (uses Ai''(z) = z Ai(z))",
             *is_zero(sp.simplify(sp.diff(w, tau) - w**2 - tau), (tau,)))
mp.mp.dps = 30
a1 = mp.airyaizero(1)
print("      first zero of Ai: a_1 =", a1, " => tau_* = -a_1 =", -a1)
ok &= report("(6b) tau_* = -a_1 = 2.33811 (DLMF 9.9)", abs(-a1 - mp.mpf('2.338107410459767')) < 1e-12, "")
# t_collapse -> 2.338 beta nu~^(2/3)
sys.exit(0 if ok else 1)
