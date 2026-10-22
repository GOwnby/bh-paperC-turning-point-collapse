"""Numbers quoted in Paper C Secs. 2.1, 2.3 and 3.6 about the equation of state (complements tov_solve.py and misc_checks.py).  ~2 s.
 (1) polytropic limits of the Chandrasekhar EOS: Gamma1 -> 5/3 (x -> 0) and 4/3 (x -> oo), to four significant figures;
 (2) the closed-form Newtonian Chandrasekhar mass M = (w_3 sqrt(3 pi)/2) (hbar c/G)^(3/2)/(mu_e m_u)^2 with w_3 = xi_1^2 |theta'(xi_1)| from the
     n = 3 Lane-Emden solution (computed here), compared with the 1.456 M_sun asymptote of the Newtonian sequence of the paper (Sec. 2.1);
 (3) at the fiducial rho_c = 5e13 kg/m^3 and at rho_max = 2.3106e13: P/(rho c^2), eps/(rho c^2) - 1 (first-law inconsistency of eps = rho c^2),
     Gamma1 - 4/3 (the margin of the instability), Sec. 2.3 and 3.6.
"""
import numpy as np
from scipy.integrate import solve_ivp
import background as bg
c = bg.c
xrho = lambda x: bg.mu_e*bg.m_u*bg.A_n*x**3
# f(x) = x(2x^2-3)sqrt(1+x^2) + 3 asinh(x) cancels catastrophically for x << 1 (relative error ~1e-16/x^4; harmless for x >~ 0.05, which covers the
# whole TOV/pulsation integration, whose surface is at P = 1e-10 P_c, x ~ 0.1), so the small-x limit is evaluated with the convergent binomial series
# f(x) = 8 sum_k binom(-1/2,k) x^(2k+5)/(2k+5)  (|x| < 1), the same series as in audit_pulsation/eos.py.
from scipy.special import binom
def Gamma1_series(x, nt=40):
    k = np.arange(nt); cP = 8*binom(-0.5, k)/(2*k + 5)
    f = np.sum(cP*x**(2*k + 5)); fp = np.sum(cP*(2*k + 5)*x**(2*k + 4))
    return x*fp/(3*f)
print('Gamma1(x=1e-5) = %.6f (series; 5/3 = %.6f)   Gamma1(x=1e5) = %.6f (closed form; 4/3 = %.6f)' % (Gamma1_series(1e-5), 5/3, bg.Gamma1_of_rho(xrho(1e5)), 4/3))
print('Gamma1(x=0.1): closed form %.8f, series %.8f (agree: no cancellation problem in the range used)' % (bg.Gamma1_of_rho(xrho(0.1)), Gamma1_series(0.1)))
f = lambda x, y: [y[1], -max(y[0], 0)**3 - 2*y[1]/x]
ev = lambda x, y: y[0]; ev.terminal = True
x0 = 1e-6
s = solve_ivp(f, [x0, 20], [1 - x0**2/6, -x0/3], events=ev, rtol=1e-13, atol=1e-15)
xi1 = s.t_events[0][0]; thp = s.y_events[0][0][1]; w3 = xi1**2*(-thp)
mu = bg.mu_e*bg.m_u
Mch = w3*np.sqrt(3*np.pi)/2*(bg.hbar*c/bg.G)**1.5/mu**2/bg.M_sun
print('Lane-Emden n=3: xi_1 = %.5f, w_3 = %.5f ; Newtonian Chandrasekhar mass = %.5f M_sun (mu_e = 2, constants of background.py)' % (xi1, w3, Mch))
for rho in (2.3106e13, 5e13):
    P = bg.P_of_rho(rho); e = bg.eps_over_c2_of_rho(rho)
    print('rho = %.4e: P/(rho c^2) = %.4e ; eps/(rho c^2) - 1 = %.4e ; Gamma1 - 4/3 = %.4e (centre value)' % (rho, P/(rho*c**2), e/rho - 1, bg.Gamma1_of_rho(rho) - 4/3))
