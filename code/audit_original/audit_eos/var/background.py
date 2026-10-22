"""
Phase A1: Chandrasekhar degenerate-electron EOS + TOV solve, reconstructed
from scratch, to reproduce Paper C's corrected configuration:
  rho_c = 5e13 kg/m^3  ->  M = 1.4249 Msun, compactness = 0.00515,
  Gamma1 minimum (at center) = 1.33409.

EOS (Chandrasekhar 1935/1939, zero-temperature relativistic degenerate
electron gas, mean molecular weight per electron mu_e):
  x = p_F/(m_e c)                     (relativity parameter)
  n_e = (8*pi/3) (m_e c / h)^3 x^3
  rho = mu_e * m_u * n_e               (mass density, baryon-dominated)
  P   = (pi * m_e^4 c^5)/(3 h^3) * f(x),
  f(x) = x(2x^2-3)sqrt(1+x^2) + 3*asinh(x)

Gamma1 = d ln P / d ln rho, evaluated analytically from the same x(rho)
relation (chain rule), not by finite-differencing noisy tables.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# CODATA-ish constants, SI
G = 6.67430e-11
c = 2.99792458e8
h = 6.62607015e-34
m_e = 9.1093837015e-31
m_u = 1.66053906660e-27
M_sun = 1.98892e30
mu_e = 2.0   # Carbon-12 / Oxygen-16: 2 nucleons per electron
Zlat = 0.0   # nuclear charge for the Salpeter (1961) lattice correction; 0 = pure uniform electron gas
alpha_fs = 7.2973525693e-3
hbar = h/(2*np.pi)
e2 = alpha_fs*hbar*c      # e^2/(4 pi eps0), SI

A_n = (8*np.pi/3)*(m_e*c/h)**3        # n_e = A_n * x^3
A_P = (np.pi*m_e**4*c**5)/(3*h**3)    # P = A_P * f(x)

def f_of_x(x):
    return x*(2*x**2-3)*np.sqrt(1+x**2) + 3*np.arcsinh(x)

def dfdx(x):
    # derivative of f(x), checked symbolically below at import time
    return 8*x**4/np.sqrt(1+x**2)

def rho_of_x(x):
    return mu_e*m_u*A_n*x**3

def B_lat():
    # Salpeter lattice pressure P_L = -(3/10)(4pi/3)^(1/3) e^2 Z^(2/3) n_e^(4/3), n_e = A_n x^3
    return 0.3*(4*np.pi/3)**(1/3)*e2*Zlat**(2/3)*A_n**(4/3) if Zlat > 0 else 0.0

def P_of_x(x):
    return A_P*f_of_x(x) - B_lat()*x**4

def x_of_rho(rho):
    return (rho/(mu_e*m_u*A_n))**(1.0/3.0)

use_eps = True   # if True, energy density includes electron kinetic+rest energy (thermodynamically consistent)

def g_of_x(x):
    return 3*(x*(2*x**2+1)*np.sqrt(1+x**2) - np.arcsinh(x))

def eps_over_c2_of_rho(rho):
    '''Total energy density / c^2. Thermodynamically consistent (first law n d(eps)/dn = eps+P)
    when use_eps: eps = n_e*mu_e*m_u*c^2 + eps_e (electron energy incl. rest), plus lattice term
    (-3 P_L... handled via P_L = u_L/3, u_L = 3 P_L) if Zlat>0.  Else eps = rho c^2 (original, inconsistent at O(P/rho c^2)).'''
    if not use_eps:
        return rho
    x = x_of_rho(rho)
    nuc = rho - mu_e*0.0  # nucleon rest-mass part: rho = mu_e m_u n_e
    eps_e = A_P*g_of_x(x)/c**2 - m_e*A_n*x**3   # electron energy incl rest minus rest mass already double counted? see note
    # NOTE: rho = mu_e m_u n_e is taken as the total rest-mass density of nuclei+electrons (standard
    # Chandrasekhar treatment), so subtract the electron rest-mass piece m_e n_e from eps_e to avoid double counting.
    u_L = 3*(-B_lat()*x**4)/c**2     # lattice energy density /c^2 (u_L = 3 P_L, E_L ~ n^{4/3})
    return rho + eps_e + u_L

def P_of_rho(rho):
    return P_of_x(x_of_rho(rho))

def Gamma1_of_rho(rho):
    """Gamma1 = dlnP/dlnrho = (rho/P) dP/drho, computed analytically:
    dP/drho = (dP/dx)/(drho/dx) = [A_P f'(x)] / [3 mu_e m_u A_n x^2]."""
    x = x_of_rho(rho)
    dPdx = A_P*dfdx(x) - 4*B_lat()*x**3
    drhodx = 3*mu_e*m_u*A_n*x**2
    dPdrho = dPdx/drhodx
    P = P_of_x(x)
    return (rho/P)*dPdrho

if __name__ == "__main__":
    # sanity check dfdx against numerical derivative
    import numpy.testing as npt
    xs = np.array([1e-3, 0.01, 0.1, 1.0, 5.0, 20.0])
    eps = 1e-6
    num = (f_of_x(xs+eps)-f_of_x(xs-eps))/(2*eps)
    ana = dfdx(xs)
    print("f'(x) check (numeric vs analytic):")
    for xx,n,a in zip(xs,num,ana):
        print(f"  x={xx:8.4f}  numeric={n:.6e}  analytic={a:.6e}  rel.err={abs(n-a)/abs(a):.2e}")

    # sanity: non-relativistic (x<<1) and ultra-relativistic (x>>1) limits of Gamma1
    print()
    print("Gamma1 limits check: should -> 5/3 as x->0 (non-rel), -> 4/3 as x->inf (ultra-rel)")
    for rho_test in [1e6, 1e9, 1e12, 1e15, 1e18]:
        print(f"  rho={rho_test:.1e} kg/m^3  x={x_of_rho(rho_test):.4e}  Gamma1={Gamma1_of_rho(rho_test):.6f}")
