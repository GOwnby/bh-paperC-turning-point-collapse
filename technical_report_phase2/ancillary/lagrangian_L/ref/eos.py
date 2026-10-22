import numpy as np

# Physical constants (SI)
G = 6.674e-11; c = 2.998e8; hbar=1.0546e-34; h = 2*np.pi*hbar
m_e = 9.109e-31; m_H = 1.6605e-27; Msun = 1.989e30

mu_e = 2.0   # mean molecular weight per electron (Carbon/Oxygen white dwarf)

# Chandrasekhar 1935 degenerate relativistic electron-gas EOS, parametrized by x = p_F/(m_e c)
def f_cs(x):
    return x*(2*x**2-3)*np.sqrt(1+x**2) + 3*np.arcsinh(x)
def g_cs(x):   # for the electron internal energy contribution (small correction to rest-mass density)
    return 8*x**3*(np.sqrt(1+x**2)-1)

def rho_of_x(x):
    n_e = (8*np.pi/3)*(m_e*c/h)**3 * x**3
    return mu_e*m_H*n_e     # rest-mass density dominated by baryons

def P_of_x(x):
    return (np.pi*m_e**4*c**5)/(3*h**3) * f_cs(x)

# sanity check: recover the known non-relativistic (x<<1) and ultra-relativistic (x>>1) polytropic limits
x_small = 1e-3; x_large = 1e3
P_s, rho_s = P_of_x(x_small), rho_of_x(x_small)
P_l, rho_l = P_of_x(x_large), rho_of_x(x_large)
# Gamma = dlnP/dlnrho numerically
x1,x2 = 1e-3, 1.1e-3
Gamma_NR = (np.log(P_of_x(x2))-np.log(P_of_x(x1)))/(np.log(rho_of_x(x2))-np.log(rho_of_x(x1)))
x1,x2 = 1e3, 1.1e3
Gamma_UR = (np.log(P_of_x(x2))-np.log(P_of_x(x1)))/(np.log(rho_of_x(x2))-np.log(rho_of_x(x1)))
print(f"Non-relativistic limit (x={x_small}): effective Gamma = {Gamma_NR:.4f}  (expect 5/3 = {5/3:.4f})")
print(f"Ultra-relativistic limit (x={x_large}): effective Gamma = {Gamma_UR:.4f}  (expect 4/3 = {4/3:.4f})")

# central density / pressure at a representative white-dwarf-core value
x_test = 10.0
print(f"\nAt x={x_test}: rho={rho_of_x(x_test):.4e} kg/m^3, P={P_of_x(x_test):.4e} Pa")
print(f"Nuclear density for reference: ~2.3e17 kg/m^3")
print(f"x needed to reach nuclear density:")
from scipy.optimize import brentq
x_nuc = brentq(lambda x: rho_of_x(x)-2.3e17, 1, 1000)
print(f"  x = {x_nuc:.3f}  (this is well past the ultra-relativistic regime; Chandrasekhar's electron-only")
print(f"   EOS is not physically trustworthy this far -- neutronization occurs first, around")
print(f"   rho ~ 1e12-1e13 kg/m^3 for typical WD compositions -- flagged explicitly, see notes below)")
