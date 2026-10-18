import numpy as np
from scipy.special import binom
# CODATA 2018 (SI)
h=6.62607015e-34; c=299792458.0; me=9.1093837015e-31; mu=1.66053906660e-27
G=6.67430e-11; Msun=1.98841e30; mue=2.0
K=np.pi*me**4*c**5/(3*h**3)      # J/m^3
A=(8*np.pi/3)*(me*c/h)**3        # n_e = A x^3
XS=0.3; NT=22
kk=np.arange(NT)
cP=8*binom(-0.5,kk)/(2*kk+5)      # P/K = sum cP x^(2k+5)
cE=24*binom(0.5,kk)/(2*kk+3)      # Ee/K = sum cE x^(2k+3)
def n_e(x): return A*x**3
def rho(x): return mue*mu*A*x**3
def x_of_rho(r): return (r/(mue*mu*A))**(1/3)
def P(x):
    x=np.asarray(x,float)
    big=K*(x*(2*x*x-3)*np.sqrt(1+x*x)+3*np.arcsinh(x))
    xs=np.minimum(x,XS)
    ser=K*sum(cP[k]*xs**(2*k+5) for k in range(NT))
    return np.where(x<XS,ser,big)
def Ee(x):
    x=np.asarray(x,float)
    big=K*3*(x*(2*x*x+1)*np.sqrt(1+x*x)-np.arcsinh(x))
    xs=np.minimum(x,XS)
    ser=K*sum(cE[k]*xs**(2*k+3) for k in range(NT))
    return np.where(x<XS,ser,big)
def eps(x): return rho(x)*c**2 + Ee(x) - me*n_e(x)*c**2
def dPdx(x): return K*8*x**4/np.sqrt(1+x*x)
def Gamma1(x): return x/(3*P(x))*dPdx(x)   # dlnP/dln rho  (rho ∝ x^3)
