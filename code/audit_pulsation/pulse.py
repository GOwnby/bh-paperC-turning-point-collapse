import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from eos import *
L=1.0e6   # length unit (m)
fG=G/c**4 # J/m^3 -> m^-2
def egeo(x): return fG*eps(x)*L**2     # energy density in L^-2
def pgeo(x): return fG*P(x)*L**2
def dpdx_geo(x): return fG*dPdx(x)*L**2

def nu1(r,m,Pg): return (m+4*np.pi*r**3*Pg)/(r*(r-2*m))

def bg_rhs(r,y):
    m,nu,x=y
    e=egeo(x); Pg=pgeo(x)
    n1=nu1(r,m,Pg)
    dPdr=-(e+Pg)*n1
    return [4*np.pi*r*r*e, n1, dPdr/dpdx_geo(x)]

def background(rho_c,x_end=1e-3,r0=1e-6,rtol=1e-12,atol=1e-14):
    xc=x_of_rho(rho_c); e=egeo(xc); Pg=pgeo(xc)
    m0=4*np.pi/3*e*r0**3
    x0=xc-(4*np.pi/3)*(e+Pg)*(e+3*Pg)*r0**2/2/dpdx_geo(xc)*(-1)*(-1)
    # dP/dr = -(4pi/3)(e+P)(e+3P) r  => P(r0)=Pc-(2pi/3)(e+P)(e+3P) r0^2
    x0=xc-(2*np.pi/3)*(e+Pg)*(e+3*Pg)*r0**2/dpdx_geo(xc)
    nu0=(2*np.pi/3)*(e+3*Pg)*r0**2
    ev=lambda r,y: y[2]-x_end
    ev.terminal=True; ev.direction=-1
    sol=solve_ivp(bg_rhs,[r0,50.0],[m0,nu0,x0],method='DOP853',rtol=rtol,atol=atol,events=ev,dense_output=True)
    R=sol.t_events[0][0]; yR=sol.y_events[0][0]
    return sol,R,yR,xc

def surf_shift(R,M,nuR):
    return 0.5*np.log(1-2*M/R)-nuR   # s

def sigma2_scale(R,M,nuR):
    s=surf_shift(R,M,nuR)
    # omega_true^2 = omega_tmp^2 * e^{2s}; omega_tmp^2 in L^-2 ; sigma^2[s^-2]=omega^2 c^2 / L^2
    return np.exp(2*s)*c**2/L**2

def coeffsA(r,m,nu,x):
    e=egeo(x); Pg=pgeo(x); G1=Gamma1(x)
    e2l=1/(1-2*m/r); el=np.sqrt(e2l)
    n1=nu1(r,m,Pg)
    Pk=G1*Pg*el*np.exp(3*nu)/r**2
    W=(e+Pg)*el**3*np.exp(nu)/r**2
    Q=el*np.exp(3*nu)*(e+Pg)*(n1**2+4*n1/r-8*np.pi*e2l*Pg)/r**2
    return Pk,W,Q

def rhsA(w2):
    def f(r,y):
        m,nu,x,z,eta=y
        d=bg_rhs(r,y[:3])
        Pk,W,Q=coeffsA(r,m,nu,x)
        return d+[eta/Pk, -(Q+w2*W)*z]
    return f

def shootA(rho_c,w2,x_end=1e-3,r0=1e-6,rtol=1e-12,atol=1e-14,full=False):
    xc=x_of_rho(rho_c); Pg=pgeo(xc); G1=Gamma1(xc); e=egeo(xc)
    m0=4*np.pi/3*e*r0**3
    x0=xc-(2*np.pi/3)*(e+Pg)*(e+3*Pg)*r0**2/dpdx_geo(xc)
    nu0=(2*np.pi/3)*(e+3*Pg)*r0**2
    z0=r0**3; eta0=3*G1*Pg
    ev=lambda r,y: y[2]-x_end
    ev.terminal=True; ev.direction=-1
    sol=solve_ivp(rhsA(w2),[r0,50.0],[m0,nu0,x0,z0,eta0],method='DOP853',rtol=rtol,atol=atol,events=ev,dense_output=full)
    yR=sol.y_events[0][0]; R=sol.t_events[0][0]
    F=yR[4]/(3*G1*Pg)   # eta(r_end)/eta0 : regular solution -> 0, singular -> const
    return (F,R,yR,sol) if full else F
