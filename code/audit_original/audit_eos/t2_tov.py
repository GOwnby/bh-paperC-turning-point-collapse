# Independent TOV re-implementation (written from scratch from the Chandrasekhar EOS formulas; does NOT import the code under audit)
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
G=6.67430e-11; c=2.99792458e8; h=6.62607015e-34; me=9.1093837015e-31; mu=1.66053906660e-27
Msun_code=1.98892e30; Msun_GM=1.32712440018e20/G
mue=2.0
lam_c=h/(me*c)             # Compton wavelength (h-bar-less)
nfac=8*np.pi/(3*lam_c**3)    # n_e = nfac x^3
Pfac=np.pi*me*c**2/(3*lam_c**3)   # = pi me^4 c^5/(3 h^3)
def f(x): return x*(2*x*x-3)*np.sqrt(1+x*x)+3*np.arcsinh(x)
def fser(x):  # series for small x to avoid cancellation: f = 8/5 x^5 - 4/7 x^7 + 1/9 x^9 ... (use only if x<0.05)
    return 8/5*x**5-4/7*x**7+1/9*x**9-5/88*x**11
def Pf(x): return Pfac*np.where(x<0.05,fser(x),f(x))
def g(x): return 3*(x*(2*x*x+1)*np.sqrt(1+x*x)-np.arcsinh(x))
def eps_nokin(x): return mue*mu*nfac*x**3*c**2
def eps_tot(x, mode='kin'):
    ne=nfac*x**3
    if mode=='rho': return mue*mu*ne*c**2
    if mode=='kin':   # atomic-mass rho incl. electron rest mass; add electron KINETIC energy only
        return mue*mu*ne*c**2 + Pfac*g(x) - me*c**2*ne
    if mode=='double': # naive: rho c^2 + full electron energy incl. rest mass (double counts m_e)
        return mue*mu*ne*c**2 + Pfac*g(x)
    if mode=='nuc':   # rho_atomic - m_e ne  (nuclei only) + full electron energy (equivalent to 'kin' to O(1)) 
        return (mue*mu*ne - me*ne)*c**2 + Pfac*g(x)
def dPdx(x): return Pfac*8*x**4/np.sqrt(1+x*x)
def solve(rhoc, mode='kin', Pcut=1e-10, rtol=1e-12, return_profile=False, mu_=None):
    xc=(rhoc/(mue*mu*nfac))**(1/3); Pc=Pf(xc)
    r0=1.0
    def rhs(r,y):
        x,m=y
        e=eps_tot(x,mode)/c**2; P=Pf(x)
        dP=-G*(e+P/c**2)*(m+4*np.pi*r**3*P/c**2)/(r**2*(1-2*G*m/(r*c**2)))
        return [dP/dPdx(x), 4*np.pi*r*r*e]
    xs=(Pcut*Pc/ (Pfac*8/5))**(1/5)   # NR estimate for surface, refine by event on P itself
    ev=lambda r,y: Pf(y[0])-Pcut*Pc; ev.terminal=True
    m0=4/3*np.pi*r0**3*eps_tot(xc,mode)/c**2
    s=solve_ivp(rhs,[r0,5e7],[xc,m0],events=ev,rtol=rtol,atol=[1e-14,1e-6],method='DOP853')
    return s.t_events[0][0], s.y_events[0][0][1]
if __name__=="__main__":
    print("M_sun: code 1.98892e30 ; GM_sun/G =",Msun_GM)
    R,M=solve(5e13)
    print("rho_c=5e13: R=%.4f km  M=%.5f Msun(code Msun) = %.5f (GM-based Msun)  compactness 2GM/Rc2=%.5f"%(R/1e3,M/Msun_code,M/Msun_GM,2*G*M/(R*c**2)))
    for cut in (1e-8,1e-10,1e-12,1e-14,1e-16):
        R,M=solve(5e13,Pcut=cut); print(" Pcut=%g R=%.4f km M=%.6f"%(cut,R/1e3,M/Msun_code))
    for mode in ('rho','kin','double','nuc'):
        res=minimize_scalar(lambda lr:-solve(10**lr,mode)[1],bounds=(13.0,13.8),method='bounded',options={'xatol':1e-5})
        Mm=-res.fun
        R5,M5=solve(5e13,mode)
        print("mode %-6s: M_max=%.5f Msun(code) / %.5f (GM Msun) at rho_c=%.4e ; at 5e13: M=%.5f R=%.3f km"%(mode,Mm/Msun_code,Mm/Msun_GM,10**res.x,M5/Msun_code,R5/1e3))
