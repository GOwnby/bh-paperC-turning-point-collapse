import sys; sys.path.insert(0,'.')
from eos import *
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import numpy as np
# SI Newtonian: state [m, x, xi, dP] with xi=Delta r / r, dP = Delta P ; use standard first-order form
def run(rc,w2,x_end=1e-3,r0=1.0):
    xc=x_of_rho(rc); Pc=P(xc)
    def f(r,y):
        m,x,xi,dP=y
        rho_=rho(x); Pp=-G*rho_*m/r**2; G1=Gamma1(x); Pn=P(x)
        dx=Pp/dPdx(x)
        dxi=-(3*xi+dP/(G1*Pn))/r
        ddP=xi*(w2*rho_*r-4*Pp)*1.0 + 0  # Newtonian: dDP/dr = xi(w2 rho r - 4 P')  (+ 0*dP)
        return [4*np.pi*r*r*rho_,dx,dxi,ddP]
    rho0=rho(xc)
    m0=4*np.pi/3*rho0*r0**3
    x0=xc-(2*np.pi/3)*G*rho0**2*r0**2/dPdx(xc)
    ev=lambda r,y:y[1]-x_end; ev.terminal=True; ev.direction=-1
    sol=solve_ivp(f,[r0,1e9],[m0,x0,1.0,-3*Gamma1(xc)*Pc],method='DOP853',rtol=1e-12,atol=1e-30,events=ev)
    y=sol.y_events[0][0]
    return y[3]/(3*Gamma1(xc)*Pc), sol.t_events[0][0], y[0]
if __name__=="__main__":
    for rc in [float(a) for a in sys.argv[1:]]:
        g=[-3e3,-300,-30,-3,0,3,30,300,3e3,1e4]
        # scale
        v=[run(rc,s)[0] for s in g]
        s0=np.sign(v[0]); 
        for i in range(len(g)-1):
            if s0*v[i]>0 and s0*v[i+1]<=0:
                r=brentq(lambda s:s0*run(rc,s)[0],g[i],g[i+1],xtol=1e-11)
                print(rc,"Newtonian sigma2=",r, run(rc,r)[1:]); break
        else: print(rc,v)
