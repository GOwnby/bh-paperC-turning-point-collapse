# Independent derivation: linearised radial Euler eq + Einstein constraint (X, DeltaP) system
import sys; sys.path.insert(0,'.')
from pulse import *
def rhsB(w2):
    def f(r,y):
        m,nu,x,X,dP=y
        d=bg_rhs(r,y[:3])
        e=egeo(x); Pg=pgeo(x); G1=Gamma1(x)
        e2l=1/(1-2*m/r)
        n1=nu1(r,m,Pg)
        lam1=e2l*(4*np.pi*r*e-m/r**2)
        Pp=-(e+Pg)*n1
        ep=(e+Pg)/(G1*Pg)*Pp
        N=m+4*np.pi*r**3*Pg; D=r*(r-2*m)
        Np=4*np.pi*r*r*e+12*np.pi*r*r*Pg+4*np.pi*r**3*Pp
        Dp=2*r-2*m-8*np.pi*r**3*e
        n2=(Np*D-N*Dp)/D**2
        Ppp=-(ep+Pp)*n1-(e+Pg)*n2
        Xp=-dP/(G1*Pg)-(2/r-n1)*X   # Delta n/n = -(X+2X/r-nu X) incl. perturbed g_rr (delta lambda)
        dm=-4*np.pi*r*r*(e+Pg)*X
        dPe=dP-Pp*X
        dn1=e2l*(2*n1*dm/r+dm/r**2+4*np.pi*r*dPe)
        dPp=Ppp*X+Pp*Xp+w2*(e+Pg)*e2l*np.exp(-2*nu)*X-n1*((e+Pg)/(G1*Pg)+1)*dP+n1*(ep+Pp)*X-(e+Pg)*dn1
        return d+[Xp,dPp]
    return f
def shootB(rho_c,w2,x_end=1e-3,r0=1e-6,rtol=1e-12,atol=1e-14):
    xc=x_of_rho(rho_c); Pg=pgeo(xc); G1=Gamma1(xc); e=egeo(xc)
    m0=4*np.pi/3*e*r0**3
    x0=xc-(2*np.pi/3)*(e+Pg)*(e+3*Pg)*r0**2/dpdx_geo(xc)
    nu0=(2*np.pi/3)*(e+3*Pg)*r0**2
    ev=lambda r,y: y[2]-x_end
    ev.terminal=True; ev.direction=-1
    sol=solve_ivp(rhsB(w2),[r0,50.0],[m0,nu0,x0,r0,-3*G1*Pg],method='DOP853',rtol=rtol,atol=atol,events=ev)
    yR=sol.y_events[0][0]
    return yR[4]/(3*G1*Pg)
