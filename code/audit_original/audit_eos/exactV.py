# Exact-GR energy functional for homologous deformations of a TOV star at fixed rest-mass labelling.
# V(lam) = M_grav(lam) c^2, where M(lam,.) solves the Hamiltonian constraint for the configuration r_f = lam*r
# with rest-mass conservation dm_rest = 4 pi r_f^2 rho (1-2GM/(r_f c^2))^(-1/2) dr_f.
import sys; sys.path.insert(0,'orig')
import numpy as np
from scipy.integrate import solve_ivp
import background as bg
G,c=bg.G,bg.c
class Star:
    def __init__(self, rhoc, Pcut=1e-10, rtol=1e-13):
        self.rhoc=rhoc; xc=bg.x_of_rho(rhoc); Pc=bg.P_of_x(xc)
        def rhs(r,y):
            x,m=y
            rho=bg.rho_of_x(x); eps=bg.eps_over_c2_of_rho(rho); P=bg.P_of_x(x)
            dP=-G*(eps+P/c**2)*(m+4*np.pi*r**3*P/c**2)/(r**2*(1-2*G*m/(r*c**2)))
            dPdx=bg.A_P*bg.dfdx(x)-4*bg.B_lat()*x**3
            return [dP/dPdx,4*np.pi*r*r*eps]
        ev=lambda r,y: bg.P_of_x(y[0])-Pcut*Pc; ev.terminal=True
        r0=1.0
        s=solve_ivp(rhs,[r0,5e7],[xc,4/3*np.pi*r0**3*bg.eps_over_c2_of_rho(rhoc)],events=ev,rtol=rtol,atol=[1e-15,1e-8],method='DOP853',dense_output=True)
        self.R=s.t_events[0][0]; self.M=s.y_events[0][0][1]; self.sol=s; self.Pc=Pc; self.Pcut=Pcut
    def bgfun(self,r):
        x,m=self.sol.sol(r)
        rho=bg.rho_of_x(x)
        return rho,m,bg.P_of_x(x)
    def M_of_lam(self, lam, rtol=1e-13, full=False):
        """total gravitational mass of the homologous deformation r_f=lam*r (rest mass of each shell preserved)."""
        def rhs(r,y):
            M=y[0]
            rho0,m0,_=self.bgfun(r)
            eb0=1/np.sqrt(1-2*G*m0/(r*c**2))
            eb=1/np.sqrt(1-2*G*M/(lam*r*c**2))
            rho=rho0*eb0/(lam**3*eb)
            eps=bg.eps_over_c2_of_rho(rho)
            return [4*np.pi*lam**3*r*r*eps]
        r0=1.0
        rho0,m0,_=self.bgfun(r0)
        s=solve_ivp(rhs,[r0,self.R],[4/3*np.pi*lam**3*r0**3*bg.eps_over_c2_of_rho(rho0/lam**3)],rtol=rtol,atol=1e-10,method='DOP853',dense_output=full)
        return (s.y[0,-1],s) if full else s.y[0,-1]
if __name__=="__main__":
    bg.mu_e=2.0; bg.Zlat=0.0; bg.use_eps=True
    st=Star(5e13)
    print("R=%.4f km M=%.6e kg"%(st.R/1e3,st.M))
    print("M(lam=1) from constraint ODE = %.8e   rel diff to TOV mass %.2e"%(st.M_of_lam(1.0),st.M_of_lam(1.0)/st.M-1))
    for h in (4e-3,2e-3,1e-3,5e-4):
        Vp=(st.M_of_lam(1+h)-st.M_of_lam(1-h))/(2*h)*c**2
        Vpp=(st.M_of_lam(1+h)-2*st.M_of_lam(1.0)+st.M_of_lam(1-h))/h**2*c**2
        print(" h=%g  V'(1)=%.4e J  V''(1)=%.6e J   (Mc^2=%.3e)"%(h,Vp,Vpp,st.M*c**2))
