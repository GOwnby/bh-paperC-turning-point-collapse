# Exact-GR energy functional for homologous deformations of a TOV star at fixed rest-mass labelling.
# V(lam) = M_grav(lam) c^2, where M(lam,.) solves the Hamiltonian constraint for the configuration r_f = lam*r
# with rest-mass conservation dm_rest = 4 pi r_f^2 rho (1-2GM/(r_f c^2))^(-1/2) dr_f.
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

def run_exact_model(rho_c=5e13, d0=1e-6, n_lam=140):
    """Calibration-free cross-check: lam'' = -V'(lam)/I_rel with V(lam)=M_grav(lam)c^2 exact (Hamiltonian constraint,
    rest-mass-preserving homologous deformation) and the relativistic inertia I_rel (Chandrasekhar weight) held at its lam=1 value."""
    from scipy.interpolate import CubicSpline
    from scipy.integrate import solve_ivp
    from metric import build_metric
    from gr_rayleigh_ritz import build_RR_matrices
    from scipy.linalg import eigh
    st=Star(rho_c); d=build_metric(rho_c)
    r,a,b,eps,P=d['r'],d['a'],d['b'],d['eps'],d['P']
    I_rel=np.trapezoid(4*np.pi*np.exp(-a+3*b)*(eps+P/c**2)*r**4,r)
    h=2e-3
    Vpp=lambda h:(st.M_of_lam(1+h)-2*st.M_of_lam(1.0)+st.M_of_lam(1-h))/h**2*c**2
    Vpp_R=(4*Vpp(h)-Vpp(2*h))/3
    A,B,_,_=build_RR_matrices(d,1)
    out=dict(k_exact=-Vpp_R/I_rel, k_RR_n1=-eigh(A,B)[0][0], I_rel=I_rel)
    lg=np.linspace(np.log(0.08),np.log(1.06),n_lam)
    V=np.array([st.M_of_lam(np.exp(l)) for l in lg])*c**2
    spl=CubicSpline(lg,V,bc_type='natural')
    F=lambda lam:-spl(np.log(lam),1)/lam
    F1=F(1.0)
    ev=lambda t,y:y[0]-0.15; ev.terminal=True
    s=solve_ivp(lambda t,y:[y[1],(F(y[0])-F1)/I_rel],(0,40/np.sqrt(out['k_exact'])),[1-d0,0],method='Radau',events=ev,max_step=0.02,rtol=1e-9,atol=1e-13,dense_output=True)
    lam=s.y[0]; th=s.t[np.argmax(lam<0.5)]
    # dense-output crossing times
    from scipy.optimize import brentq
    t05=brentq(lambda t:s.sol(t)[0]-0.5,th-0.05,th+0.05); out.update(t_half=t05, dt_plunge=s.t[-1]-t05)
    return out

if __name__=="__main__":
    bg.mu_e=2.0; bg.Zlat=0.0; bg.use_eps=True
    print(run_exact_model())
