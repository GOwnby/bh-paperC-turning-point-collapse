import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
import background as bg
from metric import build_metric
from gr_rayleigh_ritz import build_RR_matrices, make_basis
from tov_solve import rho_of_P
G,c=bg.G,bg.c
bg.mu_e=2.0; bg.Zlat=0.0; bg.use_eps=True
out='/tmp/'
# ---- Newtonian structure + Newtonian RR (c->inf limit of same variational form) ----
def newtonian_s2(rc,nb=5,n=4000):
    Pc=bg.P_of_rho(rc)
    def rhs(r,y):
        P,M=y
        if P<=1e-10*Pc: return [0,0]
        rho=rho_of_P(P,rc); return [-G*rho*M/r**2, 4*np.pi*r**2*rho]
    ev=lambda r,y:y[0]-1e-10*Pc; ev.terminal=True
    r0=1.0; s=solve_ivp(rhs,[r0,5e7],[Pc,4/3*np.pi*r0**3*rc],events=ev,rtol=1e-10,atol=[Pc*1e-14,1e-3],dense_output=True,max_step=2e4)
    R=s.t_events[0][0]; r=np.linspace(r0,R,n); P=np.clip(s.sol(r)[0],0,None); M=s.sol(r)[1]
    rho=np.array([rho_of_P(p,rc) for p in P]); dPdr=-G*rho*M/r**2
    Gam=np.array([bg.Gamma1_of_rho(x) if x>1 else 5/3 for x in rho])
    fns=make_basis(R,nb); xis=[f(r) for f in fns]; zs=[r**2*x for x in xis]; dz=[np.gradient(z,r) for z in zs]
    Pi=Gam*P/r**2; Q=4*dPdr*r; S=rho*r**2
    A=np.array([[np.trapezoid(Pi*dz[i]*dz[j]+Q*xis[i]*xis[j],r) for j in range(nb)] for i in range(nb)])
    B=np.array([[np.trapezoid(S*xis[i]*xis[j],r) for j in range(nb)] for i in range(nb)])
    return eigh(A,B)[0][0], 2*G*M[-1]/(R*c**2)
