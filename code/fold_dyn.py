import numpy as np, json, sys, time
from scipy.integrate import solve_ivp
import background as bg
from homologous_calibrated import CalibratedHomologousModel
bg.mu_e=2.0
G,c=bg.G,bg.c
d=json.load(open('fold_step1.json')); rm=d['rho_max']
kap=18.733; a=0.004027       # exact-GR fold coefficients (fold_step1.json), a in Msun
yr=3.15576e7
m=CalibratedHomologousModel(rm); m.calibrate(target_k=0.0)
comp_max=np.max(2*G*m.Mr_bg[1:]/(m.r[1:]*c**2))
# smooth surrogate for F(lam) near lam=1 to beat round-off: F/eps^2 = c2 + c3 eps + c4 eps^2 + c5 eps^3 (fit on 0.004<|eps|<0.04)
ee=np.concatenate([np.linspace(-0.04,-0.004,25),np.linspace(0.004,0.04,25)])
Ff=np.array([m.lambda_ddot(1-e) for e in ee]); A=np.vander(ee,4,increasing=True)
coef=np.linalg.lstsq(A,Ff/ee**2,rcond=None)[0]; resid=np.max(np.abs(A@coef-Ff/ee**2)/np.abs(A@coef))
c2m=-coef[0]; print('model quadratic coefficient',c2m,'vs exact-GR 1.5*kappa',1.5*kap,'ratio',c2m/(1.5*kap),'fit resid',resid,flush=True)
EPSC=0.01
def F(lam):
    e=1-lam
    if abs(e)<EPSC: return e*e*(coef[0]+coef[1]*e+coef[2]*e*e+coef[3]*e**3)
    return m.lambda_ddot(lam)
def run_model(Mdot, S=30.0, c2=None, lam_end=0.15):
    # Mdot in Msun/yr; forcing keeps equilibrium on the exact-GR curve DeltaM = a*y^2
    md=Mdot/yr; beta=(4*a/(kap**2*md))**0.2
    cc=-coef[0] if c2 is None else c2
    fz=lambda t: cc*(-md*t)/(9*a)          # +cc*(Mmax-M)/(9a)
    t0=-S*beta; e0=-np.sqrt(max(fz(t0),0)/cc)  # lam>1 branch
    def rhs(t,y):
        return [y[1], F(y[0])+fz(t)]
    # slow-manifold start: equilibrium velocity ~ d/dt eps_eq; neglect (check S-convergence)
    ev15=lambda t,y: y[0]-lam_end; ev15.terminal=True
    ev5=lambda t,y: y[0]-0.5
    evy=lambda t,y: (y[0]**-3-1)-1e-2      # y = rho/rho_max-1 crosses 1e-2
    s=solve_ivp(rhs,(t0,t0+60*beta),[1-e0,0.0],method='DOP853',events=[ev5,ev15,evy],rtol=1e-11,atol=1e-14,max_step=beta/20)
    t5=s.t_events[0][0] if len(s.t_events[0]) else np.nan; t15=s.t_events[1][0] if len(s.t_events[1]) else np.nan
    ty=s.t_events[2][0] if len(s.t_events[2]) else np.nan
    return dict(Mdot=Mdot,beta=beta,alpha=2/(kap*beta**2),S=S,t_y1e2=ty,t_half=t5,t_015=t15,dt_plunge=t15-t5,dM_half=md*t5,t0=t0)
def run_NF(Mdot,S=30.0,uh=None):
    # normal form u''=u^2+s (dimensionless), start at equilibrium u=-sqrt(-s) at rest
    out={}
    f=lambda s,u:[u[1],u[0]**2+s]
    ev=lambda s,u:u[0]-200.0; ev.terminal=True
    sol=solve_ivp(f,(-S,40),[-np.sqrt(S),-0.5/np.sqrt(S)*0],method='DOP853',events=ev,rtol=1e-12,atol=1e-12,max_step=0.05)
    return sol
if __name__=="__main__":
    t0=time.time()
    res=[]
    for Md in (1e-9,1e-8,1e-7,1e-6,1e-5):
        r=run_model(Md); r['t_sec']=time.time()-t0; res.append(r); print(r,flush=True)
    json.dump(dict(res=res,c2m=c2m,kap=kap,a=a,comp_max=comp_max),open('fold_dyn.json','w'))
