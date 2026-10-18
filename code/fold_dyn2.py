import numpy as np, json
from scipy.integrate import solve_ivp
import fold_dyn as fd
from fold_dyn import *
def run2(Mdot,S=30.0,scale=1.0,lam_end=0.15,sm=True):
    md=Mdot/yr; beta=(4*a/(kap**2*md))**0.2
    cc=-coef[0]*scale
    fz=lambda t: cc*(-md*t)/(9*a)
    t0=-S*beta; kk=md/(9*a); e0=-np.sqrt(-kk*t0)
    edot0=0.5*kk/np.sqrt(-kk*t0) if sm else 0.0
    Fs=lambda lam: scale*F(lam)
    rhs=lambda t,y:[y[1],Fs(y[0])+fz(t)]
    ev15=lambda t,y:y[0]-lam_end; ev15.terminal=True
    ev5=lambda t,y:y[0]-0.5
    evy=lambda t,y:(y[0]**-3-1)-1e-2
    s=solve_ivp(rhs,(t0,t0+60*beta),[1-e0,-edot0],method='DOP853',events=[ev5,ev15,evy],rtol=1e-11,atol=1e-14,max_step=beta/20)
    ty=s.t_events[2][0]; t5=s.t_events[0][0]; t15=s.t_events[1][0]
    return dict(Mdot=Mdot,S=S,scale=scale,sm=sm,beta=beta,t_half=t5,t_half_over_beta=t5/beta,t_y1e2=ty,dt_plunge=t15-t5,dM_half=md*t5)
out=[]
for kw in (dict(S=20),dict(S=30),dict(S=45),dict(S=30,sm=False),dict(S=30,scale=1.5*kap/c2m),dict(S=45,scale=1.5*kap/c2m)):
    r=run2(1e-7,**kw); out.append(r); print(r,flush=True)
for Md in (1e-9,1e-8,1e-6,1e-5):
    for sc in (1.0,1.5*kap/c2m):
        r=run2(Md,S=45,scale=sc); out.append(r); print(r,flush=True)
# universal PI solution
f=lambda s,u:[u[1],u[0]**2+s]
ev=lambda s,u:u[0]-1e6; ev.terminal=True
pis=[]
for S in (20,30,45,60):
    sol=solve_ivp(f,(-S,60),[-np.sqrt(S),-0.5/np.sqrt(S)],method='DOP853',events=ev,rtol=1e-12,atol=1e-12,max_step=0.02)
    # pole: u ~ 6/(s-s*)^2  -> s* = s_end + sqrt(6/u)
    s_end=sol.t[-1]; u_end=sol.y[0,-1]; pis.append((S,s_end+np.sqrt(6/u_end))); print('PI',S,pis[-1],flush=True)
json.dump(dict(out=out,PI=pis),open('fold_dyn2.json','w'))
