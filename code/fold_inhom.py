import numpy as np, json, sys, time
from scipy.integrate import solve_ivp
import background as bg; bg.mu_e=2.0
from fold_dyn import F, kap, a, yr, coef
from lag_gr import Lag
d=json.load(open('fold_step1.json')); rm=d['rho_max']
c2m=-coef[0]; sc=1.5*kap/c2m
def handoff(scale, Mdot=1e-7, S=30.0, yh=1e-2):
    md=Mdot/yr; beta=(4*a/(kap**2*md))**0.2; cc=c2m*scale; kk=md/(9*a)
    fz=lambda t: cc*(-md*t)/(9*a); t0=-S*beta; e0=-np.sqrt(-kk*t0); ed=0.5*kk/np.sqrt(-kk*t0)
    rhs=lambda t,y:[y[1],scale*F(y[0])+fz(t)]
    ev=lambda t,y:(y[0]**-3-1)-yh; ev.terminal=True
    s=solve_ivp(rhs,(t0,t0+60*beta),[1-e0,-ed],method='DOP853',events=ev,rtol=1e-11,atol=1e-14,max_step=beta/20)
    return s.t_events[0][0], s.y_events[0][0]
def model_from(state, scale, lam_end=0.15):
    ev5=lambda t,y:y[0]-0.5; ev15=lambda t,y:y[0]-lam_end; ev15.terminal=True
    s=solve_ivp(lambda t,y:[y[1],scale*F(y[0])],(0,60),list(state),method='DOP853',events=[ev5,ev15],rtol=1e-11,atol=1e-14,max_step=0.02)
    return s.t_events[0][0], s.t_events[1][0]
if __name__=="__main__":
    N=int(sys.argv[1]); 
    YH=float(sys.argv[2]); th,state=handoff(sc,yh=YH); print('handoff t',th,'lam,lamdot',state,flush=True)
    for scl in (sc,1.0):
        t5,t15=model_from(state,scl); print('model scale %.3f: t(0.5)=%.4f t(0.15)=%.4f plunge %.4f'%(scl,t5,t15,t15-t5),flush=True)
    L=Lag(rm,N=N,pert=(state[0],state[1]))
    y=L.y0.copy(); t=0.0; rows=[]; nxt=0; t0=time.time(); k=0
    while True:
        r,U,m=L.unpack(y); jh=np.searchsorted(L.mu,0.5*L.mu_tot); lam=r[jh]/L.r0[jh]
        if t>=nxt:
            Gam,rho,P,e=L.zone(r,U,m)
            rows.append([t,lam,rho[0]/L.rho_c,np.max(2*bg.G*m[1:]/(r[1:]*bg.c**2)),np.max(np.abs(np.diff(m)-e*4*np.pi/3*(r[1:]**3-r[:-1]**3)/bg.c**2))/m[-1]]); nxt+=(0.05 if lam>0.7 else 0.0005)
            if len(rows)%10==0: print(k,'t=%.3f lam_h=%.5f rho_c/rho_c0=%.4f maxcomp=%.5f viol=%.1e wall=%.0f'%(t,*rows[-1][1:],time.time()-t0),flush=True)
        if lam<0.12 or t>40: break
        dt=L.dt_cfl(y)
        k1=L.f(t,y);k2=L.f(t+dt/2,y+dt/2*k1);k3=L.f(t+dt/2,y+dt/2*k2);k4=L.f(t+dt,y+dt*k3)
        y=y+dt/6*(k1+2*k2+2*k3+k4); t+=dt; k+=1
        rows_fine=None
    # fine record of crossing times: re-interpolate from rows
    rows=np.array(rows); np.save(f'fold_inhom_N{N}_y{YH}.npy',rows)
    lam=rows[:,1]; tt=rows[:,0]
    def cross(v):
        i=np.argmax(lam<v); 
        return np.interp(v,[lam[i],lam[i-1]],[tt[i],tt[i-1]]) if i>0 else np.nan
    print('Lag N=%d: t(0.5)=%.4f t(0.15)=%.4f'%(N,cross(0.5),cross(0.15)),flush=True)
