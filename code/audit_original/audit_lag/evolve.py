# independent driver: imports lag_gr / fold_inhom (unmodified), integrates to apparent horizon with crossing interpolation
import sys, os, json, time, numpy as np
SRC='/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve'
os.chdir(SRC); sys.path.insert(0,SRC)
import background as bg; bg.mu_e=2.0
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from fold_inhom import handoff, sc, rm
from lag_gr import Lag
N=int(sys.argv[1]); cfl=float(sys.argv[2]); YH=float(sys.argv[3]); stop=sys.argv[4] if len(sys.argv)>4 else 'AH'
tag=sys.argv[5] if len(sys.argv)>5 else ''
G,c=bg.G,bg.c
th,state=handoff(sc,yh=YH)
L=Lag(rm,N=N,pert=(state[0],state[1]))
mu=L.mu; mtot=L.mu_tot
jh=np.searchsorted(mu,0.5*mtot)
def lam_mass(r):   # r at exact half rest mass, interpolated; compared to initial
    return np.interp(0.5*mtot,mu,r)/np.interp(0.5*mtot,mu,L.r0)
def diag(t,y):
    r,U,m=L.unpack(y); Gam,rho,P,e=L.zone(r,U,m)
    V=4*np.pi/3*(r[1:]**3-r[:-1]**3)
    cum=m[1:]-np.cumsum(e*V/c**2)
    comp=2*G*m[1:]/(r[1:]*c**2); jm=int(np.argmax(comp))
    Gam2=1+U[1:]**2/c**2-2*G*m[1:]/(r[1:]*c**2)
    return dict(t=t,mN=m[-1],cum_max_over_M=float(np.max(np.abs(cum))/m[-1]),cum_at_N=float(cum[-1]/m[-1]),
                zone_max_over_M=float(np.max(np.abs(np.diff(m)-e*V/c**2))/m[-1]),
                comp=float(comp.max()),mAH_over_M=float(m[jm+1]/m[-1]),jm=jm+1,minGam2=float(Gam2.min()),minDr=float(np.diff(r).min()),
                cum_at_jm_rel=float(cum[jm]/m[jm+1]),rhoc=float(rho[0]/L.rho_c),rN=float(r[-1]/L.R0),lam_h=float(r[jh]/L.r0[jh]),lam_mass=float(lam_mass(r)),
                alpha_c_over_s=float((1/(((e+P)/rho)[0]/c**2))))
y=L.y0.copy(); t=0.0; k=0; cross={}; thr=[0.5,0.15]; out=[]; nxt=0.0
lam_prev=lam_mass(L.unpack(y)[0]); lamh_prev=L.unpack(y)[0][jh]/L.r0[jh]; comp_prev=0; t0=time.time()
out.append(diag(t,y)); mN0=out[0]['mN']
tmax=float(stop) if stop!='AH' else 40
while True:
    dt=L.dt_cfl(y,cfl)
    k1=L.f(t,y);k2=L.f(t+dt/2,y+dt/2*k1);k3=L.f(t+dt/2,y+dt/2*k2);k4=L.f(t+dt,y+dt*k3)
    yn=y+dt/6*(k1+2*k2+2*k3+k4); tn=t+dt; k+=1
    r,U,m=L.unpack(yn)
    lam=lam_mass(r); lamh=r[jh]/L.r0[jh]
    for v in thr:
        if ('m%g'%v) not in cross and lam<v<=lam_prev: cross['m%g'%v]=t+dt*(lam_prev-v)/(lam_prev-lam)
        if ('h%g'%v) not in cross and lamh<v<=lamh_prev: cross['h%g'%v]=t+dt*(lamh_prev-v)/(lamh_prev-lamh)
    comp=np.max(2*G*m[1:]/(r[1:]*c**2))
    if comp>=1.0 and 'AH' not in cross:
        cross['AH']=t+dt*(1-comp_prev)/(comp-comp_prev); d=diag(tn,yn); cross['AHdiag']=d
        if stop=='AH': y=yn; t=tn; break
    if t>=nxt:
        out.append(diag(t,y)); nxt+= (0.1 if lam>0.7 else 0.002)
    lam_prev=lam; lamh_prev=lamh; comp_prev=comp; y=yn; t=tn
    if t>=tmax or np.isnan(y).any(): break
    if stop!='AH' and False: pass
d=diag(t,y); out.append(d)
res=dict(N=N,cfl=cfl,YH=YH,steps=k,wall=time.time()-t0,cross=cross,final=d,series=out,mN0=mN0,mN_end=float(L.unpack(y)[2][-1]))
json.dump(res,open('/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/audit_lag/ev_N%d_cfl%g_y%g%s.json'%(N,cfl,YH,tag),'w'),default=float)
print('N=%d cfl=%g steps=%d wall=%.0f'%(N,cfl,k,time.time()-t0)); print(' cross',{k_:(v if k_!='AHdiag' else '') for k_,v in cross.items()})
print(' final',d)
