import sys, numpy as np, time
sys.path.insert(0,'/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve')
import background as bg; bg.mu_e=2.0
from fold_inhom import handoff, sc, rm
from lag_gr import Lag
N=int(sys.argv[1]); G,c=bg.G,bg.c
th,state=handoff(sc,yh=0.1)
L=Lag(rm,N=N,pert=(state[0],state[1]))
y=L.y0.copy(); t=0.0; k=0; t0=time.time()
Msun=bg.M_sun
prev=None; times={}
while True:
    r,U,m=L.unpack(y)
    comp=2*G*m[1:]/(r[1:]*c**2); mc=comp.max()
    jh=np.searchsorted(L.mu,0.5*L.mu_tot); lam=r[jh]/L.r0[jh]
    if prev is not None:
        for v in (0.5,0.15):
            if v not in times and lam<v: times[v]=np.interp(v,[lam,prev[1]],[t,prev[0]])
    if mc>=1.0:
        # interpolate time for comp=1
        tAH=np.interp(1.0,[prev[2],mc],[prev[0],t]); 
        break
    prev=(t,lam,mc,y.copy())
    if mc>0.97 and 'f' not in times: times['f']=1
    dt=L.dt_cfl(y)
    k1=L.f(t,y);k2=L.f(t+dt/2,y+dt/2*k1);k3=L.f(t+dt/2,y+dt/2*k2);k4=L.f(t+dt,y+dt*k3)
    y=y+dt/6*(k1+2*k2+2*k3+k4); t+=dt; k+=1
r,U,m=L.unpack(y)
Gam,rho,P,e=L.zone(r,U,m)
comp=2*G*m[1:]/(r[1:]*c**2); j=np.argmax(comp)+1
print('N',N,'t(0.5)=%.4f t(0.15)=%.4f t_AH(comp=1, interp)=%.4f'%(times[0.5],times[0.15],tAH),'steps',k,'wall',time.time()-t0)
print('maxcomp',mc,'at node',j,'m/M=%.4f'%(m[j]/m[-1]),'m_AH/Msun=%.4f'%(m[j]/Msun),'M/Msun=%.4f'%(m[-1]/Msun),'rho_c/rho_c(handoff lattice)=%.3e'%(rho[0]/L.rho_c),'rho_c/rho_max=%.3e'%(rho[0]/rm))
print('central density abs [g/cc]: %.3e'%(rho[0]*1e-3), 'central rho/rho_c0',rho[0]/L.rho_c)
# initial (handoff) central density in the code
r0,U0,m0=L.unpack(L.y0); Gam0,rho0,P0,e0=L.zone(r0,U0,m0)
print('initial central rho / rho_max',rho0[0]/rm,' => rise factor from handoff state',rho[0]/rho0[0])
fr=[0.05,0.2,0.5,0.8,0.95,0.99]
mu=L.mu; 
out=[]
for f in fr:
    jj=np.searchsorted(mu,f*L.mu_tot); 
    # interpolate
    ra=np.interp(f*L.mu_tot,mu,r/L.r0.clip(1e-30)*(L.r0>0)+0)
    out.append(ra)
print('radius ratios at AH for rest-mass fractions',fr,np.round(out,4),'surface',r[-1]/L.r0[-1])
