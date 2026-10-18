import numpy as np, json, sys, time
import background as bg; bg.mu_e=2.0
from fold_inhom import handoff, sc, rm
from lag_gr import Lag
N=int(sys.argv[1]); YH=float(sys.argv[2])
th,state=handoff(sc,yh=YH)
L=Lag(rm,N=N,pert=(state[0],state[1])); G,c=bg.G,bg.c
y=L.y0.copy(); t=0.0; rows=[]; t0=time.time(); k=0; last=-1; thr=[0.05,0.1,0.2,0.3,0.5,0.7,0.9,0.99]; hit={}
while True:
    r,U,m=L.unpack(y)
    comp=2*G*m[1:]/(r[1:]*c**2); mc=comp.max(); jm=np.argmax(comp)+1
    Gam2=1+U[1:]**2/c**2-2*G*m[1:]/(r[1:]*c**2)
    if k%2000==0 or mc>last*1.02:
        last=mc; Gam,rho,P,e=L.zone(r,U,m)
        jh=np.searchsorted(L.mu,0.5*L.mu_tot)
        rows.append([t,r[jh]/L.r0[jh],r[1]/L.r0[1],rho[0]/L.rho_c,mc,jm,m[jm]/m[-1],Gam2.min(),np.min(np.diff(r)),U.min()/c])
        if k%20000==0: print(k,'t=%.5f lam_h=%.4f lam_in=%.4f rhoc/rho0=%.3e maxcomp=%.4f at j=%d (m/M=%.3f) minGam2=%.2e minDr=%.2e Umin/c=%.3f wall=%.0f'%tuple([t]+rows[-1][1:4]+[mc,jm,m[jm]/m[-1],Gam2.min(),np.min(np.diff(r)),U.min()/c,time.time()-t0]),flush=True)
    for th_ in thr:
        if mc>=th_ and th_ not in hit: hit[th_]=(t,rows[-1]); print('comp>=',th_,'at t=%.5f'%t,'lam_h=%.4f'%rows[-1][1],flush=True)
    if mc>0.99 or Gam2.min()<1e-5 or np.min(np.diff(r))<=0 or t>40 or np.isnan(y).any(): break
    dt=L.dt_cfl(y)
    k1=L.f(t,y);k2=L.f(t+dt/2,y+dt/2*k1);k3=L.f(t+dt/2,y+dt/2*k2);k4=L.f(t+dt,y+dt*k3)
    y=y+dt/6*(k1+2*k2+2*k3+k4); t+=dt; k+=1
np.save(f'fold_BH_N{N}.npy',np.array(rows)); print('END t=%.5f steps=%d maxcomp=%.4f minGam2=%.2e minDr=%.2e'%(t,k,mc,Gam2.min(),np.min(np.diff(r))),flush=True)
