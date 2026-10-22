import sys, os, numpy as np, io, contextlib
SRC='/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve'
os.chdir(SRC); sys.path.insert(0,SRC)
import background as bg; bg.mu_e=2.0
from lag_gr import Lag
G,c=bg.G,bg.c
rc=2.3106026e13
def acc(N,lam,lamdot):
    L=Lag(rc,N=N,pert=(lam,lamdot))
    a=L.f(0.0,L.y0,subtract=False,bracket=True); S=L.S
    r,U,m=L.unpack(L.y0); Gam,rho,P,e=L.zone(r,U,m)
    return L,a,S,rho,P
for lam,ld in ((1.0,0.0),(0.9687,-0.0247*0+(-2.0)),(0.9,-30.0)):
    print('state lam=%.4f lamdot=%.1f /s'%(lam,ld))
    Nref=3200; Lr,ar,Sr,rhor,Pr=acc(Nref,lam,ld)
    gref=G*Lr.M/Lr.R0**2
    for N in (100,200,400,800):
        L,a,S,rho,P=acc(N,lam,ld); f=Nref//N
        ia=np.arange(1,N+1)*f-1       # acc index j (interface j+1)
        ref=ar[ia]; d=a-ref
        w=np.abs(ref)+1e-30
        bulk=slice(0,int(0.9*N)); edge=slice(int(0.9*N),N)
        print('  N=%4d  max|acc-ref|/gsurf: bulk %.2e edge %.2e last node %.2e | rho(last zone) %.3e kg/m3'%(N,np.max(np.abs(d[bulk]))/gref,np.max(np.abs(d[edge]))/gref,abs(d[-1])/gref,rho[-1]))
    L,a,S,rho,P=acc(200,lam,ld)
    print('  N=200 static residual |S|/gsurf  max %.2e  at node %d ; last-node %.2e; |S|/|g_local| max %.2e'%(np.max(np.abs(S))/gref,np.argmax(np.abs(S))+1,abs(S[-1])/gref,np.max(np.abs(S)/np.abs(G*L.unpack(L.y0)[2][1:]/L.unpack(L.y0)[0][1:]**2))))
    print('  rho first/last zone', rho[0],rho[-1],'P last',P[-1])
