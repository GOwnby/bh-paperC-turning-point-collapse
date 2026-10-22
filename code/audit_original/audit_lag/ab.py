import sys, numpy as np, time
sys.path.insert(0,'/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve')
import background as bg; bg.mu_e=2.0
from lag_gr import Lag
from lag_jac import jac_eigs
def run(rc,N):
    t=time.time(); L=Lag(rc,N=N,pert=None); ev,V=jac_eigs(L); return L,ev,V,time.time()-t
mode=sys.argv[1]; Ns=[int(x) for x in sys.argv[2:]]
for N in Ns:
    if mode=='a':
        L,ev,V,dt=run(5e13,N)
        print('a N=%d max Re=%.6f  (#Re>1e-3: %d) time %.0fs'%(N,ev.real.max(),(ev.real>1e-3).sum(),dt),flush=True)
        print('   max|Re| of rest',np.sort(np.abs(ev.real))[-3:],flush=True)
    elif mode=='b':
        rc=2.3106026e13
        L,ev,V,dt=run(rc,N)
        s2=-(ev**2)   # sigma^2 = -lambda^2
        sr=s2.real
        idx=np.argsort(np.abs(ev))
        print('b N=%d  smallest |lambda| list:'%N,flush=True)
        for i in idx[:8]: print('    lam=%.5g%+.5gj  sigma2=%.5f'%(ev[i].real,ev[i].imag,s2[i].real))
        nz=[i for i in idx if abs(ev[i])>1e-3]
        print('   lowest nonzero sigma2 = %.5f'%s2[nz[0]].real, 'time %.0fs'%dt,flush=True)
    elif mode=='c':  # stable branch
        L,ev,V,dt=run(1e13,N)
        w=np.abs(ev.imag); w=np.sort(w[w>1e-3])
        print('c N=%d lowest freq %.5f  maxRe=%.2e'%(N,w[0],np.abs(ev.real).max()),flush=True)
