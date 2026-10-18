import numpy as np, json
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from twomode import TwoMode
import background as bg
G,c=bg.G,bg.c
T=TwoMode(); K=T.jac(); w,v=np.linalg.eig(K); i=np.argmax(w.real); vec=v[:,i].real; vec/=vec[0]
from gr_rayleigh_ritz import build_RR_matrices
from scipy.linalg import eigh
A,B,_,_=build_RR_matrices(T.m.bgd,2); s2n2=-eigh(A,B)[0][0]
out=dict(sigma2=w[i].real,second=w[1-i].real,mu_over_lam=vec[1],RR_n2=s2n2,f_cal=T.f_cal)
rhs=lambda t,y:[y[2],y[3],*T.acc(y[:2])] if T.valid(y[:2]) else [0,0,0,0]
for d0 in (1e-3,1e-6):
    s=solve_ivp(rhs,(0,6),np.concatenate([np.array([1.,0])-d0*vec,[0,0]]),method='Radau',max_step=0.02,rtol=1e-9,atol=1e-13,dense_output=True,events=[lambda t,y:y[0]-0.17])
    s.t_events
    lam=s.y[0]
    t05=brentq(lambda t:s.sol(t)[0]-0.5,0,s.t[-1])
    out[str(d0)]=dict(t_half=t05,t_end=s.t[-1],lam_end=lam[-1],dt_plunge=s.t[-1]-t05, mu_lam_05=s.sol(t05)[1]/0.5, status=int(s.status))
    print(d0,out[str(d0)],flush=True)
json.dump(out,open('twomode_final.json','w'),default=float)
