import sys
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np, functools
import background as bg, tov_solve, metric
from scipy.linalg import eigh
from scipy.optimize import brentq
import gr_rayleigh_ritz as rr
def s2(rc,nb=5,n=3000):
    metric.solve_tov=functools.partial(tov_solve.solve_tov,n_eval=n)
    b=metric.build_metric(rc)
    A,B,_,_=rr.build_RR_matrices(b,nb)
    return eigh(A,B)[0][0], b['M']/bg.M_sun
for rc in [1e13,2e13,2.31e13,3e13,5e13]:
    print(rc, "proj(3000pts,5fn)=",s2(rc), "proj(12000pts)=",s2(rc,n=12000)[0],flush=True)
f=lambda lr: s2(10**lr)[0]
print("proj zero crossing (3000):",10**brentq(f,13.30,13.45,xtol=1e-9))
