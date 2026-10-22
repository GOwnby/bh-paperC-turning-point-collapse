import numpy as np, background as bg
import metric, tov_solve
from scipy.linalg import eigh
import gr_rayleigh_ritz as rr
bg.mu_e=2.0; bg.Zlat=0.0
orig=tov_solve.solve_tov
for N in (3000,6000,12000):
    metric.solve_tov=lambda rc,N=N: orig(rc,n_eval=N)
    for rc in (5e13,2.31e13*1.0):
        d=metric.build_metric(rc)
        A,B,_,_=rr.build_RR_matrices(d,5); print(N,rc,eigh(A,B)[0][0],flush=True)
# sign flip
metric.solve_tov=lambda rc: orig(rc,n_eval=3000)
import types
src=open('gr_rayleigh_ritz.py').read()
for sgn in (+1,0,-1):
    rr.INCLUDE_8PI = (sgn!=0)
    # patch sign
    def mk(sgn):
        def Qf(d,n=5):
            r=d['r'];P=d['P'];a=d['a'];b=d['b'];dPdr=d['dPdr'];eps=d['eps']
            return None
    pass
