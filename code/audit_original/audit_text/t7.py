import numpy as np, background as bg
import metric, tov_solve
from scipy.linalg import eigh
import gr_rayleigh_ritz as rr
bg.mu_e=2.0; bg.Zlat=0.0
orig=tov_solve.solve_tov
rp=2.310519549261331e13
for N in (3000,12000,24000,48000):
    metric.solve_tov=lambda rc,N=N: orig(rc,n_eval=N)
    d=metric.build_metric(rp)
    out=[]
    for nb in (1,2,3,5,7):
        A,B,_,_=rr.build_RR_matrices(d,nb); out.append(round(eigh(A,B)[0][0],4))
    print(N,out,flush=True)
