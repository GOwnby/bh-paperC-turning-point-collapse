import sys
sys.path.insert(0,'/mnt/user-data/outputs/consolidated/paperC_core_instability')
import numpy as np, functools
import background as bg, tov_solve, metric
from scipy.linalg import eigh
import gr_rayleigh_ritz as rr
orig=tov_solve.solve_tov
for n in [3000,6000,12000,24000]:
    metric.solve_tov=functools.partial(orig,n_eval=n)
    b=metric.build_metric(5e13)
    out=[]
    for nb in [1,5]:
        A,B,_,_=rr.build_RR_matrices(b,nb)
        out.append(eigh(A,B)[0][0])
    print(n,"R km",b['R']/1e3,"sigma2(n=1,5)=",out,flush=True)
