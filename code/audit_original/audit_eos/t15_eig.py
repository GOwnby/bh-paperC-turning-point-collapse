import sys; sys.path.insert(0,'orig')
import numpy as np
from scipy.linalg import eigh
import background as bg
from metric import build_metric
from gr_rayleigh_ritz import build_RR_matrices
bg.mu_e=2.0; bg.Zlat=0.0
d=build_metric(5e13)
for n in (2,3,5):
    A,B,_,_=build_RR_matrices(d,n); w,v=eigh(A,B); c0=v[:,0]
    print("n=%d sigma2=%.4f  xi/r(center)/xi/r(surface)=%.4f"%(n,w[0],c0.sum()/c0[0]))
# Newtonian-limit sanity: sigma^2 for n=2 vs two-mode
