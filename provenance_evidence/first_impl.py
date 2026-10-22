import numpy as np, background as bg, metric, gr_rayleigh_ritz as rr
from scipy.linalg import eigh
bg.mu_e=2.0; bg.Zlat=0.0
for use_eps in (False,):
    bg.use_eps=use_eps
    rr.INCLUDE_8PI=False
    d=metric.build_metric(5e13)
    print('use_eps',use_eps,'8PI off: sigma2(n=1..5) =',[round(eigh(*rr.build_RR_matrices(d,n)[:2])[0][0],4) for n in range(1,6)], 'M',d['M']/bg.M_sun)
