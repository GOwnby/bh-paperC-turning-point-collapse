import numpy as np, background as bg
import metric, gr_rayleigh_ritz as rr
from scipy.linalg import eigh
from scipy.optimize import brentq
bg.mu_e=2.0; bg.Zlat=0.0
src=open('gr_rayleigh_ritz.py').read()
for use_eps in (True,False):
  bg.use_eps=use_eps
  for sgn in (1.0,0.0,-1.0):
    s2src=src.replace("if INCLUDE_8PI: Q = Q + 8*np.pi*G","Q = Q + SGN*8*np.pi*G")
    ns={'SGN':sgn,'__name__':'x'}
    exec(compile(s2src,'rrmod','exec'),ns)
    f=lambda rc: eigh(*ns['build_RR_matrices'](metric.build_metric(rc),5)[:2])[0][0]
    try:
        z=brentq(f,0.3e13,8e13,xtol=1e10)
    except Exception as e: z=str(e)
    print("use_eps",use_eps,"sign",sgn,"zero",z,flush=True)
