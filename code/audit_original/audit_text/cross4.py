import numpy as np, background as bg
import gr_rayleigh_ritz as rr
from metric import build_metric
from scipy.linalg import eigh
from scipy.optimize import brentq
bg.mu_e=2.0; bg.Zlat=0.0
for flag in (False,True):
  bg.use_eps=flag
  s2=lambda rc: eigh(*rr.build_RR_matrices(build_metric(rc),5)[:2])[0][0]
  print("use_eps",flag,[ (rc, round(s2(rc),2)) for rc in (1e13,2e13,3e13,4e13,5e13)],flush=True)
  try:
    z=brentq(s2,1e13,6e13,xtol=1e10); print("  sigma^2=0 at",z)
  except Exception as e: print("  no root",e)
  from scipy.optimize import minimize_scalar
  r=minimize_scalar(lambda lr:-build_metric(10**lr)['M'],bounds=(13.0,14.2),method='bounded',options={'xatol':1e-3}); print("  M peak at",10**r.x, -r.fun/bg.M_sun)
