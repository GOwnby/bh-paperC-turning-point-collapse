import numpy as np, background as bg
from metric import build_metric
from scipy.optimize import minimize_scalar
FMT = {'He4':(2,2,1.56e10,1.40906),'C12':(6,2,2.12e10,1.38603),'O16':(8,2,1.94e10,1.38024),'Fe56':(26,56/26,1.18e9,1.10618)}
for name,(Z,mue,rf,mf) in FMT.items():
    bg.Zlat=float(Z); bg.mu_e=mue
    res=minimize_scalar(lambda lr:-build_metric(10**lr)['M'],bounds=(12.8,14.2),method='bounded',options={'xatol':1e-3})
    rc=10**res.x; print(f"{name}: Salpeter-GR rho_max={rc/1e3:.3e} g/cc  Mmax={-res.fun/bg.M_sun:.5f} | FMT-GR (Rotondo) rho={rf:.3e} M={mf}",flush=True)
