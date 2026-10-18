import numpy as np, json, background as bg
from homologous_calibrated import CalibratedHomologousModel
bg.mu_e=2.0
d=json.load(open('fold_step1.json')); rm=d['rho_max']
m=CalibratedHomologousModel(rm)
fc,kN,kF=m.calibrate(target_k=0.0)   # exact turning-point theorem: sigma^2=0 at dM/drho_c=0
print('f_cal(k=0)',fc,flush=True)
kap=18.733
for eps in (-0.1,-0.03,-0.01,-0.003,-0.001,0.001,0.003,0.01,0.03,0.1,0.2,0.4):
    F=m.lambda_ddot(1-eps); print(eps,F, F/(-1.5*kap*eps**2),flush=True)
