import sys
sys.dont_write_bytecode=True
sys.path.insert(0,'/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve')
import numpy as np, json, background as bg
from homologous_calibrated import CalibratedHomologousModel
bg.mu_e=2.0
rm=23106017651.354734e3
m=CalibratedHomologousModel(rm)
fc,kN,kF=m.calibrate(target_k=0.0); print('f_cal',fc,flush=True)
kap=18.7335
for eps in (-0.003,-0.001,0.001,0.003,0.01):
    F=m.lambda_ddot(1-eps); print(eps,F,F/(-1.5*kap*eps**2),flush=True)
