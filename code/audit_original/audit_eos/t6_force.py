import sys; sys.path.insert(0,'orig')
import numpy as np
import background as bg
from homologous_calibrated import CalibratedHomologousModel
from exactV import Star
G,c=bg.G,bg.c
bg.mu_e=2.0; bg.Zlat=0.0; bg.use_eps=True
st=Star(5e13)
m=CalibratedHomologousModel(5e13); fc,kN,kF=m.calibrate_to_RR()
I=m.I
print("f_cal",fc)
print(" lam    F_exact/I  [s^-2]   code(f_cal)   code(f=1,R0 sub)  code(f=0 Newt)   ratio code/exact")
F0=None
for lam in (1.01,0.99,0.98,0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2):
    h=min(2e-3,0.05*(1-lam) if lam<0.99 else 2e-3)
    h=max(h,1e-3) if lam<0.99 else 1e-3
    Fex=-(st.M_of_lam(lam+h)-st.M_of_lam(lam-h))/(2*h)*c**2/I
    cal=m.lambda_ddot(lam,fc); f1=m.lambda_ddot(lam,1.0); f0=m.lambda_ddot(lam,0.0)
    print("%5.2f %14.5e %14.5e %14.5e %14.5e   %.4f"%(lam,Fex,cal,f1,f0,cal/Fex))
