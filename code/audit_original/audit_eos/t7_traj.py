import sys; sys.path.insert(0,'orig')
import numpy as np, time, pickle
import background as bg
from homologous_calibrated import CalibratedHomologousModel
from exactV import Star
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp
G,c=bg.G,bg.c
bg.mu_e=2.0; bg.Zlat=0.0; bg.use_eps=True
st=Star(5e13)
m=CalibratedHomologousModel(5e13); fc,kN,kF=m.calibrate_to_RR()
# tabulate V(lam) exactly on a grid (in x = ln lam for smoothness)
t0=time.time()
lg=np.linspace(np.log(0.08),np.log(1.06),160)
V=np.array([st.M_of_lam(np.exp(l),rtol=1e-13) for l in lg])*c**2
print("tabulated V in %.1fs"%(time.time()-t0))
Vspl=CubicSpline(lg,V,bc_type='natural')
# F = -dV/dlam = -(dV/dlnlam)/lam
def F_exact(lam): return -Vspl(np.log(lam),1)/lam
# check vs finite diff
print("F_exact/I_code at 0.5:",F_exact(0.5)/m.I," at 0.99:",F_exact(0.99)/m.I)
# static eq shift: subtract F(1) (spline noise) -> impose F(1)=0 by subtracting constant
F1=F_exact(1.0); print("F_exact(1)/I_code (spline residual) = %.3e s^-2"%(F1/m.I))
I_rel=2.59565e41
def run(accel,d0=1e-6,smax=40/np.sqrt(27.81)):
    ev=lambda t,y:y[0]-0.15; ev.terminal=True
    s=solve_ivp(lambda t,y:[y[1],accel(y[0])],(0,smax),[1-d0,0],method='Radau',events=ev,max_step=0.02,rtol=1e-9,atol=1e-13)
    lam=s.y[0]; th=s.t[np.argmax(lam<0.5)]
    # time from 0.5 to 0.15 using interpolation
    return th, s.t[-1]-th, s
res={}
res['code(f_cal), I_code']=run(lambda l:m.lambda_ddot(l))
res['exact V, I_rel']=run(lambda l:(F_exact(l)-F1)/I_rel)
res['exact V, I_code']=run(lambda l:(F_exact(l)-F1)/m.I)
res['code f=1 (uncal), I_code']=run(lambda l:m.lambda_ddot(l,1.0))
res['code f=0 (Newt)']=None
for k,v in res.items():
    if v: print("%-28s t_half=%.5f s  plunge(0.5->0.15)=%.5f s"%(k,v[0],v[1]))
