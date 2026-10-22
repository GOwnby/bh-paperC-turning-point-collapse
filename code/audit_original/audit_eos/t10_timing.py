import sys; sys.path.insert(0,'orig')
import numpy as np
import background as bg
from homologous_calibrated import CalibratedHomologousModel
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
bg.mu_e=2.0; bg.Zlat=0.0
m=CalibratedHomologousModel(5e13); fc,kN,kF=m.calibrate_to_RR(); print("f_cal",fc)
def go(ms,rtol,d0=1e-6,method='Radau'):
    ev=lambda t,y:y[0]-0.15; ev.terminal=True
    s=solve_ivp(lambda t,y:[y[1],m.lambda_ddot(y[0])],(0,40/np.sqrt(m.target_k)),[1-d0,0],method=method,events=ev,max_step=ms,rtol=rtol,atol=1e-13,dense_output=True)
    th_sample=s.t[np.argmax(s.y[0]<0.5)]
    i=np.argmax(s.y[0]<0.5); th_exact=brentq(lambda t:s.sol(t)[0]-0.5,s.t[i-1],s.t[i])
    # time to 0.15 (event) and 0.173
    j=np.argmax(s.y[0]<0.173); t173=brentq(lambda t:s.sol(t)[0]-0.173,s.t[j-1],s.t[j])
    return th_sample,th_exact,s.t[-1],s.t[-1]-th_sample,s.t[-1]-th_exact,s.t[-1]-t173,s.t[-1]-th_exact
for ms,rt,meth in [(0.02,1e-9,'Radau'),(0.005,1e-10,'Radau'),(0.02,1e-12,'DOP853'),(0.002,1e-12,'DOP853')]:
    o=go(ms,rt,method=meth)
    print("max_step=%g rtol=%g %s: t_half(sample)=%.5f t_half(exact)=%.5f t(0.15)=%.5f | plunge(sample)=%.5f plunge(exact 0.5->0.15)=%.5f ; 0.173->0.15: %.5f"%(ms,rt,meth,o[0],o[1],o[2],o[3],o[4],o[5]))
