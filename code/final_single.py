import numpy as np, json
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from homologous_calibrated import CalibratedHomologousModel
import background as bg
G,c=bg.G,bg.c
m=CalibratedHomologousModel(5e13); f,kN,kF=m.calibrate_to_RR()
f5,_,_=m.calibrate(target_k=28.5534); m.f_cal=f
comp_max=np.max(2*G*m.Mr_bg[1:]/(m.r[1:]*c**2))
out=dict(f_cal=f,kN=kN,kF=kF,target=m.target_k,f_cal_n5=f5,lam_crit=comp_max,Tdyn=np.sqrt(m.R**3/(G*m.M)),M=m.M/bg.M_sun,R=m.R,C=2*G*m.M/(m.R*c**2),comp_max=comp_max)
def traj(d0,fc,lam_stop=0.15):
    ev=lambda t,y:y[0]-lam_stop; ev.terminal=True
    s=solve_ivp(lambda t,y:[y[1],m.lambda_ddot(y[0],fc)],(0,12),[1-d0,0],method='Radau',events=ev,max_step=0.02,rtol=1e-9,atol=1e-13,dense_output=True)
    t05=brentq(lambda t:s.sol(t)[0]-0.5,0,s.t[-1]); return t05,s.t[-1]-t05,s
res={}
for d0 in (1e-3,1e-6,1e-9,1e-12):
    t05,dt,s=traj(d0,f); res[d0]=(t05,dt); print(d0,t05,dt,flush=True)
out['traj']={str(k):v for k,v in res.items()}
# stop at 0.173 for two-mode comparison
t05,dt173,_=traj(1e-6,f,0.173); out['dt_0.5_to_0.173']=dt173
# n5 target
t05b,dtb,_=traj(1e-6,f5); out['n5_target']=(t05b,dtb)
# long plunge to the wall
ev=lambda t,y:y[0]-1.02*comp_max; ev.terminal=True
s=solve_ivp(lambda t,y:[y[1],m.lambda_ddot(y[0],f)],(0,6),[1-1e-6,0],method='Radau',events=ev,max_step=0.02,rtol=1e-9,atol=1e-13,dense_output=True)
t05=brentq(lambda t:s.sol(t)[0]-0.5,0,s.t[-1]); out['end']=dict(status=int(s.status),t_end=s.t[-1],lam_end=s.y[0][-1],plunge_to_end=s.t[-1]-t05,rho_amp=1/s.y[0][-1]**3)
np.save('sol_t.npy',s.t); np.save('sol_y.npy',s.y)
print(json.dumps(out,default=float,indent=1))
json.dump(out,open('final_single.json','w'),default=float)
