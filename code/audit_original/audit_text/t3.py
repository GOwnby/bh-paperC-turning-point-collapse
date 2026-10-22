import numpy as np, background as bg
from scipy.integrate import solve_ivp
from homologous_calibrated import CalibratedHomologousModel
G,c=bg.G,bg.c
m=CalibratedHomologousModel(5e13)
fc,kN,kF=m.calibrate_to_RR()
print("f_cal",fc,"kN",kN,"kF",kF,"target",m.target_k, "Gamma1avg",m.Gamma1avg)
def run(d0,f_cal,sig2,stops=(0.5,0.173,0.15)):
    m.f_cal=f_cal
    lamcrit=2*G*m.M/(m.R*c**2)
    hor=lambda t,y: (1-2*G*m.M/(y[0]*m.R*c**2))-0.02
    hor.terminal=True; hor.direction=-1
    small=lambda t,y:y[0]-1e-3
    small.terminal=True;small.direction=-1
    smax=(np.log(1/d0)+8)/np.sqrt(sig2)
    s=solve_ivp(lambda t,y:[y[1],m.lambda_ddot(y[0])],(0,smax),[1-d0,0],method='Radau',events=[hor,small],max_step=0.02,rtol=1e-9,atol=1e-13,dense_output=True)
    lam=s.y[0]
    def tcross(L):
        i=np.argmax(lam<L)
        if lam[i]>=L: return np.nan
        # interpolate
        t0,t1=s.t[i-1],s.t[i]; l0,l1=lam[i-1],lam[i]
        return t0+(L-l0)/(l1-l0)*(t1-t0)
    th=tcross(0.5)
    res={L:tcross(L)-th for L in (0.173,0.15,0.1,0.01)}
    return th,s.t[-1]-th,lam[-1],res,s.status,[len(e) for e in s.t_events]
for d0 in (1e-3,1e-6,1e-9,1e-12):
    print(d0, run(d0,fc,m.target_k))
# sensitivity to target -28.51
fc2,_,_=m.calibrate(target_k=28.510545800106442)
print("f_cal(28.51)",fc2,"rel change",(fc2-fc)/fc*100,"%")
a=run(1e-6,fc,m.target_k); b=run(1e-6,fc2,28.5105)
print(a); print(b)
print("plunge change %",(b[1]-a[1])/a[1]*100,"t_half change %",(b[0]-a[0])/a[0]*100)
print("0.5->0.15 change %",(b[3][0.15]-a[3][0.15])/a[3][0.15]*100)
print("lamcrit",2*G*m.M/(m.R*c**2),"Tdyn",np.sqrt(m.R**3/(G*m.M)))
