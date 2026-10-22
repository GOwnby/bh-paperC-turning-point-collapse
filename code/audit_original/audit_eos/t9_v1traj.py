import sys; sys.path.insert(0,'orig')
exec(open('t8_variants.py').read().split("target=27.81035")[0])
from scipy.integrate import solve_ivp
target=27.81035
def build(kw,I):
    kN=kslope(0.0,I,**kw); kF=kslope(1.0,I,**kw); fc=(target-kN)/(kF-kN)
    R0=force(1.0,fc,**kw)
    return fc,lambda lam:(force(lam,fc,**kw)-R0)/I
def run(acc,d0=1e-6):
    ev=lambda t,y:y[0]-0.15; ev.terminal=True
    s=solve_ivp(lambda t,y:[y[1],acc(y[0])],(0,40/np.sqrt(target)),[1-d0,0],method='Radau',events=ev,max_step=0.02,rtol=1e-9,atol=1e-13)
    th=s.t[np.argmax(s.y[0]<0.5)]; return th,s.t[-1]-th
for name,kw,I in [("code",{},I_code),("V1 proper volume, recalibrated",dict(proper_vol=True),I_code),("V1 + I_rel",dict(proper_vol=True),I_rel)]:
    fc,acc=build(kw,I); th,dp=run(acc); print("%-34s f_cal=%.4f t_half=%.5f plunge=%.5f"%(name,fc,th,dp))
