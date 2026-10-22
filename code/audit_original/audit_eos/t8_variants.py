import sys; sys.path.insert(0,'orig')
import numpy as np
import background as bg
from metric import build_metric
from homologous_calibrated import CalibratedHomologousModel
G,c=bg.G,bg.c
bg.mu_e=2.0; bg.Zlat=0.0; bg.use_eps=True
m=CalibratedHomologousModel(5e13)
r=m.r; Mr=m.Mr_bg; rho_bg=m.rho_bg; I_code=m.I
a=m.bgd['a']; b0=m.bgd['b']; eps0=m.bgd['eps']; P0=m.bgd['P']
eb0=np.exp(b0)
I_rel=np.trapezoid(4*np.pi*np.exp(-a+3*b0)*(eps0+P0/c**2)*r**4,r)
def force(lam, f, proper_vol=False, enthalpy=True):
    """returns generalized force (-dU/dlam - dW/dlam) [kg m^2 s^-2 units], same structure as code."""
    rf=lam*r
    if proper_vol:
        eb=1/np.sqrt(1-2*G*Mr/(rf*c**2))
        rho=rho_bg*eb0/(lam**3*eb)
        w=4*np.pi*r**2*rho_bg*eb0          # dm_rest/dr
    else:
        rho=rho_bg/lam**3
        w=4*np.pi*r**2*rho_bg
    ok=rho>1e-3
    P=np.where(ok,bg.P_of_rho(np.maximum(rho,1e-3)),0.0)
    if proper_vol:
        dlnrho_dlam_neg=3/lam - eb**2*G*Mr/(lam**2*r*c**2)   # = -dln rho/dlam with frozen M_r
        dU=np.trapezoid(np.where(ok,P/np.maximum(rho,1e-30),0)*dlnrho_dlam_neg*w,r)
    else:
        dU=np.trapezoid(np.where(ok,P/np.maximum(rho,1e-30),0)*3/lam*w,r)
    gN=-G*Mr/rf**2
    epsc=bg.eps_over_c2_of_rho(np.maximum(rho,1e-3))
    hfac=(epsc+P/c**2)/np.maximum(rho,1e-30) if enthalpy else 1.0
    gF=-G*hfac*(Mr+4*np.pi*rf**3*P/c**2)/(rf**2*(1-2*G*Mr/(rf*c**2)))
    g=gN+f*(gF-gN)
    dW=np.trapezoid(g*r*w,r)
    return dU+dW
def kslope(f,I,**kw):
    e=1e-5
    return (force(1+e,f,**kw)-force(1-e,f,**kw))/(2*e)/I
target=27.81035
for name,I,kw in [("original (code)",I_code,{}),
                  ("original, I_rel",I_rel,{}),
                  ("gF without enthalpy factor (eps+P)/rho->1",I_code,dict(enthalpy=False)),
                  ("proper-volume rest-mass conservation (V1)",I_code,dict(proper_vol=True)),
                  ("V1 + I_rel",I_rel,dict(proper_vol=True))]:
    kN=kslope(0.0,I,**kw); kF=kslope(1.0,I,**kw)
    fc=(target-kN)/(kF-kN)
    print("%-48s k_N=%9.3f k_F=%8.3f  gap to 27.81=%7.3f  f_cal=%.4f"%(name,kN,kF,target-kF,fc))
print("I_code=%.5e I_rel=%.5e ratio %.4f"%(I_code,I_rel,I_rel/I_code))
# value of I with rest-mass weight
w=4*np.pi*r**2*rho_bg*eb0
print("I_restmass(proper vol) = %.5e ; total rest mass (code dm)=%.6e kg ; proper rest mass=%.6e ; M_grav=%.6e"%(np.trapezoid(w*r**2,r),np.trapezoid(4*np.pi*r**2*rho_bg,r),np.trapezoid(w,r),m.M))
