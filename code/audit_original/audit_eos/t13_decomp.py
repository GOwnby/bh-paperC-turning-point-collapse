import sys; sys.path.insert(0,'orig')
import numpy as np
import background as bg
from homologous_calibrated import CalibratedHomologousModel
from exactV import Star
G,c=bg.G,bg.c
bg.mu_e=2.0; bg.Zlat=0.0; bg.use_eps=True
st=Star(5e13)
m=CalibratedHomologousModel(5e13); fc,kN,kF=m.calibrate_to_RR()
r=m.r; I=m.I
# use code's grid, evaluate background from code arrays (eb0 from M_r)
eb0=1/np.sqrt(1-2*G*m.Mr_bg/(r*c**2)); wrest=4*np.pi*r**2*m.rho_bg*eb0
def exact_fields(lam):
    M,s=st.M_of_lam(lam,full=True)
    Mr=s.sol(r)[0]
    eb=1/np.sqrt(1-2*G*Mr/(lam*r*c**2))
    rho=m.rho_bg*eb0/(lam**3*eb)
    return M,Mr,rho
def U_exact(lam):
    M,Mr,rho=exact_fields(lam)
    u=(bg.eps_over_c2_of_rho(rho)-rho)*c**2/rho
    return np.trapezoid(u*wrest,r)
def Uc(lam): # code's internal energy (same u, but rho=rho_bg/lam^3, dm=4pi r^2 rho_bg dr)
    rho=m.rho_bg/lam**3; u=(bg.eps_over_c2_of_rho(rho)-rho)*c**2/rho
    return np.trapezoid(u*4*np.pi*r**2*m.rho_bg,r)
print(" lam | F_int exact/I | F_int code/I | F_grav exact/I | F_grav code(f_cal)/I | F_grav code(f=1)/I")
rows=[]
for lam in (1.0,0.99,1.01,0.9,0.7,0.5,0.3):
    h=1e-3
    Fint_ex=-(U_exact(lam+h)-U_exact(lam-h))/(2*h)/I
    Fint_c=-(Uc(lam+h)-Uc(lam-h))/(2*h)/I
    Ftot=-(st.M_of_lam(lam+h)-st.M_of_lam(lam-h))/(2*h)*c**2/I
    Fgr_ex=Ftot-Fint_ex
    Fint_c2=m.minus_dUint_dlam(lam)/I
    Fgr_c=m.minus_dWgrav_dlam(lam,fc)/I; Fgr_1=m.minus_dWgrav_dlam(lam,1.0)/I
    rows.append((lam,Fint_ex,Fint_c2,Fgr_ex,Fgr_c,Fgr_1))
    print(" %.2f  %13.6e %13.6e %14.6e %14.6e %14.6e   (FD int check %.6e)"%(lam,Fint_ex,Fint_c2,Fgr_ex,Fgr_c,Fgr_1,Fint_c))
# linear slopes at lam=1 from the +/-0.01 rows
d={l:(a,b,cc,dd,e) for l,a,b,cc,dd,e in rows}
sl=lambda i:(d[1.01][i]-d[0.99][i])/0.02
print("slopes d/dlam at 1 (s^-2): int exact %.3f | int code %.3f | grav exact %.3f | grav code f_cal %.3f | grav code f=1 %.3f"%(sl(0),sl(1),sl(2),sl(3),sl(4)))
print("total exact slope %.3f ; code(f_cal) %.3f ; code(f=1) %.3f"%(sl(0)+sl(2),sl(1)+sl(3),sl(1)+sl(4)))
