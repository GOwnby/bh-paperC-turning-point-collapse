import sys; sys.path.insert(0,'orig')
import numpy as np, time
import background as bg
from homologous_calibrated import CalibratedHomologousModel
from gr_rayleigh_ritz import build_RR_matrices
from scipy.linalg import eigh
G,c=bg.G,bg.c
t=time.time()
m=CalibratedHomologousModel(5e13)
fc,kN,kF=m.calibrate_to_RR()
print("target k (RR n=1) = %.5f  (sigma^2=%.5f)"%(m.target_k,-m.target_k))
print("k_N=%.4f k_F=%.4f f_cal=%.5f   R0(f_cal)=%.4f s^-2"%(kN,kF,fc,m.lambda_ddot_raw(1.0,fc)))
print("k(f_cal) =",m.linearized_k(fc))
# linearity in f
for f in (0,0.5,1.0,1.5,2.0): print(" f=%.1f k=%.5f  linear-pred=%.5f"%(f,m.linearized_k(f),kN+f*(kF-kN)))
# FD step dependence
for e in (1e-3,1e-4,1e-5,1e-6): print(" eps",e,"kN",m.linearized_k(0.0,e),"kF",m.linearized_k(1.0,e))
# analytic linearization of U_int  (code W_P = int P 4 pi r^2 dr, no factor 3)
Gam=bg.Gamma1_of_rho(m.rho_bg)
Gavg=np.trapezoid(Gam*m.P_bg*4*np.pi*m.r**2,m.r)/m.W_P
dU=lambda lam: m.minus_dUint_dlam(lam)
e=1e-5; num=(dU(1+e)-dU(1-e))/(2*e)
print("dUint'(1): numeric %.6e analytic 3 W_P (2-3 Gam) = %.6e ratio %.8f ; Gamma1avg=%.6f"%(num,3*m.W_P*(2-3*Gavg),num/(3*m.W_P*(2-3*Gavg)),Gavg))
# anchors
Wvir=3*m.W_P
AN=-np.trapezoid(m.g_eff_newtonian(1.0)*m.r*4*np.pi*m.r**2*m.rho_bg,m.r)
AF=-np.trapezoid(m.g_eff_full(1.0)*m.r*4*np.pi*m.r**2*m.rho_bg,m.r)
print("3W_P=%.6e  A_N=%.6e (ratio %.5f)  A_GR=%.6e (ratio %.8f)"%(Wvir,AN,AN/Wvir,AF,AF/Wvir))
print("compactness GM/Rc^2=%.5f"%(G*m.M/(m.R*c**2)))
# time scale
print("I=%.4e  GM/R^3=%.2f s^-2   Newtonian-only virial: k_N analytic check"%(m.I,G*m.M/m.R**3))
print("sensitivity dk/df_cal =%.3f  => 1%% of f_cal = %.3f s^-2 = %.1f%% of target"%(kF-kN,(kF-kN)*0.01,(kF-kN)*0.01/m.target_k*100))
print("gap to be calibrated: %.3f s^-2 = %.1f%% of target; f_cal-1=%.4f"%(m.target_k-kF,(m.target_k-kF)/m.target_k*100,fc-1))
print("runtime",time.time()-t)
