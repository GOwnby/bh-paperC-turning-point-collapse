import sys; sys.path.insert(0,'orig')
import numpy as np
from scipy.optimize import brentq
import background as bg
from metric import build_metric
c=bg.c
# 1) direct lattice-energy construction (independent of B_lat): per-ion Wigner-Seitz energy -0.9 Z^2 e^2/r_i
alpha=1/137.035999084; hbar=1.054571817e-34; e2=alpha*hbar*c
e2_SI=(1.602176634e-19)**2/(4*np.pi*8.8541878128e-12)
print("e^2/(4 pi eps0) from e and eps0: %.6e ; alpha hbar c: %.6e (code e2=%.6e) rel diff %.1e"%(e2_SI,e2,bg.e2,e2_SI/bg.e2-1))
for Z,mue in [(2,2.0),(6,2.0),(8,2.0),(26,56/26)]:
    bg.Zlat=float(Z); bg.mu_e=mue
    ne=3e36   # arbitrary test density
    n_ion=ne/Z; ri=(3/(4*np.pi*n_ion))**(1/3)
    uL_direct=-0.9*Z**2*e2*n_ion/ri   # energy density
    x=(ne/bg.A_n)**(1/3)
    PL_code=-bg.B_lat()*x**4
    # pressure from E(V) derivative: P=-dE/dV, E_cell = -0.9 Z^2 e2/ri(V),  V_cell=Z/ne ; numerically
    Ecell=lambda V:-0.9*Z**2*e2/(3*V/(4*np.pi))**(1/3)
    V0=Z/ne; hh=1e-6*V0
    P_num=-(Ecell(V0+hh)-Ecell(V0-hh))/(2*hh)
    print("Z=%2d: u_L(direct)=%.6e  3*P_L(code)=%.6e  | P_L from -dE/dV=%.6e  P_L(code)=%.6e  ratios %.10f %.10f"%(Z,uL_direct,3*PL_code,P_num,PL_code,uL_direct/(3*PL_code),P_num/PL_code))
    # Pressure zero crossing and density at surface cutoff
    rz=brentq(lambda lr:bg.P_of_rho(np.exp(lr)),np.log(1e3),np.log(1e9))
    print("     P=0 at rho=%.3e kg/m^3 = %.2f g/cc ; ratio P_L/P_e at ultrarel limit = %.4f (=-0.6252*alpha*Z^(2/3)=%.4f)"%(np.exp(rz),np.exp(rz)/1e3,-0.3*(4*np.pi/3)**(1/3)*e2*Z**(2/3)*bg.A_n**(4/3)*1e0/(bg.A_P*2) ,-(6/5)*(4/(9*np.pi))**(1/3)*alpha*Z**(2/3)))
    d=build_metric(5e13)
    P=d['P']; rho=d['rho']
    print("     surface rho(outermost)=%.3e kg/m^3, P_outer/Pc=%.2e; P(rho) monotonic over rho range [rho_surf,rho_c]: %s"%(rho[-1],P[-1]/P[0],bool(np.all(np.diff(bg.P_of_rho(np.logspace(np.log10(rho[-1]),np.log10(5e13),400)))>0))))
