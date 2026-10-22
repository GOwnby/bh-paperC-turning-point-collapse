import sys; sys.path.insert(0,'orig')
import numpy as np
import background as bg
from twomode import TwoMode
from scipy.integrate import quad
G,c=bg.G,bg.c
bg.mu_e=2.0; bg.Zlat=0.0
t=TwoMode()
m=t.m
print("f_cal used in TwoMode:",t.f_cal," (paper text says 1.272)")
print("Mmat",t.Mmat.tolist())
# 1) independent internal-energy potential: u(rho) numerically = int_0^rho P/rho'^2 drho' ; check against first-law (eps/c^2-rho)c^2/rho
def u_num(rho):
    f=lambda lr: bg.P_of_rho(np.exp(lr))/np.exp(lr)  # d(rho)/rho^2 = dlnrho/rho
    return quad(f,np.log(1e-2),np.log(rho),epsabs=0,epsrel=1e-12)[0]   # lower limit 1e-2 kg/m3 (negligible)
for rho in (1e9,1e12,2e13,5e13):
    uf=(bg.eps_over_c2_of_rho(rho)-rho)*c**2/rho
    print(" rho=%.0e  u_quad=%.8e  u_firstlaw=(eps/c2-rho)c2/rho=%.8e rel diff %.2e"%(rho,u_num(rho),uf,u_num(rho)/uf-1))
# 2) generalized force from FD of U_int(q), using u from first law (exact), dm=w dr
r,w=t.r,t.w
def Uint(q):
    s,J=t.fields(q); rho=m.rho_bg/(s**2*J)
    u=(bg.eps_over_c2_of_rho(rho)-rho)*c**2/rho
    return np.trapezoid(u*w,r)
def Qint_code(q):
    # replicate Q_raw int-part (without gravity, without R0)
    s,J=t.fields(q); rho=m.rho_bg/(s**2*J); P=np.where(rho>1e-3,bg.P_of_rho(np.maximum(rho,1e-3)),0.)
    out=np.zeros(2)
    for k in range(2):
        phik=t.phi[k]; rd=0 if k==0 else -2*t.x2
        out[k]=np.trapezoid(P/rho*(2*phik/s+(phik+rd)/J)*w,r)
    return out
for q in ([1.0,0.0],[0.9,0.03],[0.7,-0.1],[0.5,0.05]):
    q=np.array(q); e=1e-6
    fd=np.array([-(Uint(q+e*np.eye(2)[i])-Uint(q-e*np.eye(2)[i]))/(2*e) for i in range(2)])
    cd=Qint_code(q)
    print(" q=",q," -dU/dq (FD, first-law u) =",fd," code Qint =",cd," rel diff",(cd-fd)/np.abs(fd))
# 3) mass matrix from T = 1/2 int |rdot|^2 dm with rdot = r*(lamdot*1 + mudot*(1-x^2))
M=np.array([[np.trapezoid(r**2*t.phi[i]*t.phi[j]*w,r) for j in range(2)] for i in range(2)])
print("Mmat recomputed:",np.allclose(M,t.Mmat), "  M00 == I_single?",M[0,0],m.I)
# 4) single-mode slice: Q_0(lam,0)/M00 vs lambda_ddot_raw (with f_cal) ; both raw (before R0)
t.R0_save=t.R0.copy(); t.R0=np.zeros(2)
for lam in (1.0,0.9,0.5,0.2):
    Q=t.Q_raw(np.array([lam,0.0]))
    print(" lam=%.2f  Q0/M00=%.8e  single-mode raw=%.8e   (Q1(lam,0)=%.3e  -> mu is forced even if seeded at 0)"%(lam,Q[0]/t.Mmat[0,0],m.lambda_ddot_raw(lam,t.f_cal),Q[1]))
t.R0=t.R0_save
print("R0 vector (two-mode):",t.R0, " single R0*I:",m.lambda_ddot_raw(1.0,t.f_cal)*m.I)
K=t.jac(1e-5); w_,v=np.linalg.eig(K); print("K=",K.tolist()); print("eig",w_); print("vec",v)
print("single-mode k",m.linearized_k(t.f_cal)," K00",K[0,0])
