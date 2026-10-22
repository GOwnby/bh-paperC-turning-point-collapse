import sys; sys.path.insert(0,'/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/audit_eos/orig')
import numpy as np
from fractions import Fraction
import background as bg
c=bg.c
def check(Z,mue,label):
    bg.Zlat=float(Z); bg.mu_e=mue; bg.use_eps=True
    print("==",label,"Z=",Z,"mu_e=",mue)
    for rho in [1e3,1e6,1e8,1e9,1e10,1e11,1e12,1e13,2.3e13,5e13,1e14,1e15]:
        h=1e-4
        eps=lambda r: bg.eps_over_c2_of_rho(r)
        # rho d(eps/c2)/drho  vs  eps/c2 + P/c2  (n d eps/dn = eps + P with n propto rho)
        d=(eps(rho*(1+h))-eps(rho*(1-h)))/(2*h)   # = rho d eps/drho  (central diff in ln rho)
        d4=(-eps(rho*(1+2*h))+8*eps(rho*(1+h))-8*eps(rho*(1-h))+eps(rho*(1-2*h)))/(12*h)
        rhs=eps(rho)+bg.P_of_rho(rho)/c**2
        # Gamma1 analytic vs numeric
        P=bg.P_of_rho
        Gn=(np.log(P(rho*(1+h)))-np.log(P(rho*(1-h))))/(np.log(1+h)-np.log(1-h)) if P(rho)>0 else np.nan
        Ga=bg.Gamma1_of_rho(rho)
        print(f" rho={rho:8.1e} P/rhoc2={P(rho)/rho/c**2: .3e} first-law rel.err={(d4-rhs)/rhs: .2e}  Gamma1 ana={Ga:.8f} num={Gn:.8f} diff={Ga-Gn: .1e}")
for Z,mue,l in [(0,2.0,'ideal'),(6,2.0,'C12'),(26,56/26,'Fe56')]:
    check(Z,mue,l)
# Also: use_eps=False inconsistency
bg.Zlat=0.;bg.mu_e=2.;bg.use_eps=False
print("use_eps False: eps=rho c^2 -> first law residual = -P/c^2 (relative P/(rho c^2)):")
for rho in [2.3e13,5e13]: print(rho, bg.P_of_rho(rho)/rho/c**2)
