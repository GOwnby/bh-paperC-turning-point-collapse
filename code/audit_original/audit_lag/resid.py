import sys, os, numpy as np
SRC='/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve'
os.chdir(SRC); sys.path.insert(0,SRC)
import background as bg; bg.mu_e=2.0
from lag_gr import Lag
G,c=bg.G,bg.c
rc=2.3106026e13
fr=[0.01,0.1,0.3,0.5,0.7,0.9,0.99,0.999]
print('static residual S_j / g_local at mass fractions (g=G m_j/r_j^2), and |S|/g at last 3 nodes')
for N in (100,200,400,800,1600):
    L=Lag(rc,N=N,pert=None); r,U,m=L.unpack(L.y0); g=G*m[1:]/r[1:]**2
    out=[]
    for f in fr:
        j=int(np.argmin(abs(L.mu[1:]/L.mu_tot-f)))
        out.append(L.S[j]/g[j])
    print(N,' '.join('%.1e'%x for x in out),'| last3',' '.join('%.2f'%(L.S[j]/g[j]) for j in (N-3,N-2,N-1)),' mass of outer zone/M %.1e'%(L.dmu[-1]/L.mu_tot))
print('fractions',fr)
# zone densities near surface and TOV-profile
L=Lag(rc,N=200,pert=None); r,U,m=L.unpack(L.y0); Gam,rho,P,e=L.zone(r,U,m)
print('N=200 last 5 zone rho',rho[-5:],' min rho over run initial',rho.min(),' floor 1e-3')
print('Delta r last zones',np.diff(r)[-5:],' R0',L.R0)
print('m_N(discrete)/M_TOV',m[-1]/L.M, ' rest mass mu_tot/M',L.mu_tot/L.M)
