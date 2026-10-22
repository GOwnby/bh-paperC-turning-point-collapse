import sys, numpy as np
sys.path.insert(0,'/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve')
import background as bg; bg.mu_e=2.0
from metric import build_metric
from gr_rayleigh_ritz import build_RR_matrices
from scipy.linalg import eigh
for rc in (5e13,1e13,2.3106026e13):
    d=build_metric(rc)
    for nb in (3,5,7):
        A,B,_,_=build_RR_matrices(d,nb); print(rc,nb,eigh(A,B)[0][0], flush=True)
    print('M/Msun',d['M']/bg.M_sun,'R km',d['R']/1e3,'2GM/Rc2',2*bg.G*d['M']/(d['R']*bg.c**2))
