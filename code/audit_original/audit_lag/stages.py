import sys, os, json, numpy as np, io, contextlib
SRC='/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/collapse_solve'
os.chdir(SRC); sys.path.insert(0,SRC)
import background as bg; bg.mu_e=2.0
with contextlib.redirect_stdout(io.StringIO()):
    from fold_inhom import handoff, sc
    from fold_dyn import F
from scipy.integrate import solve_ivp
th,state=handoff(sc,yh=0.1)
levels=[0.9,0.8,0.7,0.6,0.5,0.3,0.15]
def hom(scale):
    evs=[(lambda t,y,v=v:y[0]-v) for v in levels]
    for e in evs: e.terminal=False
    s=solve_ivp(lambda t,y:[y[1],scale*F(y[0])],(0,10),list(state),method='DOP853',events=evs,rtol=1e-11,atol=1e-14,max_step=0.01)
    return [s.t_events[i][0] for i in range(len(levels))]
A='/tmp/claude-0/-home-claude/717c308b-4c66-5933-99e1-d42d1016a975/scratchpad/audit_lag/'
def lag(N):
    d=json.load(open(A+'ev_N%d_cfl0.4_y0.1.json'%N)); s=d['series']
    t=np.array([x['t'] for x in s]); l=np.array([x['lam_mass'] for x in s])
    o=[]
    for v in levels:
        i=np.argmax(l<v); o.append(np.interp(v,[l[i],l[i-1]],[t[i],t[i-1]]))  # coarse (sample spacing 0.1/0.002 s)
    return o
print('levels   ',levels)
print('hom resc ',np.round(hom(sc),4))
print('hom raw  ',np.round(hom(1.0),4))
for N in (400,800): print('lag N=%d '%N,np.round(lag(N),4))
