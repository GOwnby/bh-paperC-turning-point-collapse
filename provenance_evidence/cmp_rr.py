import sys, numpy as np, importlib
def load(path):
    sys.path.insert(0,path)
    for m in ['background','metric','tov_solve','gr_rayleigh_ritz']:
        sys.modules.pop(m,None)
    import background as bg, metric, gr_rayleigh_ritz as rr
    bg.mu_e=2.0; bg.Zlat=0.0
    d=metric.build_metric(5e13); out=[rr.build_RR_matrices(d,n)[:2] for n in (1,3,5)]
    sys.path.remove(path)
    return out
a=load('old_code'); b=load('newcopy/code')
for (A1,B1),(A2,B2) in zip(a,b): print('max |dA|,|dB| =',np.max(np.abs(A1-A2)),np.max(np.abs(B1-B2)))
