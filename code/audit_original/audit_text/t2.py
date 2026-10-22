import numpy as np
from figs_nl import *
for rc in [1e8,3e8,1e9,3e9,1e10,3e10,1e11,3e11,1e12,3e12,1e13]:
    sN,C=newtonian_s2(rc)
    d=build_metric(rc); A,B,_,_=build_RR_matrices(d,5); sG=eigh(A,B)[0][0]
    print(f"rc={rc:.0e} compact={2*G*d['M']/(d['R']*c**2):.3e} sN={sN:.5g} sG={sG:.5g} ratio={sG/sN:.5f}",flush=True)
