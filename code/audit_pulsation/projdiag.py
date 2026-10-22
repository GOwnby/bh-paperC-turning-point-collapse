import sys
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import numpy as np
import background as bg
from metric import build_metric
import gr_rayleigh_ritz as rr
b=build_metric(5e13)
print("proj R km",b['R']/1e3,"M/Msun(1.98892e30)",b['M']/bg.M_sun, "M/(1.98841e30)",b['M']/1.98841e30, "M kg",b['M'])
r=b['r']; print("n pts",len(r),"r0",r[0],"dr",r[1]-r[0])
print("a(0),a(R)",b['a'][0],b['a'][-1])
A,B,xis,rr_=rr.build_RR_matrices(b,5)
print("A diag",np.diag(A)); print("B diag",np.diag(B))
print("A asym rel", np.max(np.abs(A-A.T))/np.max(np.abs(A)))
print(A[0,:], A[:,0])
