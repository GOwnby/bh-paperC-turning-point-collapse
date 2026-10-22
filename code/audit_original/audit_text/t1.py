import time, numpy as np, background as bg
from metric import build_metric
from gr_rayleigh_ritz import build_RR_matrices
from scipy.linalg import eigh
t=time.time()
bg.mu_e=2.0; bg.Zlat=0.0
d=build_metric(5e13)
print("time",time.time()-t)
print("M",d['M']/bg.M_sun,"R",d['R'],"C",2*bg.G*d['M']/(d['R']*bg.c**2))
G1=np.array([bg.Gamma1_of_rho(x) if x>1 else 5/3 for x in d['rho']])
print("Gamma1 min",G1.min(),"center",G1[0],"surf",G1[-5], "argmin",G1.argmin())
print("Gamma1 max", G1.max())
vals=[]
for n in range(1,6):
    A,B,_,_=build_RR_matrices(d,n); vals.append(eigh(A,B)[0][0])
print(vals)
# two trial functions: xi=r and xi=r e^a
r=d['r']; a=d['a']; 
import gr_rayleigh_ritz as rr
A,B,xis,r=rr.build_RR_matrices(d,1)
print("xi=r",A[0,0]/B[0,0])
# xi=r e^a
b=d['b']; P=d['P']; eps=d['eps']; dPdr=d['dPdr']; rho=d['rho']
Pi=G1*P*np.exp(3*a+b)/r**2; Q=4*dPdr*r*np.exp(a+b)-dPdr**2*r**2*np.exp(a+b)/(P+eps*bg.c**2)+8*np.pi*bg.G/bg.c**4*(P+eps*bg.c**2)*P*r**2*np.exp(a+3*b)
S=np.exp(-a+3*b)*(1+P/(eps*bg.c**2))*eps*r**2
xi=r*np.exp(a); z=np.exp(-a)*r**2*xi; dz=np.gradient(z,r)
print("xi=r e^a", np.trapezoid(Pi*dz**2+Q*xi**2,r)/np.trapezoid(S*xi**2,r))
np.save('d_r.npy',r)
# eigenfunction
A,B,xis,r=rr.build_RR_matrices(d,5)
w,v=eigh(A,B); c0=v[:,0]; xi=sum(c0[i]*xis[i] for i in range(5)); ratio=xi/r
print("xi/r center, surface-ish", ratio[1], ratio[-2], ratio[1]/ratio[-2], ratio[1]/ratio[-1])
print("tail",ratio[-5:])
