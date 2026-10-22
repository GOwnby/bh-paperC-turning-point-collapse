import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pulse import *
from scipy.optimize import minimize_scalar
def Mof(lr):
    sol,R,yR,xc=background(10**lr,x_end=1e-3)
    return yR[0]*L*c**2/G/Msun, R*L/1e3
if __name__=="__main__":
    for lr in np.arange(12.5,14.51,0.25):
        print(lr, Mof(lr))
    res=minimize_scalar(lambda lr:-Mof(lr)[0],bracket=(13.2,13.36,13.5),tol=1e-10)
    print("Mmax at log10 rho_c=",res.x, 10**res.x, -res.fun, Mof(res.x))
