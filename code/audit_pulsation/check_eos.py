from eos import *
import numpy as np
for x in [1e-3,0.1,1,3,10,30]:
    # n d eps/dn = eps+P ; n ∝ x^3 -> n d/dn = (x/3) d/dx
    hh=x*1e-5
    de=(eps(x+hh)-eps(x-hh))/(2*hh)
    lhs=x/3*de; rhs=eps(x)+P(x)
    # Gamma numeric
    g=(np.log(P(x+hh))-np.log(P(x-hh)))/(np.log(rho(x+hh))-np.log(rho(x-hh)))
    # dP/dx analytic check
    dp=(P(x+hh)-P(x-hh))/(2*hh)
    print(x, lhs/rhs-1, g/Gamma1(x)-1, dp/dPdx(x)-1)
# small-x P check: P -> (8/15)K x^5
x=1e-2; print(P(x)/(8/15*K*x**5))
print("x for 5e13:",x_of_rho(5e13), "rho_c range x:",[x_of_rho(r) for r in (1e13,2e13,2.31e13,3e13)])
