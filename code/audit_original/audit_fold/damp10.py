import numpy as np
from scipy.integrate import solve_ivp
exec(open('damp_chk.py').read().split("print('nu")[0])
for S,st,m in ((3000,'first','Radau'),(6050,'first','Radau'),(4000,'rest','Radau')):
    print(10,S,st,pole(10,S,st,m,1e-10),flush=True)
