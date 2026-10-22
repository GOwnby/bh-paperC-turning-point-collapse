import sys; sys.path.insert(0,'var')
import numpy as np
import background as bg
import tov_solve
from metric import build_metric
from homologous_calibrated import CalibratedHomologousModel
from gr_rayleigh_ritz import build_RR_matrices
from scipy.linalg import eigh
bg.mu_e=2.0; bg.Zlat=0.0
print("NEVAL  PCUT   R[km]   M      s2_n1     s2_n5     k_N     k_F    f_cal")
for ne,pc in [(3000,1e-10),(1500,1e-10),(6000,1e-10),(12000,1e-10),(3000,1e-8),(3000,1e-12),(6000,1e-12),(12000,1e-14)]:
    tov_solve.NEVAL=ne; tov_solve.PCUT=pc
    m=CalibratedHomologousModel(5e13)
    A1,B1,_,_=build_RR_matrices(m.bgd,1); s1=eigh(A1,B1)[0][0]
    A5,B5,_,_=build_RR_matrices(m.bgd,5); s5=eigh(A5,B5)[0][0]
    fc,kN,kF=m.calibrate(target_k=-s1)
    print("%5d %.0e %8.3f %.6f %9.4f %9.4f %9.3f %8.4f %.5f"%(ne,pc,m.R/1e3,m.M/bg.M_sun,s1,s5,kN,kF,fc))
