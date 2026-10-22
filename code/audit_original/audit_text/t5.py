import numpy as np, background as bg
from metric import build_metric
bg.mu_e=2.0; bg.Zlat=0.0
rp=2.3105195e13
Mp=build_metric(rp)['M']/bg.M_sun
for ratio in (0.5,0.667,0.8,0.9,1.0,1.1,1.25,1.5,2.0):
    M=build_metric(rp*ratio)['M']/bg.M_sun
    print(ratio, M, (Mp-M)/Mp)
