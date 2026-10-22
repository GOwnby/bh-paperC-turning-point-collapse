import numpy as np, background as bg
exec(open('t3.py').read().split("for d0 in (1e-3")[0])
fc2,_,_=m.calibrate(target_k=28.510545800106442)
print("f_cal(28.51)",fc2,"rel change",(fc2-fc)/fc*100,"%")
a=run(1e-6,fc,27.81034874841159); b=run(1e-6,fc2,28.5105458)
print(a); print(b)
print("plunge(to end) change %",(b[1]-a[1])/a[1]*100,"t_half change %",(b[0]-a[0])/a[0]*100)
print("0.5->0.15 change %",(b[3][0.15]-a[3][0.15])/a[3][0.15]*100, "0.5->0.173",(b[3][0.173]-a[3][0.173])/a[3][0.173]*100)
