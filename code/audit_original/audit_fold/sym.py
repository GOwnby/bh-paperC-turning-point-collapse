import sympy as sp
# (1) normal form scaling
y,t,k,a,Md,al,be,u,s=sp.symbols('y t kappa a Mdot alpha beta u s',positive=True)
# y=al*u(s), t=be*s : y'' = al/be^2 u''
# eq: al/be^2 u'' = k/2 al^2 u^2 + k/(2a) Md be s
c1=sp.simplify((k/2*al**2)/(al/be**2)); c2=sp.simplify((k/(2*a)*Md*be)/(al/be**2))
sol=sp.solve([sp.Eq(c1,1),sp.Eq(c2,1)],[al,be],dict=True); print(sol)
b5=sp.simplify(be**5); print('beta^5 check', sp.simplify(c2.subs(al,2/(k*be**2))-k**2*Md*be**5/(4*a)))
# (2) slow manifold series for u''+nu u' = u^2+s,  s=-X^2 ; u = -X + sum d_j X^{-m}
X,nu=sp.symbols('X nu',positive=True)
# s=-X^2, d/ds = -1/(2X) d/dX
D=lambda f: -sp.diff(f,X)/(2*X)
ds=sp.symbols('d1:6')
# ansatz in powers: for nu=0: terms X^-4, X^-9 ... ; with nu: X^-2 (nu), X^-4,...  general: powers X^-2,-3,-4,-5.. use all
cs=sp.symbols('c2:9')
U=-X+sum(c*X**(-j) for c,j in zip(cs,range(2,9)))
eq=sp.expand(D(D(U))+nu*D(U)-U**2+X**2)
eq=sp.series(eq,X,sp.oo,10).removeO() if False else eq
# collect in powers of 1/X
w=sp.symbols('w',positive=True)
e2=sp.expand(eq.subs(X,1/w))
sols={}
for p in range(-1,8):
    co=sp.expand(e2).coeff(w,p)
    co=co.subs(sols)
    if co==0: continue
    fr=[c for c in cs if co.has(c)]
    if not fr: print('residual at w^',p,co); continue
    c_=sorted(fr,key=lambda q:int(str(q)[1:]))[0] if False else fr[0]
    r=sp.solve(co,c_)
    if r: sols[c_]=sp.simplify(r[0])
print(sols)
