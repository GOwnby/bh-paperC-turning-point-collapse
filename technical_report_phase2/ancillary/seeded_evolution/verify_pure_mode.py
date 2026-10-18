import sympy as sp

t, sigma, C1, C2 = sp.symbols('t sigma C1 C2', real=True, positive=False)
sigma = sp.symbols('sigma', positive=True)

xi_t = C1*sp.exp(sigma*t) + C2*sp.exp(-sigma*t)
v_t = sp.diff(xi_t, t)

# Check: does this solve d^2(xi)/dt^2 = sigma^2 * xi exactly?
residual = sp.simplify(sp.diff(xi_t, t, 2) - sigma**2*xi_t)
print("ODE residual (should be 0):", residual)

xi0 = xi_t.subs(t,0)
v0 = v_t.subs(t,0)
print("\nxi(0) =", xi0)
print("v(0)  =", v0)

# Solve for C1, C2 given xi(0), v(0)
xi0_sym, v0_sym = sp.symbols('xi0 v0')
sol = sp.solve([sp.Eq(xi0, xi0_sym), sp.Eq(v0, v0_sym)], [C1, C2])
print("\nC1, C2 in terms of xi(0), v(0):")
print(sol)

print("\nPure growing mode (C2=0) requires:")
C2_eq0 = sp.Eq(sol[C2], 0)
print(C2_eq0, " => v0 =", sp.solve(C2_eq0, v0_sym))
