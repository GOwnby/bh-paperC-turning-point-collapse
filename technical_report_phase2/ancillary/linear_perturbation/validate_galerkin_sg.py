import numpy as np
from scipy.linalg import eigh
from galerkin_fem import assemble
from validate_polytrope import build_polytrope, G

def count_nodes(xi):
    s = np.sign(xi[2:-2]); s = s[s!=0]
    return np.sum(np.abs(np.diff(s))>0)

print("Galerkin WITH self-gravity (exact l=0 closed form) vs Pekeris Gamma1=4/3:\n")
K = 1e10
for Gamma in [1.20, 1.25, 1.30, 4/3, 1.40, 1.50, 5/3]:
    sol = build_polytrope(K, Gamma, rho_c=1e3, rmax=5e8)
    r_raw, m_raw, rho_raw = sol.t, sol.y[0], np.maximum(sol.y[1], 1e-10)
    N = 300
    r = np.concatenate([[0.0], np.linspace(r_raw[1], r_raw[-1]*0.9995, N-1)])
    rho0 = np.interp(r, r_raw, rho_raw); rho0[0] = rho_raw[0]
    P0 = K*rho0**Gamma
    Gamma1 = np.full(len(r), Gamma)

    M, D = assemble(r, rho0, P0, Gamma1, G=G, include_self_gravity=True)
    M = M[1:,1:]; D = D[1:,1:]
    eigvals, eigvecs = eigh(M, D)
    order = np.argsort(eigvals)
    fund = None
    nodeinfo = []
    for idx in order[:25]:
        xi = eigvecs[:,idx]
        nn = count_nodes(xi)
        nodeinfo.append((eigvals[idx], nn))
        if fund is None and nn == 0:
            fund = eigvals[idx]
    verdict = "UNSTABLE" if (fund is not None and fund < 0) else "stable" if fund is not None else "NOT FOUND"
    fund_str = f"{fund:+.5e}" if fund is not None else "None"
    print(f"Gamma1={Gamma:.4f} (3G-4={3*Gamma-4:+.3f}): fundamental={fund_str}  {verdict}   nodes={[(round(float(e),12),int(n)) for e,n in nodeinfo[:5]]}")
