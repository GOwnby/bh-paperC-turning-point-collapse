import numpy as np
from scipy.optimize import brentq, minimize_scalar
import eos as E

# ============================================================
# Neutronization via the semi-empirical mass formula (SEMF), following the same
# physics as Baym-Pethick-Sutherland (1971) / Harrison-Wheeler, in a bounded,
# single-A simplification: fix the mass number A (representative of C/O white
# dwarf material) and let the equilibrium charge Z(n_B) float continuously,
# found by minimizing total energy per baryon (nuclear + electron) at fixed
# baryon density. This captures the correct threshold density (from real
# nuclear binding energy differences) and the correct softening of the EOS,
# without the full multi-nuclide staircase of the complete BPS treatment.
# ============================================================

MeV = 1.602176634e-13   # J
u = 1.66053907e-27      # atomic mass unit, kg
m_p_c2 = 938.272*MeV; m_n_c2 = 939.565*MeV

# SEMF coefficients (MeV), standard values
aV, aS, aC, aA = 15.75, 17.8, 0.711, 23.7

def nuclear_mass_c2(Z, A):
    # binding energy (MeV), semi-empirical mass formula (no pairing term -- Z continuous)
    Zc = Z
    B = aV*A - aS*A**(2/3) - aC*Zc**2/A**(1/3) - aA*(A-2*Zc)**2/A
    Mc2 = Zc*m_p_c2 + (A-Zc)*m_n_c2 - B*MeV
    return Mc2   # Joules

def total_energy_per_baryon(Z, A, n_B):
    # n_B: baryon number density (1/m^3). Charge neutrality: n_e = Z/A * n_B.
    n_e = (Z/A)*n_B
    if n_e <= 0:
        e_elec_per_baryon = 0.0
    else:
        x = (3*n_e/(8*np.pi))**(1/3) * E.h/(E.m_e*E.c)
        e_elec = E.g_cs(x) * (np.pi*E.m_e**4*E.c**5)/(3*E.h**3)   # J/m^3, electron energy density (kinetic+rest)
        e_elec_per_baryon = e_elec/n_B
    e_nuclear_per_baryon = nuclear_mass_c2(Z, A)/A
    return e_nuclear_per_baryon + e_elec_per_baryon

def equilibrium_Z(A, n_B):
    res = minimize_scalar(lambda Z: total_energy_per_baryon(Z, A, n_B),
                           bounds=(0.05*A, 0.55*A), method='bounded',
                           options={'xatol':1e-6})
    return res.x

A_rep = 12.0   # representative mass number (Carbon-Oxygen WD material)

# sanity check: at very low density (n_e -> 0), equilibrium Z should match the
# ordinary (zero-pressure) beta-stability value from the SEMF alone (Z_0 = A/(2+aC A^(2/3)/(2 aA)))
Z0 = A_rep/(2 + aC*A_rep**(2/3)/(2*aA))
print(f"SEMF zero-density beta-stability Z for A={A_rep}: {Z0:.3f}  (Carbon-12 has Z=6, close as expected)")

print("\nScanning baryon density for the equilibrium Z(n_B) / Y_e(n_B):")
n_B_scan = np.geomspace(1e28, 1e33, 20)   # 1/m^3;  rho ~ n_B * m_u roughly
for nB in n_B_scan:
    Zeq = equilibrium_Z(A_rep, nB)
    Ye = Zeq/A_rep
    rho_approx = nB * u   # kg/m^3 (baryon rest mass density, ignoring small binding energy correction)
    print(f"  n_B={nB:.3e} 1/m^3  (rho~{rho_approx:.3e} kg/m^3)   Z_eq={Zeq:.4f}   Y_e={Ye:.4f}")

print("\n\n=== Extended scan to find the actual neutronization threshold ===")
n_B_scan2 = np.geomspace(1e33, 1e41, 30)
for nB in n_B_scan2:
    Zeq = equilibrium_Z(A_rep, nB)
    Ye = Zeq/A_rep
    rho_approx = nB * u
    print(f"  n_B={nB:.3e} 1/m^3  (rho~{rho_approx:.3e} kg/m^3)   Z_eq={Zeq:.4f}   Y_e={Ye:.4f}")
