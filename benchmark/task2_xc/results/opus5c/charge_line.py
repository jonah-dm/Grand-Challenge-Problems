"""Derivation check on the exact H density (no SCF): E_x^SL[q n_H] along the charge line for
F_1 = (5/9)/s^2, PBE-form F (bounded, SCAN-like 1.174 cap shape) — exact is -q^2 * 5/16.
Also Theorem 2 inequality check with libxc PBE correlation on hydrogenic densities."""
import numpy as np
from pyscf.dft import libxc
r = np.linspace(1e-6, 40, 400001); dr = r[1]-r[0]; w = 4*np.pi*r**2*dr
nH = np.exp(-2*r)/np.pi
def s2(n, g): return g**2/(4*(6*np.pi**2)**(2/3)*n**(8/3))
def Ex(n, F):
    g = 2*n; ss = s2(n, g)
    return np.sum(w*(-(3/4)*(6/np.pi)**(1/3)*n**(4/3))*F(ss))
A = 5/9; s0sq = (6*np.pi)**(-2/3)
F1 = lambda ss: np.where(ss >= s0sq, A/ss, A/s0sq*(2-ss/s0sq))
Fpbe = lambda ss: 1.804 - 0.804/(1+0.2195*ss/0.804)
Fscan0 = lambda ss: 1.174*(1-np.exp(-4.9479/ss**0.25))   # SCAN alpha=0 branch shape
print(" q     exact     F1(A/s^2)    PBE-form    SCAN(alpha=0)")
for q in [0.1, 0.25, 0.5, 0.75, 1.0]:
    print(f"{q:4.2f} {-q*q*5/16:10.6f} {Ex(q*nH,F1):10.6f} {Ex(q*nH,Fpbe):10.6f} {Ex(q*nH,Fscan0):10.6f}")
# curvature measure: E(q)-q*E(1) at q=1/2 relative (0 for linear Hxc)
J = 5/16
for name, F in [("F1", F1), ("PBE", Fpbe), ("SCAN0", Fscan0)]:
    hx = lambda q: q*q*J + Ex(q*nH, F)   # Hartree + semilocal exchange, correlation gated to 0
    print(name, "Hx(1/2) - Hx(1)/2 =", hx(0.5) - 0.5*hx(1.0), " (0 = linear; >0 concave... <0 convex delocalization error)")
# Theorem 2: |Ec[nZ/2,nZ/2]| <= |Ec[nZ,nZ]| for PBE correlation (unpolarized input), Z=1
for Z in [1, 2, 5]:
    n1 = Z**3/np.pi*np.exp(-2*Z*r)
    out = []
    for q in [0.5, 1.0]:
        n = 2*q*n1; g = 2*Z*n
        ec = libxc.eval_xc(",PBE", np.array([n, g, 0*n, 0*n]), 0, deriv=0)[0]
        out.append(np.sum(w*n*ec))
    print(f"Z={Z}: PBE Ec[fractional-spin hydrogenic] = {out[0]:.5f}, Ec[2e hydrogenic] = {out[1]:.5f}; needed for constancy with charge-linear x: {-5/32*Z:.4f}")
