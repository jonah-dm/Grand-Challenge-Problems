"""Derivation numerics (no SCF): continuation constant d from an LCAO H2+ density (zeta=1.238,
R=2.0 bohr; Finkelstein-Horowitz), STO-6G expansion with correct scaling (1.24^2 removed)."""
import numpy as np
from pyscf import gto, dft
A = 5/9; S0SQ = (6*np.pi)**(-2/3)
E6 = [35.52322122,6.513143725,1.822142904,0.6259552659,0.2430767471,0.1001124280]
C6 = [0.00916359628,0.04936149294,0.16853830490,0.37056279970,0.41649152980,0.13033408410]
def basis(zeta):
    sc = (zeta/1.24)**2
    return gto.basis.parse("H S\n" + "\n".join(f"  {e*sc} {c}" for e,c in zip(E6,C6)))
def parts(atom, zeta, spin):
    mol = gto.M(atom=atom, unit='Bohr', basis={'H':basis(zeta)}, spin=spin, charge=(1 if atom.count('H')==2 and spin==1 else 0), verbose=0)
    S = mol.intor('int1e_ovlp'); c = np.ones(S.shape[0]); c /= np.sqrt(c@S@c); D = np.outer(c,c)
    J = 0.5*np.einsum('ij,kl,ijkl', D, D, mol.intor('int2e'))
    g = dft.gen_grid.Grids(mol); g.level = 8; g.build()
    ao = dft.numint.eval_ao(mol, g.coords, deriv=1); rho = dft.numint.eval_rho(mol, ao, D, xctype='GGA')
    n = rho[0]; sig = (rho[1:4]**2).sum(0); m = n > 1e-14; n, sig, w = n[m], sig[m], g.weights[m]
    s2 = sig/(4*(6*np.pi**2)**(2/3)*n**(8/3)); y = 1 - s2/S0SQ
    ex = -(3/4)*(6/np.pi)**(1/3)*n**(4/3)
    hi = s2 >= S0SQ
    E_hi = np.sum(w*ex*np.where(hi, A/np.where(hi,s2,1), 0))
    E_lo0 = np.sum(w*ex*np.where(hi, 0, A/S0SQ*(1+y)))
    E_lo1 = np.sum(w*ex*np.where(hi, 0, A/S0SQ*y*y))
    Fpbe = 1.804-0.804/(1+0.2195*s2/0.804)
    return J, E_hi, E_lo0, E_lo1, np.sum(w*ex*Fpbe), np.sum(w*n*(~hi))
for label, atom, zeta in [("H atom (zeta=1)", "H 0 0 0", 1.0), ("H2+ LCAO R=2", "H 0 0 -1; H 0 0 1", 1.238)]:
    J, Eh, El0, El1, Ep, wlo = parts(atom, zeta, 1)
    d = (-J - Eh - El0)/El1
    print(f"{label}: -J={-J:.5f}  F1(s>=s0)={Eh:.5f}  lo(d=0)={El0:.5f} lo-per-d={El1:.5f}  => d_norm={d:.3f}  PBE-x={Ep:.5f} weight(s<s0)={wlo:.3f}")
    for dd in [-2, -1, 0, 1]:
        print(f"    d={dd:+d}: E_x^F1 = {Eh+El0+dd*El1:.5f}")

print("--- C0 continuation p(x) = x^m (x = s/s0), F dropping below s0 ---")
from scipy.optimize import brentq
def lo_m(atom, zeta, m):
    mol = gto.M(atom=atom, unit='Bohr', basis={'H':basis(zeta)}, spin=1, charge=(1 if atom.count('H')==2 else 0), verbose=0)
    S = mol.intor('int1e_ovlp'); c = np.ones(S.shape[0]); c /= np.sqrt(c@S@c); D = np.outer(c,c)
    g = dft.gen_grid.Grids(mol); g.level = 8; g.build()
    ao = dft.numint.eval_ao(mol, g.coords, deriv=1); rho = dft.numint.eval_rho(mol, ao, D, xctype='GGA')
    n = rho[0]; sig = (rho[1:4]**2).sum(0); msk = n > 1e-14; n, sig, w = n[msk], sig[msk], g.weights[msk]
    s2 = sig/(4*(6*np.pi**2)**(2/3)*n**(8/3)); ex = -(3/4)*(6/np.pi)**(1/3)*n**(4/3); hi = s2>=S0SQ
    return np.sum(w*ex*np.where(hi, 0, A/S0SQ*(s2/S0SQ)**(m/2)))
J, Eh, *_ = parts("H 0 0 -1; H 0 0 1", 1.238, 1)
f = lambda m: Eh + lo_m("H 0 0 -1; H 0 0 1", 1.238, m) + J
for m in [0, 1, 2, 4, 8]: print(m, f(m))
