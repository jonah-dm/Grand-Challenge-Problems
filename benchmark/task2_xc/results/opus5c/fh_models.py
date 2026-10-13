import numpy as np
from pyscf import gto, dft
from h2p_model import basis
from potexact_test import make_F
_, FH = make_F(0.2795)
S0SQ = (6*np.pi)**(-2/3)
def FHx(s2): return np.where(s2 < S0SQ, FH(np.array([S0SQ]))[0], FH(s2))
def run(label, atom, zeta, spin, charge):
    mol = gto.M(atom=atom, unit='Bohr', basis={'H':basis(zeta)}, spin=spin, charge=charge, verbose=0)
    S = mol.intor('int1e_ovlp'); c = np.ones(S.shape[0]); c /= np.sqrt(c@S@c); D = np.outer(c,c)
    J = 0.5*np.einsum('ij,kl,ijkl', D, D, mol.intor('int2e'))
    g = dft.gen_grid.Grids(mol); g.level = 8; g.build()
    ao = dft.numint.eval_ao(mol, g.coords, deriv=1); rho = dft.numint.eval_rho(mol, ao, D, xctype='GGA')
    n = rho[0]; sig = (rho[1:4]**2).sum(0); m = n > 1e-14; n, sig, w = n[m], sig[m], g.weights[m]
    s2 = sig/(4*(6*np.pi**2)**(2/3)*n**(8/3)); ex = -(3/4)*(6/np.pi)**(1/3)*n**(4/3)
    Fp = 1.804-0.804/(1+0.2195*s2/0.804)
    print(f"{label:22s} -J={-J:.5f}  F_H: {np.sum(w*ex*FHx(s2)):.5f}   PBE-form: {np.sum(w*ex*Fp):.5f}")
run("H (zeta=1, STO-6G)", "H 0 0 0", 1.0, 1, 0)
run("H2+ LCAO R=2", "H 0 0 -1; H 0 0 1", 1.238, 1, 1)
for R, z in [(1.4, 1.166), (2.0, 1.10), (3.0, 1.0)]:
    run(f"H2 channel R={R}", f"H 0 0 {-R/2}; H 0 0 {R/2}", z, 0, 0)
