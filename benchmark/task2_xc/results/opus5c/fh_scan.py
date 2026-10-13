import numpy as np
from pyscf import gto, dft
from scipy.optimize import minimize
from h2p_model import basis
from potexact_test import make_F
from ipth import energy, E, dE, zs
S0SQ = (6*np.pi)**(-2/3)
def dens(atom, zeta, spin, charge):
    mol = gto.M(atom=atom, unit='Bohr', basis={'H':basis(zeta)}, spin=spin, charge=charge, verbose=0)
    S = mol.intor('int1e_ovlp'); c = np.ones(S.shape[0]); c /= np.sqrt(c@S@c); D = np.outer(c,c)
    J = 0.5*np.einsum('ij,kl,ijkl', D, D, mol.intor('int2e'))
    g = dft.gen_grid.Grids(mol); g.level = 8; g.build()
    ao = dft.numint.eval_ao(mol, g.coords, deriv=1); rho = dft.numint.eval_rho(mol, ao, D, xctype='GGA')
    n = rho[0]; sig = (rho[1:4]**2).sum(0); m = n > 1e-14
    return J, n[m], sig[m], g.weights[m]
Dh2p = dens("H 0 0 -1; H 0 0 1", 1.238, 1, 1)
Dh2 = dens("H 0 0 -0.7; H 0 0 0.7", 1.166, 0, 0)
def ex(D, F):
    J, n, sig, w = D; s2 = sig/(4*(6*np.pi**2)**(2/3)*n**(8/3))
    return -J, np.sum(w*(-(3/4)*(6/np.pi)**(1/3)*n**(4/3))*F(s2))
for c in [0.10, 0.15, 0.20, 0.25, 0.2795]:
    F0, FH = make_F(c)
    Fx = lambda s2: np.where(s2 < S0SQ, F0, FH(s2))
    c0 = np.zeros(len(zs)); c0[3] = 1.0
    e0 = energy(c0@E, c0@dE, None, lambda s2,p=None: Fx(s2))
    res = minimize(lambda cc: energy(cc@E, cc@dE, None, lambda s2,p=None: Fx(s2)), c0, method='BFGS', options={'maxiter':300})
    print(f"c={c:.4f} F(s0)={F0:.3f}  H2+ (-J, F_H) = {ex(Dh2p, Fx)}  H2-chan = {ex(Dh2, Fx)}  H var lowering {1000*(res.fun-e0):.3f} mEh")
