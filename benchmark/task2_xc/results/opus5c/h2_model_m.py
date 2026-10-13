import numpy as np
from pyscf import gto, dft
from h2p_model import basis, A, S0SQ
def run(atom, zeta, ms):
    mol = gto.M(atom=atom, unit='Bohr', basis={'H':basis(zeta)}, spin=0, verbose=0)
    S = mol.intor('int1e_ovlp'); c = np.ones(S.shape[0]); c /= np.sqrt(c@S@c); D = np.outer(c,c)  # one orbital, one electron per channel
    J = 0.5*np.einsum('ij,kl,ijkl', D, D, mol.intor('int2e'))
    g = dft.gen_grid.Grids(mol); g.level = 8; g.build()
    ao = dft.numint.eval_ao(mol, g.coords, deriv=1); rho = dft.numint.eval_rho(mol, ao, D, xctype='GGA')
    n = rho[0]; sig = (rho[1:4]**2).sum(0); msk = n > 1e-14; n, sig, w = n[msk], sig[msk], g.weights[msk]
    s2 = sig/(4*(6*np.pi**2)**(2/3)*n**(8/3)); ex = -(3/4)*(6/np.pi)**(1/3)*n**(4/3); hi = s2>=S0SQ
    out = {m: np.sum(w*ex*np.where(hi, A/np.where(hi,s2,1), A/S0SQ*(s2/S0SQ)**(m/2))) for m in ms}
    Fp = 1.804-0.804/(1+0.2195*s2/0.804)
    print(atom, "zeta", zeta, "-J =", round(-J,5), {m: round(v,5) for m,v in out.items()}, "PBE", round(np.sum(w*ex*Fp),5))
for R in [1.4, 2.0, 3.0]:
    run(f"H 0 0 {-R/2}; H 0 0 {R/2}", 1.166 if R==1.4 else (1.10 if R==2.0 else 1.0), [1.75, 2.0])
