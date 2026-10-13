"""Derivation check (no SCF): homogeneous iso-orbital exchange e = -K n^4/|grad n|^2 (F = A s^-2)
on one-orbital densities: H atom (norm), Wang H2 sigma_g (zeta=1.166, R=1.4), compared with exact -J[n_sigma].
J computed with PySCF integrals on an STO-6G expansion of the Slater functions (no SCF)."""
import numpy as np
from pyscf import gto, dft
# K for spin-polarized one-electron channel: E_x = -K int n^4/|grad n|^2 with F = A/s^2
cK = 0.375*(3/np.pi)**(1/3)*2**(4/3)*4*(6*np.pi**2)**(2/3)
A = (5/16)/(cK/(32*np.pi)); s0 = (6*np.pi)**(-1/3)
print("A =",A," s0 =",s0," F(s0)=",A/s0**2)
def sfun(n, g): return g/(2*(6*np.pi**2*n)**(1/3)*n)
def run(atom, zeta, label):
    mol = gto.M(atom=atom, basis={'H': gto.basis.parse('''H S
  %s 0.00916359628
  %s 0.04936149294
  %s 0.16853830490
  %s 0.37056279970
  %s 0.41649152980
  %s 0.13033408410''' % tuple(e*zeta**2 for e in [35.52322122,6.513143725,1.822142904,0.6259552659,0.2430767471,0.1001124280]))}, spin=None if mol_spin(atom)==0 else 1, verbose=0)
    S = mol.intor('int1e_ovlp'); nb = S.shape[0]
    c = np.ones(nb); c /= np.sqrt(c@S@c)
    D = np.outer(c,c)
    eri = mol.intor('int2e')
    J = 0.5*np.einsum('ij,kl,ijkl', D, D, eri)
    g = dft.gen_grid.Grids(mol); g.level = 7; g.build()
    ao = dft.numint.eval_ao(mol, g.coords, deriv=1)
    rho = dft.numint.eval_rho(mol, ao, D, xctype='GGA')
    n = rho[0]; gn = np.linalg.norm(rho[1:4], axis=0)
    m = n > 1e-12
    s = sfun(n[m], gn[m]); w = g.weights[m]
    exlda_pol = -(3/4)*(6/np.pi)**(1/3)*n[m]**(4/3)
    out = {}
    for name, F in [("A/s^2 (cap at 2-x^2 below s0)", np.where(s>=s0, A/s**2, A/s0**2*(2-(s/s0)**2))),
                    ("A/s^2 switch-to-0 below s0", np.where(s>=s0, A/s**2, 0.0)),
                    ("PBE F (mu=0.2195,k=0.804)", 1.804-0.804/(1+0.2195*s**2/0.804))]:
        out[name] = np.sum(w*exlda_pol*F)
    print(label, "N=%.4f"%np.sum(w*n[m]), "exact -J = %.5f"%(-J), {k:"%.5f"%v for k,v in out.items()},
          "frac weight s<s0: %.3f"%(np.sum(w*n[m]*(s<s0))))
def mol_spin(a): return 1 if a.count('H')==1 else 0
run('H 0 0 0', 1.0, 'H atom ')
run('H 0 0 0; H 0 0 0.7408481', 1.166, 'H2 Wang')
