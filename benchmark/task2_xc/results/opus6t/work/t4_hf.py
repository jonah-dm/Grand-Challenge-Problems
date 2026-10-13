"""T4' on HF densities: forced-form bound E_c^sl <= -K int_{s>=s*} n^{4/3} s^-2 for two-electron singlets."""
import numpy as np
from pyscf import gto, scf, dft
K = 5*np.pi/(4*(3*np.pi**2)**(2/3)); sstar = (3*np.pi)**(-1/3)
for name, atom, ref in [('He', 'He 0 0 0', -0.0420), ('H2', 'H 0 0 0; H 0 0 0.741', -0.0409), ('Li+', 'Li 0 0 0', -0.0435)]:
    mol = gto.M(atom=atom, basis='aug-cc-pvqz' if name != 'Li+' else 'cc-pvqz', charge=1 if name == 'Li+' else 0, verbose=0)
    mf = scf.RHF(mol).run()
    g = dft.gen_grid.Grids(mol); g.level = 6; g.build()
    ao = dft.numint.eval_ao(mol, g.coords, deriv=1)
    rho = dft.numint.eval_rho(mol, ao, mf.make_rdm1(), xctype='GGA')
    n = np.maximum(rho[0], 1e-30); gr = np.sqrt((rho[1:4]**2).sum(0))
    s = gr/(2*(3*np.pi**2*n)**(1/3)*n); m = (s >= sstar) & (n > 1e-12)
    b = -K*np.sum((g.weights*n**(4/3)/s**2)[m])
    print(f"{name}: T4' bound E_c^sl <= {b:.4f} Ha   (near-exact E_c {ref}); fraction of density in s>=s*: {np.sum((g.weights*n)[m])/2:.3f}")
