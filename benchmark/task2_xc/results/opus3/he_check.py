# Validation (not used for any parameter): exchange of the He HF density with the z=1 (H-normed) form vs exact.
import numpy as np, sys
sys.path.insert(0, '.')
import submission as S
from pyscf import gto, scf, dft
mol = gto.M(atom='He', basis='aug-cc-pv5z', verbose=0)
mf = scf.RHF(mol).run()
dm = mf.make_rdm1()
vk = mf.get_k(dm); exx = -0.25 * np.einsum('ij,ji', dm, vk)
g = dft.gen_grid.Grids(mol); g.level = 6; g.build()
ao = dft.numint.eval_ao(mol, g.coords, deriv=1)
r = dft.numint.eval_rho(mol, ao, dm, xctype='MGGA')
n, sig, tau = r[0], (r[1:4]**2).sum(0), r[-1]
m = n > 1e-12
for name, mu1 in [('H-normed mu1', S.MU1), ('PBE mu0', S.MU0)]:
    ex = S._ex_unpol(n[m], sig[m], np.maximum(tau[m], sig[m]/(8*n[m])), mu1)
    print(name, 'Ex =', (ex * g.weights[m]).sum(), ' exact HF =', exx)
