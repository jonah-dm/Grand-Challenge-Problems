# Generate He Hartree-Fock density + exact exchange (our own norm data).
import numpy as np
from pyscf import gto, scf, dft
mol = gto.M(atom="He 0 0 0", basis="aug-cc-pv5z", verbose=0)
mf = scf.RHF(mol).run()
dm = mf.make_rdm1()
vk = mf.get_k(dm)
ex = -0.25 * np.einsum("ij,ji", dm, vk)
r = np.linspace(1e-4, 15, 30001)
coords = np.zeros((r.size, 3)); coords[:, 2] = r
ao = dft.numint.eval_ao(mol, coords, deriv=1)
rho = dft.numint.eval_rho(mol, ao, dm, xctype="GGA")
np.savez("he_hf.npz", r=r, rho=rho[0], drho=rho[3], ex=ex, e=mf.e_tot)
print(mf.e_tot, ex, np.trapezoid(rho[0] * 4 * np.pi * r * r, r))
