"""Cache UHF densities on grids for all dev species (non-SCF evaluation of functional pieces)."""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../kit')))
from pyscf import gto, scf, dft, lib
from geom import to_str
import dev_set as S
W = 0.5767436061911404
basis = sys.argv[1] if len(sys.argv) > 1 else 'cc-pvtz'
od = os.path.join(os.path.dirname(__file__), 'cache_' + basis); os.makedirs(od, exist_ok=True)
for name, (atoms, charge, spin) in S.SPECIES.items():
    fn = os.path.join(od, name + '.npz')
    if os.path.exists(fn): continue
    mol = gto.M(atom=to_str(atoms), basis=basis, charge=charge, spin=spin, verbose=0)
    mf = scf.UHF(mol); mf.conv_tol = 1e-9; mf.max_cycle = 200; mf.run()
    if not mf.converged: mf = scf.newton(mf).run()
    dm = mf.make_rdm1()
    Kf = mf.get_k(mol, dm); Kl = mf.get_k(mol, dm, omega=W)
    exs = [-0.5*np.einsum('ij,ji', dm[s], Kf[s]) for s in range(2)]
    exl = [-0.5*np.einsum('ij,ji', dm[s], Kl[s]) for s in range(2)]
    g = dft.gen_grid.Grids(mol); g.level = 4; g.build()
    ni = dft.numint.NumInt()
    R = []
    for s in range(2):
        rr = []
        for ao, mask, wt, coords in ni.block_loop(mol, g, mol.nao, deriv=1):
            rr.append(ni.eval_rho(mol, ao, dm[s], xctype='MGGA', with_lapl=False))
        R.append(np.concatenate(rr, axis=1))
    np.savez(fn, rho=np.array(R), w=g.weights, coords=g.coords, E_hf=mf.e_tot, Ex=np.array(exs), Exl=np.array(exl),
             s2=mf.spin_square()[0], atom0=mol.atom_coords()[0])
    print(name, mf.e_tot, flush=True)
