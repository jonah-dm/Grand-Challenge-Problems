"""Second-row transfer check (not used for selection of constants): CCSD(T)/cc-pVTZ (fc) reference + UHF/cc-pVTZ density proxy."""
import sys, os, json, numpy as np
sys.path.insert(0, os.path.abspath('../../../kit'))
from pyscf import gto, scf, cc, dft
from geom import to_str, diatomic, bent, pyramid, tetra
SP = {'H_atom': ([('H', (0, 0, 0))], 0, 1), 'P_atom': ([('P', (0, 0, 0))], 0, 3), 'S_atom': ([('S', (0, 0, 0))], 0, 2),
      'Cl_atom': ([('Cl', (0, 0, 0))], 0, 1), 'Si_atom': ([('Si', (0, 0, 0))], 0, 2),
      'HCl': (diatomic('H', 'Cl', 1.275), 0, 0), 'H2S': (bent('S', 'H', 1.336, 92.1), 0, 0), 'PH3': (pyramid('P', 'H', 1.420, 93.3), 0, 0),
      'SiH4': (tetra('Si', 'H', 1.480), 0, 0), 'Cl2': (diatomic('Cl', 'Cl', 1.988), 0, 0),
      'Cl_cation': ([('Cl', (0, 0, 0))], 1, 2), 'S_cation': ([('S', (0, 0, 0))], 1, 3)}
od = 'cache_row2'; os.makedirs(od, exist_ok=True)
W = [0.5767436061911404, 0.5256070335369912]
NC = {1: 0, 14: 5, 15: 5, 16: 5, 17: 5}
for nm, (atoms, ch, sp) in SP.items():
    fn = os.path.join(od, nm + '.npz')
    if os.path.exists(fn): continue
    mol = gto.M(atom=to_str(atoms), basis='cc-pvtz', charge=ch, spin=sp, verbose=0)
    mf = scf.UHF(mol); mf.conv_tol = 1e-10; mf.max_cycle = 300; mf.run()
    if not mf.converged: mf = scf.newton(mf).run()
    nc = sum(NC[gto.charge(a)] for a, _ in atoms)
    ec = 0.0
    if mol.nelectron > 1:
        m = cc.CCSD(mf, frozen=nc if nc else None); m.run(); ec = m.e_corr + m.ccsd_t()
    dm = mf.make_rdm1()
    Ex = [-0.5*np.einsum('ij,ji', dm[s], mf.get_k(mol, dm)[s]) for s in range(2)]
    Exl = [[-0.5*np.einsum('ij,ji', dm[s], K[s]) for s in range(2)] for K in [mf.get_k(mol, dm, omega=w) for w in W]]
    g = dft.gen_grid.Grids(mol); g.level = 4; g.build(); ni = dft.numint.NumInt()
    R = []
    for s in range(2):
        rr = [ni.eval_rho(mol, ao, dm[s], xctype='MGGA', with_lapl=False) for ao, mask, wt, c in ni.block_loop(mol, g, mol.nao, deriv=1)]
        R.append(np.concatenate(rr, axis=1))
    np.savez(fn, rho=np.array(R), w=g.weights, E_hf=mf.e_tot, Ex=np.array(Ex), ExlW=np.array(Exl), ws=np.array(W), ec=ec)
    print(nm, mf.e_tot, ec, flush=True)
