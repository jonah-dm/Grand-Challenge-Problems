import sys, os
import numpy as np
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../kit')))
from pyscf import gto, scf
from geom import to_str
import dev_set as S
WS = [0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 1.2]
od = os.path.join(os.path.dirname(__file__), 'cache_cc-pvtz')
for name, (atoms, charge, spin) in S.SPECIES.items():
    fn = os.path.join(od, name + '_w.npz')
    if os.path.exists(fn): continue
    mol = gto.M(atom=to_str(atoms), basis='cc-pvtz', charge=charge, spin=spin, verbose=0)
    mf = scf.UHF(mol); mf.conv_tol = 1e-9; mf.max_cycle = 200; mf.run()
    if not mf.converged: mf = scf.newton(mf).run()
    dm = mf.make_rdm1()
    E = []
    for w in WS:
        K = mf.get_k(mol, dm, omega=w); E.append([-0.5*np.einsum('ij,ji', dm[s], K[s]) for s in range(2)])
    np.savez(fn, ws=np.array(WS), Exl=np.array(E), dm=dm, mo_occ=mf.mo_occ, mo_e=mf.mo_energy)
    print(name, flush=True)
