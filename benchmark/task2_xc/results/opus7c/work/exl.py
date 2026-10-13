"""Add LR exact exchange at extra omegas to the cache (from stored UHF dm)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../kit')))
from pyscf import gto, scf
from geom import to_str
import dev_set as S
od = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache_cc-pvtz')
def get(name, w):
    fn = os.path.join(od, f'{name}_exl_{w:.6f}.npy')
    if os.path.exists(fn): return np.load(fn)
    atoms, charge, spin = S.SPECIES[name]
    mol = gto.M(atom=to_str(atoms), basis='cc-pvtz', charge=charge, spin=spin, verbose=0)
    dm = np.load(os.path.join(od, name + '_w.npz'))['dm']
    K = scf.hf.get_jk(mol, dm, hermi=1, with_j=False, omega=w)[1]
    e = np.array([-0.5*np.einsum('ij,ji', dm[s], K[s]) for s in range(2)])
    np.save(fn, e); return e
if __name__ == '__main__':
    for w in map(float, sys.argv[1:]):
        for nm in S.SPECIES: get(nm, w)
        print('done', w)
