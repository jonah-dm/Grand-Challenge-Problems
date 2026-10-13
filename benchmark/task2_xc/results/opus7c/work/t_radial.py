"""Radial decomposition of the SR exchange IP error: exact conventional SR-HF energy density vs semilocal SR (opus6c) energy density."""
import sys, os, numpy as np
sys.path.insert(0, os.path.abspath('../../../kit'))
from pyscf import gto, dft
from geom import to_str
import dev_set as S
from diag import parts, M
H2K = 627.5094740631; W = M.OMEGA
od = 'cache_cc-pvtz'
bins = np.array([0, 0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 5.0, 50])
def run(nm):
    atoms, charge, spin = S.SPECIES[nm]
    mol = gto.M(atom=to_str(atoms), basis='cc-pvtz', charge=charge, spin=spin, verbose=0)
    dm = np.load(os.path.join(od, nm + '_w.npz'))['dm']
    g = dft.gen_grid.Grids(mol); g.level = 3; g.build()
    ni = dft.numint.NumInt(); ao = ni.eval_ao(mol, g.coords, deriv=1)
    rho = [ni.eval_rho(mol, ao, dm[s], xctype='MGGA', with_lapl=False) for s in range(2)]
    ra = np.maximum(rho[0][0], 1e-14); rb = np.maximum(rho[1][0], 1e-14); ga, gb = rho[0][1:4], rho[1][1:4]
    with np.errstate(all='ignore'):
        P = parts(ra, rb, (ga**2).sum(0), (ga*gb).sum(0), (gb**2).sum(0), np.maximum(rho[0][4], 1e-20), np.maximum(rho[1][4], 1e-20))
    esl = np.nan_to_num(sum((1-P['w_'+l])*P['base_'+l] + P['w_'+l]*P['br_'+l] for l in 'ab'))
    ex = np.zeros(len(g.weights)); chi = ao[0]
    for p0 in range(0, len(g.weights), 2000):
        c = g.coords[p0:p0+2000]
        with mol.with_range_coulomb(-W):
            V = mol.intor('int1e_grids', grids=c)          # (ng, nao, nao) SR
        for s in range(2):
            Dc = chi[p0:p0+2000] @ dm[s]                         # (ng, nao)
            ex[p0:p0+2000] += -0.5*np.einsum('gi,gij,gj->g', Dc, V, Dc)
    r = np.linalg.norm(g.coords - mol.atom_coords()[0], axis=1)
    idx = np.digitize(r, bins) - 1
    out = np.array([[np.sum((g.weights*ex)[idx == i]), np.sum((g.weights*esl)[idx == i])] for i in range(len(bins)-1)])
    return out
res = {nm: run(nm) for nm in sys.argv[1:]}
for a, b in [(sys.argv[i], sys.argv[i+1]) for i in range(1, len(sys.argv), 2)]:
    print(f'{b} - {a}: per radial bin (bohr) {bins.tolist()}')
    dexact = (res[b][:, 0]-res[a][:, 0])*H2K; dsl = (res[b][:, 1]-res[a][:, 1])*H2K
    print('  d exact SR-HF   ', np.round(dexact, 1), 'sum', round(dexact.sum(), 2))
    print('  d semilocal SR  ', np.round(dsl, 1), 'sum', round(dsl.sum(), 2))
    print('  IP error (sl-ex)', np.round(dsl-dexact, 1), 'sum', round((dsl-dexact).sum(), 2))
