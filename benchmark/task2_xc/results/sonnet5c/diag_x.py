"""Exact-exchange partition diagnostic (non-dev species): HF density, exact HF exchange vs semilocal exchange
forms evaluated on the same density. Usage: python diag_x.py [basis]"""
import sys, pickle, numpy as np
from pyscf import gto, scf, dft
from pyscf.dft import numint
import importlib, os
BASIS = sys.argv[1] if len(sys.argv) > 1 else 'def2-qzvp'
sys.path.insert(0, '/home/user/Grand-Challenge-Problems/benchmark/task2_xc/kit')
from geom import diatomic, bent, pyramid, tetra, linear
# non-dev training species: (atoms, charge, spin)
SP = {
 'Be':([('Be',(0,0,0))],0,0), 'Ne':([('Ne',(0,0,0))],0,0), 'Mg':([('Mg',(0,0,0))],0,0), 'Ar':([('Ar',(0,0,0))],0,0),
 'Na':([('Na',(0,0,0))],0,1), 'Si':([('Si',(0,0,0))],0,2), 'P':([('P',(0,0,0))],0,3), 'S':([('S',(0,0,0))],0,2), 'Cl':([('Cl',(0,0,0))],0,1),
 'He':([('He',(0,0,0))],0,0),
 'HCl':(diatomic('Cl','H',1.275),0,0), 'H2S':(bent('S','H',1.336,92.1),0,0), 'PH3':(pyramid('P','H',1.42,93.3),0,0),
 'SiH4':(tetra('Si','H',1.48),0,0), 'LiF':(diatomic('Li','F',1.564),0,0), 'BeH2':(linear(('H',0),('Be',1.334),('H',1.334)),0,0),
 'NaCl':(diatomic('Na','Cl',2.361),0,0), 'Li2':(diatomic('Li','Li',2.673),0,0), 'BH':(diatomic('B','H',1.232),0,0),
 'HeH+':(diatomic('He','H',0.774),1,0),
}
def run(name):
    atoms, ch, sp = SP[name]
    mol = gto.M(atom=atoms, charge=ch, spin=sp, basis=BASIS, verbose=0)
    mf = scf.UHF(mol) if sp else scf.RHF(mol)
    mf.conv_tol = 1e-10; mf.kernel()
    dm = mf.make_rdm1()
    g = dft.gen_grid.Grids(mol); g.level = 5; g.build()
    if sp:
        dma, dmb = dm; ka = mf.get_k(mol, dma); kb = mf.get_k(mol, dmb)
        Ex = -0.5*(np.einsum('ij,ji', dma, ka) + np.einsum('ij,ji', dmb, kb))
    else:
        k = mf.get_k(mol, dm); Ex = -0.25*np.einsum('ij,ji', dm, k)
        dma = dmb = dm/2
    ni = numint.NumInt()
    ao = ni.eval_ao(mol, g.coords, deriv=1)
    ra = ni.eval_rho(mol, ao, dma, xctype='MGGA', with_lapl=False)
    rb = ni.eval_rho(mol, ao, dmb, xctype='MGGA', with_lapl=False)
    return dict(name=name, E=mf.e_tot, Ex=Ex, w=g.weights, ra=ra, rb=rb)
if __name__ == '__main__':
    out = {n: run(n) for n in SP}
    pickle.dump(out, open('diag_x_%s.pkl' % BASIS, 'wb'))
    print({n: (round(v['E'], 5), round(v['Ex'], 5)) for n, v in out.items()})
