"""Correlation energies of He, H2 (HF densities, aug-cc-pVQZ) and H for the gate variants."""
import os, importlib.util, numpy as np
from pyscf import gto, scf, dft
def load(gate,pw='2'):
    os.environ['S3_GATE']=gate; os.environ['O4_PW']=pw; os.environ['O4_OMEGA']='0'
    spec=importlib.util.spec_from_file_location('s','submission.py'); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); m.HYB=1.0
    return m
def ec(m, atom, spin=0):
    mol=gto.M(atom=atom,basis='aug-cc-pvqz',spin=spin,verbose=0)
    mf=scf.UHF(mol).run(); g=dft.gen_grid.Grids(mol); g.level=6; g.build()
    ni=dft.numint.NumInt(); dm=mf.make_rdm1(); ao=ni.eval_ao(mol,g.coords,deriv=1)
    rho=np.array([ni.eval_rho(mol,ao,dm[s],xctype='MGGA') for s in (0,1)])
    exc=m.eval_xc('',rho,1)[0]
    return np.dot(g.weights, exc*(rho[0,0]+rho[1,0]))
for gate,pw in [('none','2'),('stoll','2'),('pairz','2'),('pairz','1'),('pairt','2')]:
    m=load(gate,pw)
    print(gate,pw,'H %.5f He %.5f H2 %.5f'%(ec(m,'H 0 0 0',1),ec(m,'He 0 0 0'),ec(m,'H 0 0 0; H 0 0 0.741')),flush=True)
print('exact: H 0, He -0.04204, H2 -0.0409')
