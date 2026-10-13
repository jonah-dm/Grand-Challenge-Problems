"""Optimal tuning of omega from the exact IP theorem for the H atom: -eps_HOMO(H) = 0.5 Ha."""
import os, sys, importlib.util, numpy as np
from pyscf import gto, dft
def run(om, basis='aug-cc-pvqz', atom='H 0 0 0', spin=1):
    os.environ['O4_OMEGA']=str(om)
    spec=importlib.util.spec_from_file_location('s','submission.py'); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    mol=gto.M(atom=atom,basis=basis,spin=spin,verbose=0)
    mf=dft.UKS(mol); mf.grids.level=4
    rsh=getattr(m,'RSH',(0,0,0))
    mf.xc='wB97X' if rsh[0] else 'HF'
    mf=mf.define_xc_(m.eval_xc,'MGGA',hyb=m.HYB,rsh=rsh); mf.conv_tol=1e-10; mf.run()
    ea=mf.mo_energy[0][mf.mo_occ[0]>0].max()
    return mf.e_tot, ea
ATOM=os.environ.get('TUNE_ATOM','H 0 0 0'); SP=int(os.environ.get('TUNE_SPIN','1'))
for om in [float(x) for x in sys.argv[1:]]:
    e,h=run(om,atom=ATOM,spin=SP); print(f'omega={om:.4f} E(H)={e:.6f}  -eps_HOMO={-h:.5f}', flush=True)
