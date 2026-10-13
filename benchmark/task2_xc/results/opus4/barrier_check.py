"""Self-generated validation (not dev set): H + H2 -> H2 + H barrier, exact 9.6 kcal/mol.
Linear symmetric TS r=0.930 A, H2 r=0.741 A, aug-cc-pVTZ. Usage: env vars select candidate."""
import os, sys, importlib.util
from pyscf import gto, dft
spec=importlib.util.spec_from_file_location('s','submission.py'); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
def E(atom, spin):
    mol=gto.M(atom=atom,basis='aug-cc-pvtz',spin=spin,verbose=0)
    mf=(dft.RKS if spin==0 else dft.UKS)(mol); mf.grids.level=4; rsh=getattr(m,'RSH',(0,0,0))
    mf.xc='wB97X' if rsh[0] else 'HF'
    mf=mf.define_xc_(m.eval_xc,'MGGA',hyb=m.HYB,rsh=rsh); mf.conv_tol=1e-9; mf.run(); return mf.e_tot
eh=E('H 0 0 0',1); eh2=E('H 0 0 0; H 0 0 0.741',0); ets=E('H 0 0 -0.93; H 0 0 0; H 0 0 0.93',1)
print(os.environ.get('S3_GATE'), os.environ.get('O4_PW'), os.environ.get('O4_OMEGA'), 'barrier %.2f kcal/mol (exact 9.6)'%((ets-eh-eh2)*627.5095), 'H2 De %.2f (exact 109.5)'%((2*eh-eh2)*627.5095))
