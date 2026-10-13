"""H-atom norm imposed at the self-consistent density (def2-QZVP, grid level 4, as the harness):
find C_Q such that E_SCF(H) = -0.5 exactly. Not dev-set data."""
import os, importlib.util
from scipy.optimize import brentq
from pyscf import gto, dft
spec=importlib.util.spec_from_file_location('s','../submission.py'); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
mol=gto.M(atom='H 0 0 0',basis='def2-qzvp',spin=1,verbose=0)
def EH(c):
    m.C_Q=c
    mf=dft.UKS(mol); mf.grids.level=4; mf.xc='HF'; mf=mf.define_xc_(m.eval_xc,'MGGA',hyb=m.HYB); mf.conv_tol=1e-10; mf.run()
    print(c, mf.e_tot, flush=True); return mf.e_tot+0.5
c=brentq(EH,1.2,2.2088,xtol=1e-5); print('C_Q_SCF',repr(c))
