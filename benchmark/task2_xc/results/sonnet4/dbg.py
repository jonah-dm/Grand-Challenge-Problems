import sys, importlib.util, numpy as np
from pyscf import gto, dft
spec=importlib.util.spec_from_file_location("s","submission.py"); s=importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
mol=gto.M(atom='H 0 0 0',basis='def2-svp',spin=1,verbose=4)
mf=dft.UKS(mol); mf.grids.level=4; mf.xc='HF'
mf=mf.define_xc_(s.eval_xc,'MGGA',hyb=0.25); mf.max_cycle=40; mf.run()
