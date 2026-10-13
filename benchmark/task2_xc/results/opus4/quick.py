import os, sys, importlib.util, time
sys.path.insert(0,'../../kit')
import dev_set
from pyscf import gto, dft
spec=importlib.util.spec_from_file_location('s','submission.py'); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
for name in sys.argv[2:]:
    atoms,ch,sp=dev_set.SPECIES[name]
    mol=gto.M(atom=[(a,c) for a,c in atoms],basis=sys.argv[1],charge=ch,spin=sp,verbose=0)
    mf=(dft.RKS if sp==0 else dft.UKS)(mol); mf.grids.level=4; mf.xc='HF'
    mf=mf.define_xc_(m.eval_xc,'MGGA',hyb=m.HYB); mf.conv_tol=1e-8; t=time.time(); mf.run()
    print(name, mf.e_tot, mf.converged, mf.cycles if hasattr(mf,'cycles') else '', round(time.time()-t,1), flush=True)
