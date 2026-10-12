import sys, numpy as np
sys.path.insert(0,'.'); sys.path.insert(0,'../../kit'); import submission as S
from pyscf import gto, dft
import dev_set as D; from geom import to_str
name=sys.argv[1]; basis=sys.argv[2] if len(sys.argv)>2 else 'def2-svp'
atoms,ch,sp=D.SPECIES[name]
mol=gto.M(atom=to_str(atoms),basis=basis,charge=ch,spin=sp,verbose=4)
mf=(dft.RKS if sp==0 else dft.UKS)(mol); mf.grids.level=4; mf.xc='HF'; mf=mf.define_xc_(S.eval_xc,'MGGA',hyb=S.HYB); mf.conv_tol=1e-8; mf.max_cycle=100; mf.run()
