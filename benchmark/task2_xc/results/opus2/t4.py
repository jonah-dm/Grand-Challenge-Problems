import sys, numpy as np
sys.path.insert(0,'.'); sys.path.insert(0,'../../kit'); import submission as S
for kv in sys.argv[3:]:
    k,v=kv.split('='); setattr(S,k,float(v))
S.SLOPE = 0.5 * (1.0 + 2.0 * S.EPS) / (1.0 + S.EPS); S.MU1 = S.MU - (20.0 / 27.0) * S.KAPPA0 * S.SLOPE
from pyscf import gto, dft
import dev_set as D; from geom import to_str
name=sys.argv[1]; basis=sys.argv[2]
atoms,ch,sp=D.SPECIES[name]
mol=gto.M(atom=to_str(atoms),basis=basis,charge=ch,spin=sp,verbose=0)
mf=(dft.RKS if sp==0 else dft.UKS)(mol); mf.grids.level=4; mf.xc='HF'; mf=mf.define_xc_(S.eval_xc,'MGGA',hyb=S.HYB); mf.conv_tol=1e-8; mf.max_cycle=100; mf.run()
print(sys.argv[1:], mf.converged, mf.e_tot, mf.cycles if hasattr(mf,'cycles') else '')
