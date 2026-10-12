import sys, numpy as np
sys.path.insert(0,'.'); import submission as S
from pyscf import gto, dft
# derivative check
rng=np.random.default_rng(1); n=5
rho=np.zeros((2,5,n)); 
for k in range(2):
    rho[k,0]=10**rng.uniform(-3,0,n); g=rng.normal(size=(3,n))*rho[k,0]; rho[k,1:4]=g
    rho[k,4]=(g**2).sum(0)/8/rho[k,0]+rng.uniform(0,1,n)*rho[k,0]**(5/3)
e,v=S.eval_xc('',rho,1)[:2]
E=lambda R:(S.eval_xc('',R,1)[0]*(R[0,0]+R[1,0]))
d=1e-6
for k in range(2):
    R=rho.copy(); R[k,0]+=d; print('vrho',k,(E(R)-E(rho))/d, v[0][:,k])
    R=rho.copy(); R[k,4]+=d; print('vtau',k,(E(R)-E(rho))/d, v[3][:,k])
    R=rho.copy(); R[k,1]+=d; print('vsig',k,(E(R)-E(rho))/d, v[1][:,2*k]*2*rho[k,1]+v[1][:,1]*rho[1-k,1])
mol=gto.M(atom='H 0 0 0',basis='def2-svp',spin=1,verbose=3)
mf=dft.UKS(mol); mf.grids.level=4; mf.xc='HF'; mf=mf.define_xc_(S.eval_xc,'MGGA',hyb=S.HYB); mf.run()
