import os, importlib.util, numpy as np
from pyscf import gto, dft
def load(env):
    os.environ['O5_W']=env
    spec=importlib.util.spec_from_file_location('s'+env,'../submission.py'); m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
mz, mo = load('zeta2'), load('off')
mol=gto.M(atom='H 0 0 0',basis='aug-cc-pvtz',spin=1,verbose=0)
mf=dft.UKS(mol); mf.grids.level=4; mf.xc='HF'; mf=mf.define_xc_(mo.eval_xc,'MGGA',hyb=mo.HYB); mf.run()
dm=mf.make_rdm1(); ni=mf._numint; g=mf.grids
ao=ni.eval_ao(mol,g.coords,deriv=1)
ra=ni.eval_rho(mol,ao,dm[0],xctype='MGGA'); rb=ni.eval_rho(mol,ao,dm[1],xctype='MGGA')
rho=np.array([ra,rb]); print(rho.shape)
for nm,m in (('off',mo),('zeta2',mz)):
    exc=m.eval_xc('',rho,1,0,1,None,0)[0]
    print(nm, np.sum(exc*(ra[0]+rb[0])*g.weights))
# far-tail contributions
r=np.linalg.norm(g.coords,axis=1)
exc=mz.eval_xc('',rho,1,0,1,None,0)[0]; exo=mo.eval_xc('',rho,1,0,1,None,0)[0]
for rc in (2,4,6,8,10,15):
    sel=r>rc; print(rc, np.sum((exc*(ra[0]+rb[0])*g.weights)[sel]), np.sum((exo*(ra[0]+rb[0])*g.weights)[sel]))
mf2=dft.UKS(mol); mf2.grids.level=4; mf2.xc='HF'; mf2=mf2.define_xc_(mz.eval_xc,'MGGA',hyb=mz.HYB); mf2.verbose=4
mf2.max_cycle=30
mf2.kernel(dm0=dm)
