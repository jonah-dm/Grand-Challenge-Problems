import numpy as np, sys
sys.path.insert(0,'../sonnet3')
from pyscf import gto, scf, dft
import importlib.util
spec=importlib.util.spec_from_file_location('s3','../sonnet3/submission.py'); s3=importlib.util.module_from_spec(spec); spec.loader.exec_module(s3)
def parts(atom,spin=0,basis='aug-cc-pvqz'):
    mol=gto.M(atom=atom,basis=basis,spin=spin,verbose=0)
    mf=scf.UHF(mol).run()
    g=dft.gen_grid.Grids(mol); g.level=6; g.build()
    ni=dft.numint.NumInt()
    dm=mf.make_rdm1()
    ao=ni.eval_ao(mol,g.coords,deriv=1)
    ra=ni.eval_rho(mol,ao,dm[0],xctype='MGGA'); rb=ni.eval_rho(mol,ao,dm[1],xctype='MGGA')
    w=g.weights
    na,nb=np.maximum(ra[0],1e-14),np.maximum(rb[0],1e-14)
    saa=np.sum(ra[1:4]**2,0); sbb=np.sum(rb[1:4]**2,0); sab=np.sum(ra[1:4]*rb[1:4],0)
    ta,tb=np.maximum(ra[-1],1e-20),np.maximum(rb[-1],1e-20)
    m=(na+nb)>1e-10
    n=na+nb; sig=saa+2*sab+sbb
    # total PBE corr via libxc
    tot=dft.libxc.eval_xc(',PBE',(ra[:4],rb[:4]),spin=1)[0]*n
    pa=s3._ec_pol(na[m],saa[m]); pb=s3._ec_pol(nb[m],sbb[m])
    za=s3._zsat(saa/(8*na*ta)); zb=s3._zsat(sbb/(8*nb*tb))
    Etot=np.dot(w,tot); Ess=np.dot(w[m],pa+pb)
    print(atom, 'Ec_PBE',Etot,'Ess',Ess,'Eos',Etot-Ess,'gated ss',np.dot(w[m],(1-za[m]**2)*pa+(1-zb[m]**2)*pb))
    return Etot,Ess
parts('He 0 0 0')
parts('H 0 0 0; H 0 0 0.741')
