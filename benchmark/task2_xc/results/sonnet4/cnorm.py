import numpy as np
from pyscf import gto, scf, dft
from pyscf.dft import libxc
m=gto.M(atom='He 0 0 0',basis='aug-cc-pv5z',verbose=0); mf=scf.RHF(m).run()
g=dft.gen_grid.Grids(m); g.level=9; g.build(); dm=mf.make_rdm1()
ao=dft.numint.eval_ao(m,g.coords,deriv=1); rho=dft.numint.eval_rho(m,ao,dm,xctype='GGA'); w=g.weights
n=rho[0]/2; gr=rho[1:4]/2; sg=np.einsum('ik,ik->k',gr,gr)
ok=n>1e-12
z=np.zeros_like(n)
def ec(ra,rb,ga,gb):
    e=libxc.eval_xc('GGA_C_PBE',(np.array([ra,ga,z,z]),np.array([rb,gb,z,z])),spin=1,deriv=0)[0]
    return (ra+rb)*e
ga=np.sqrt(sg)
Etot=np.sum(ec(n,n,ga,ga)*w*ok)
Ess=np.sum(ec(n,1e-14+z,ga,z)*w*ok)
print('PBE Ec He',Etot,'ss each',Ess,'opp',Etot-2*Ess)
print('c =',-0.04204/(Etot-2*Ess))
