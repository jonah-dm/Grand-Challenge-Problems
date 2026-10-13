import sys, os, numpy as np
here=os.path.dirname(os.path.abspath(__file__)); root=os.path.abspath(os.path.join(here,'..','..','..'))
sys.path.insert(0,root); sys.path.insert(0,os.path.join(root,'kit'))
import score
from pyscf import gto, dft
sub=score.load_submission(os.path.join(here,'..','submission.py')); print('cand', os.environ.get('O7_K','S5'), sub.OMEGA, sub.F_LO)
basis=sys.argv[1] if len(sys.argv)>1 else 'def2-qzvp'
mol=gto.M(atom='H 0 0 0',basis=basis,spin=1,verbose=0)
mf=dft.UKS(mol); mf.grids.level=4
mf.xc='wB97X'; mf=mf.define_xc_(sub.eval_xc,'MGGA',hyb=sub.HYB,rsh=sub.RSH); mf.run()
print('E',mf.e_tot)
r=np.array([0,0.01,0.05,0.1,0.2,0.4,0.7,1,1.5,2,3,4.])
c=np.zeros((len(r),3)); c[:,2]=r
ao=dft.numint.eval_ao(mol,c,deriv=1); dm=mf.make_rdm1()[0]
rho=dft.numint.eval_rho(mol,ao,dm,xctype='GGA')
n=rho[0]; g=np.linalg.norm(rho[1:4],axis=0); nex=np.exp(-2*r)/np.pi
b=g**2/(4*sub.OMEGA**2*n**2); k=(6*np.pi**2*n)**(1/3); F=5/9*(k/sub.OMEGA)**2*sub.Qfun(b)
for i in range(len(r)): print('%5.2f n %.4f exact %.4f  b %.3f (H exact %.3f) F %.3f' % (r[i],n[i],nex[i],b[i],1/sub.OMEGA**2,F[i]))
