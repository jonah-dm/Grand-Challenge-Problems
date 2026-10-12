import sys, numpy as np
sys.path.insert(0,'.'); import submission as S
from pyscf import gto, dft
mol=gto.M(atom='B 0 0 0',basis='def2-svp',charge=1,verbose=0)
mf=dft.RKS(mol); mf.xc='r2scan'; mf.grids.level=4; mf.run()
ni=mf._numint; g=mf.grids
ao=ni.eval_ao(mol,g.coords,deriv=1); rho=ni.eval_rho(mol,ao,mf.make_rdm1(),xctype='MGGA')
rho=rho[[0,1,2,3,5]]; print(rho.shape)
e,v=S.eval_xc('',rho,0)[:2]
E=lambda R:S.eval_xc('',R,0)[0]*R[0]
for idx,name,vv in [(0,'rho',v[0]),(4,'tau',v[3])]:
    d=1e-7*np.abs(rho[idx])+1e-12
    R=rho.copy(); R[idx]+=d; R2=rho.copy(); R2[idx]-=d
    fd=(E(R)-E(R2))/(2*d); err=np.abs(fd-vv)*g.weights
    i=np.argsort(err)[-5:]; print(name, err.sum(), [(rho[0][j],fd[j],vv[j]) for j in i])
# sigma via grad x
d=1e-7*np.abs(rho[1])+1e-12
R=rho.copy(); R[1]+=d; R2=rho.copy(); R2[1]-=d
fd=(E(R)-E(R2))/(2*d); an=2*v[1]*rho[1]; err=np.abs(fd-an)*g.weights
i=np.argsort(err)[-5:]; print('gx', err.sum(), [(rho[0][j],rho[1][j],fd[j],an[j]) for j in i])
r,sig,tau=rho[0],(rho[1:4]**2).sum(0),rho[4]
e,a,b,c=S._ex_and_derivs(r,sig,tau)
for nm,arr,an in [('r',0,a),('t',2,c)]:
    X=[r.copy(),sig.copy(),tau.copy()]; d=1e-6*X[arr]+1e-14
    Xp=[x.copy() for x in X]; Xp[arr]=Xp[arr]+d; Xm=[x.copy() for x in X]; Xm[arr]=Xm[arr]-d
    fd=(S._ex_and_derivs(*Xp)[0]-S._ex_and_derivs(*Xm)[0])/(2*d)
    err=np.abs(fd-an)*g.weights; i=np.argsort(err)[-3:]; print('X',nm,err.sum(),[(r[j],fd[j],an[j]) for j in i])
rc=rho.copy()
ec,vc=__import__('pyscf').dft.libxc.eval_xc('MGGA_C_R2SCAN',rho,0,0,1)[:2]
EC=lambda R:__import__('pyscf').dft.libxc.eval_xc('MGGA_C_R2SCAN',R,0,0,1)[0]*R[0]
for idx,vv in [(0,vc[0]),(4,vc[3])]:
    d=1e-6*np.abs(rho[idx])+1e-14
    R=rho.copy(); R[idx]+=d; R2=rho.copy(); R2[idx]-=d
    fd=(EC(R)-EC(R2))/(2*d); err=np.abs(fd-vv)*g.weights; i=np.argsort(err)[-3:]; print('C',idx,err.sum(),[(rho[0][j],fd[j],vv[j]) for j in i])
j=np.argmin(np.abs(r-2.7712621323283666e-09))
print('pt',r[j],sig[j],tau[j], 'tw',sig[j]/8/r[j])
for f in [0.999,0.9999,1,1.0001,1.001]:
    print(f, S._ex_unpol(np.array([r[j]]),np.array([sig[j]]),np.array([tau[j]*f]),S.A1))
