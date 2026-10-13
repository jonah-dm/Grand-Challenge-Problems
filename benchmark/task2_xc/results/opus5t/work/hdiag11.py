exec(open('hdiag.py').read().split('mf2=')[0])
S=mol.intor('int1e_ovlp')
c0=mf.mo_coeff[0][:,0]+0.05*mf.mo_coeff[0][:,1]; c0/=np.sqrt(c0@S@c0)
rec=[]
def ev(*a,**k):
    rec.append(np.array(a[1]).copy()); return mz.eval_xc(*a,**k)
m2=dft.UKS(mol); m2.grids.level=4; m2.xc='HF'; m2=m2.define_xc_(ev,'MGGA',hyb=0.25); m2.grids.build()
D=np.array([np.outer(c0,c0),0*S]); out=m2._numint.nr_uks(mol,m2.grids,m2.xc,D)
R=np.concatenate(rec,axis=-1); print(R.shape, 'min rho_b', R[1,0].min(), R[1,0].max(), 'rows', [float(abs(R[0,i]).max()) for i in range(R.shape[1])])
ao=ni.eval_ao(mol,m2.grids.coords,deriv=1); Rm=np.array([ni.eval_rho(mol,ao,D[0],xctype='MGGA'),ni.eval_rho(mol,ao,D[1],xctype='MGGA')])
print('my rows', [float(abs(Rm[0,i]).max()) for i in range(Rm.shape[1])], 'npts', Rm.shape[-1])
w=m2.grids.weights
exc,vxc=mz.eval_xc('',Rm,1,0,1,None,0)[:2]
vr,vs,vt=vxc[0][:,0],vxc[1][:,0],vxc[3][:,0]
phi=ao[0]; dphi=ao[1:4]
V1=np.einsum('g,gm,gn->mn',w*vr,phi,phi)
gr=Rm[0,1:4]
t=np.einsum('xg,xgm->gm',gr,dphi)
V2=2*(np.einsum('g,gm,gn->mn',w*vs,t,phi)+np.einsum('g,gm,gn->mn',w*vs,phi,t))
V3=0.5*np.einsum('g,xgm,xgn->mn',w*vt,dphi,dphi)
V=V1+V2+V3; Vp=out[2][0]
print('max diff', abs(V-Vp).max(), 'parts', abs(V1).max(),abs(V2).max(),abs(V3).max())
for nm,X in (('V1',V1),('V2',V2),('V3',V3)): print(nm, abs(X-Vp+ (V-X)).max())
print('V1+V2 vs Vp', abs(V1+V2-Vp).max(), ' V1+V2+2V3', abs(V1+V2+2*V3-Vp).max(), 'V1+V2+0.5V3',abs(V1+V2+0.5*V3-Vp).max())
print('vt range', vt.min(), vt.max(), 'argmax r', np.linalg.norm(m2.grids.coords,axis=1)[np.argmax(abs(vt))])
