exec(open('hdiag.py').read().split('mf2=')[0])
S=mol.intor('int1e_ovlp'); w=g.weights
c0=mf.mo_coeff[0][:,0]+0.05*mf.mo_coeff[0][:,1]; c0/=np.sqrt(c0@S@c0)
d=np.random.default_rng(0).normal(size=c0.size)*0.01
def rh(c):
    c=c/np.sqrt(c@S@c); return np.array([ni.eval_rho(mol,ao,np.outer(c,c),xctype='MGGA'),ni.eval_rho(mol,ao,0*S,xctype='MGGA')])
h=1e-4
Rp,R0,Rm=rh(c0+h*d),rh(c0),rh(c0-h*d)
def e(R): return mz.eval_xc('',R,1,0,1,None,0)[0]*(R[0,0]+R[1,0])
num=(e(Rp)-e(Rm))/(2*h)*w
exc,vxc=mz.eval_xc('',R0,1,0,1,None,0)[:2]
dR=(Rp-Rm)/(2*h)
lin=(vxc[0][:,0]*dR[0,0]+2*vxc[1][:,0]*np.sum(R0[0,1:4]*dR[0,1:4],0)+vxc[3][:,0]*dR[0,5])*w
print('sum num',num.sum(),'sum lin',lin.sum())
diff=np.abs(num-lin); o=np.argsort(-diff)[:8]
print(np.c_[r[o],diff[o],R0[0,0,o],np.sum(R0[0,1:4]**2,0)[o],R0[0,5,o],vxc[3][o,0]])
