exec(open('hdiag.py').read().split('mf2=')[0])
R=rho.copy()
exc,vxc=mz.eval_xc('',R,1,0,1,None,0)[:2]
vs=vxc[1][:,0]
def E(R):
    e=mz.eval_xc('',R,1,0,1,None,0)[0]; return e*(R[0,0]+R[1,0])
# perturb sigma_aa by scaling gradient components: sigma -> sigma(1+2h)
h=1e-6
Rp=R.copy(); Rp[0,1:4]*=1+h; Rm=R.copy(); Rm[0,1:4]*=1-h
sig=np.sum(R[0,1:4]**2,0)
fd=(E(Rp)-E(Rm))/(2*h*2*sig+1e-300)
bad=np.abs(fd-vs)>1e-4*np.abs(vs)+1e-8
print('n bad',bad.sum(),'of',len(vs))
w=g.weights
print('weighted mismatch', np.sum(np.abs((fd-vs)*2*sig*w)[bad]))
i=np.where(bad)[0][:10]; print(np.c_[r[i],R[0,0,i],sig[i],vs[i],fd[i]])
