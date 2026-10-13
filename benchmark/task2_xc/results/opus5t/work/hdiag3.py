exec(open('hdiag.py').read().split('mf2=')[0])
idx=np.argsort(r)[::max(1,len(r)//40)]
R=rho[:,:,idx].copy()
exc,vxc=mz.eval_xc('',R,1,0,1,None,0)[:2]
def E(R): 
    e=mz.eval_xc('',R,1,0,1,None,0)[0]; return e*(R[0,0]+R[1,0])
va=vxc[0][:,0]; vs=vxc[1][:,0]; vt=vxc[3][:,0]
h=1e-6
Rp=R.copy(); Rp[0,0]*=1+h; Rm=R.copy(); Rm[0,0]*=1-h
fd=(E(Rp)-E(Rm))/(2*h*R[0,0])
Rp=R.copy(); Rp[0,5]*=1+h; Rm=R.copy(); Rm[0,5]*=1-h
fdt=(E(Rp)-E(Rm))/(2*h*R[0,5])
for i in range(len(idx)): print('%6.2f n=%.2e vra %.4e fd %.4e | vta %.4e fd %.4e'%(r[idx][i],R[0,0,i],va[i],fd[i],vt[i],fdt[i]))
