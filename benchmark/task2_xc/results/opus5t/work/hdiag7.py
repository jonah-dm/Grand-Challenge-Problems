exec(open('hdiag.py').read().split('mf2=')[0])
R=rho.copy(); w=g.weights
exc,vxc=mz.eval_xc('',R,1,0,1,None,0)[:2]
def E(R):
    e=mz.eval_xc('',R,1,0,1,None,0)[0]; return e*(R[0,0]+R[1,0])
for nm,row,v in (('rho',0,vxc[0][:,0]),('tau',5,vxc[3][:,0])):
    for h in (1e-4,1e-6,1e-8):
        Rp=R.copy(); Rp[0,row]*=1+h; Rm=R.copy(); Rm[0,row]*=1-h
        fd=(E(Rp)-E(Rm))/(2*h*R[0,row]+1e-300)
        mm=np.abs((fd-v)*R[0,row]*w)
        print(nm,h,'weighted mismatch',mm.sum(),'worst r',r[np.argmax(mm)], R[0,0,np.argmax(mm)], R[0,5,np.argmax(mm)], v[np.argmax(mm)], fd[np.argmax(mm)])
