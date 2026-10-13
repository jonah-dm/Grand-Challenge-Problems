import numpy as np, pickle
from scipy.optimize import least_squares
allr=pickle.load(open('lib2.pkl','rb'))
sel=[0,2,3,4,5,6,9,10]; names=['H','H2+1','2','3','4.5','7','HeH2','HeH3.5']
D=[]
for i in sel:
    r=allr[i]; m=r['n']>1e-10
    n=r['n'][m]; sig=r['sig'][m]; lap=r['lap'][m]
    n2=2*n; sig2=4*sig
    kf=(3*np.pi**2*n2)**(1/3); s=np.sqrt(sig2)/(2*kf*n2)
    w=n*lap/np.maximum(sig,1e-30)
    D.append((r['w'][m],-0.75*(3/np.pi)**(1/3)*n2**(4/3)*0.5,s,w,-r['J']))
tar=np.array([d[4] for d in D])
def Ex(G):
    return np.array([np.sum(wt*exl*G(s,w)/np.maximum(s,1e-6)**2) for wt,exl,s,w,_ in D])
def mk(p):
    def G(s,w):
        t=w/(1+np.abs(w))
        return np.exp(p[0]+p[1]*t+p[2]*t**2+p[3]*t**3)
    return G
best=None
for p0 in ([-1,0,0,0],[-0.5,-1,0,0],[-1,1,0,0],[0,0,-1,0]):
    sol=least_squares(lambda p:Ex(mk(p))/tar-1,p0)
    if best is None or sol.cost<best.cost: best=sol
print(best.x,best.cost)
print(names); print(np.round(100*(Ex(mk(best.x))/tar-1),1))
print('--- held-out test: fit on H, H2+ R=1,3,7 ; test others')
idx=[0,1,3,5]
def fitidx(idx,deg=3):
    def mkp(p):
        def G(s,w):
            t=w/(1+np.abs(w)); return np.exp(sum(p[k]*t**k for k in range(len(p))))
        return G
    best=None
    for p0 in ([-1,0,0,0],[-0.5,-1,0,0],[-1,1,0,0]):
        p0=p0[:deg+1]
        sol=least_squares(lambda p:(Ex(mkp(p))/tar-1)[idx],p0)
        if best is None or sol.cost<best.cost: best=sol
    return best.x,mkp(best.x)
for deg in (2,3):
    p,G=fitidx(idx,deg); print('deg',deg,np.round(p,3)); print(np.round(100*(Ex(G)/tar-1),1))
