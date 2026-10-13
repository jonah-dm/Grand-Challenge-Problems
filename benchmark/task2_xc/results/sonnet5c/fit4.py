import numpy as np, pickle
from scipy.optimize import least_squares
allr=[r for r in pickle.load(open('lib.pkl','rb')) if r['k']==0]
sel=[0,2,3,4,5,6,9,10]
D=[]
for i in sel:
    r=allr[i]
    n=np.maximum(r['n'],1e-30); n2=2*n; sig2=4*r['sig']
    kf=(3*np.pi**2*n2)**(1/3); s=np.sqrt(sig2)/(2*kf*n2)
    exl=-0.75*(3/np.pi)**(1/3)*n2**(4/3)*0.5
    D.append((r['w'],exl,s,-r['J']))
tar=np.array([d[3] for d in D])
def Ex(F,lam=1.0): return np.array([lam**(4/3)*np.sum(w*exl*F(s*lam**(-1/3))) for w,exl,s,_ in D])
def mk(p):
    a,b,c,d=p
    return lambda s:1+(a*s**2+d*s**4)/(1+b*s**4+c*s**6)
lams=[1.0,0.75,0.5]
def res(p,wl=(3,1,1)):
    return np.concatenate([w*(Ex(mk(p),l)/(l**2*tar)-1) for l,w in zip(lams,wl)])
best=None
rng=np.random.default_rng(0)
for t in range(40):
    p0=rng.uniform(0.01,3,4)*np.array([1,1,1,1])
    try: sol=least_squares(res,p0,bounds=([0,0,1e-6,0],[100,1e3,1e3,1e3]))
    except Exception as e: continue
    if best is None or sol.cost<best.cost: best=sol
p=best.x; F=mk(p); print(p,best.cost)
for l in lams: print(l,np.round(100*(Ex(F,l)/(l**2*tar)-1),1))
s=np.linspace(0,8,17); print(np.round(F(s),3),max(F(np.linspace(0,100,10000))))
