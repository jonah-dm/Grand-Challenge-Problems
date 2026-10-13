import numpy as np, pickle
from scipy.optimize import brentq
allr=[r for r in pickle.load(open('lib.pkl','rb')) if r['k']==0]
sel=[0,2,3,4,5,6]
D=[]
for i in sel:
    r=allr[i]
    m=r['n']>1e-12
    n=r['n'][m]; n2=2*n; sig2=4*r['sig'][m]
    kf=(3*np.pi**2*n2)**(1/3); s=np.sqrt(sig2)/(2*kf*n2)
    exl=-0.75*(3/np.pi)**(1/3)*n2**(4/3)*0.5
    D.append((r['w'][m],exl,s,-r['J']))
tar=np.array([d[3] for d in D])
def Ex(F,lam=1.0): return np.array([lam**(4/3)*np.sum(w*exl*F(s*lam**(-1/3))) for w,exl,s,_ in D])
for p in (4,3,2.5):
  for b in (0.05,0.1,0.3,1,3,10):
    f=lambda a:Ex(lambda s:1+a*s**2/(1+b*s**p))[0]/tar[0]-1
    try: a=brentq(f,0.01,200)
    except: continue
    F=lambda s:1+a*s**2/(1+b*s**p)
    errs={l:np.round(100*(Ex(F,l)/(l**2*tar)-1),0) for l in (1,.5)}
    print(p,b,round(a,3),'maxF',round(max(F(np.linspace(0,30,3000))),3),'lam1',errs[1],'lam.5 H',errs[.5][0])
print('--- pure C s^-2 (exact amplitude scaling), C from H')
f=lambda C:Ex(lambda s:C/np.maximum(s,1e-8)**2)[0]/tar[0]-1
C=brentq(f,1e-4,100); F=lambda s:C/np.maximum(s,1e-8)**2
print(C,'lam1',np.round(100*(Ex(F)/tar-1),1))
for l in (.75,.5,.25): print(l,np.round(100*(Ex(F,l)/(l**2*tar)-1),1))
