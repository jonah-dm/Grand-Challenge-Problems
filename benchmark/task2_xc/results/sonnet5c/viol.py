import numpy as np, pickle
from scipy.optimize import brentq
allr=[r for r in pickle.load(open('lib.pkl','rb')) if r['k']==0]
sel=[0,2,3,4,5,6]; names=['H','H2+ R=1','R=2','R=3','R=4.5','R=7']
D=[]
for i in sel:
    r=allr[i]; m=r['n']>1e-12
    n=r['n'][m]; n2=2*n; sig2=4*r['sig'][m]
    kf=(3*np.pi**2*n2)**(1/3); s=np.sqrt(sig2)/(2*kf*n2)
    D.append((r['w'][m],-0.75*(3/np.pi)**(1/3)*n2**(4/3)*0.5,s,-r['J']))
tar=np.array([d[3] for d in D])
def Ex(F,lam=1.0): return np.array([lam**(4/3)*np.sum(w*exl*F(s*lam**(-1/3))) for w,exl,s,_ in D])
kap,mu=0.804,0.2195149727645171
CT=4/(5*kap/mu)**1.25
fs={'LDA':lambda s:1+0*s,'PBE':lambda s:1+kap-kap/(1+mu*s**2/kap),'orb s^-1/2 (kappa-capped)':lambda s:1+mu*s**2/(1+CT*s**2.5)}
f=lambda b:Ex(lambda s:1+mu*s**2/(1+b*s**4))[0]/tar[0]-1
b=brentq(f,1e-4,10); print('b_H(s^-2 tail)=',b,'maxF',1+mu/(2*np.sqrt(b)))
fs['orb s^-2 tail, H-exact']=lambda s:1+mu*s**2/(1+b*s**4)
print('%-28s'%'', *['%8s'%n for n in names])
for k,F in fs.items():
    for lam in (1,.75,.5):
        print('%-22s l=%.2f'%(k,lam),*['%+7.1f%%'%(100*x) for x in (Ex(F,lam)/(lam**2*tar)-1)])
