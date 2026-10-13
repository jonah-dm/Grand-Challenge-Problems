import numpy as np, pickle
from scipy.optimize import brentq, least_squares
l2=pickle.load(open('lib2.pkl','rb')); l3=pickle.load(open('lib3.pkl','rb'))
items=[('H',l2[0]),('H2+R1',l2[2]),('R2',l2[3]),('R3',l2[4]),('R4.5',l2[5]),('R7',l2[6]),('HeH2',l2[9]),('HeH3.5',l2[10])]+[(r['frame'],r) for r in l3]
D=[]
for name,r in items:
    m=r['n']>1e-12
    n=r['n'][m]; n2=2*n; sig2=4*r['sig'][m]
    kf=(3*np.pi**2*n2)**(1/3); s=np.sqrt(sig2)/(2*kf*n2)
    D.append((name,r['w'][m],-0.75*(3/np.pi)**(1/3)*n2**(4/3)*0.5,s,-r['J']))
names=[d[0] for d in D]; tar=np.array([d[4] for d in D])
Ex=lambda F:np.array([np.sum(w*exl*F(s)) for _,w,exl,s,_ in D])
kap,mu=0.804,0.2195149727645171
def rep(lab,F):
    e=Ex(F)/tar-1; print('%-26s'%lab,' '.join('%s %+.1f'%(n,100*x) for n,x in zip(names,e)))
rep('PBE',lambda s:1+kap-kap/(1+mu*s**2/kap))
f=lambda b:Ex(lambda s:1+mu*s**2/(1+b*s**4))[0]/tar[0]-1
bH=brentq(f,1e-3,5); print('bH',bH); rep('H-exact (a=mu)',lambda s:1+mu*s**2/(1+bH*s**4))
# H+He joint, two params a,b (form a s^2/(1+b s^4)) ; indices H=0, He=8
def joint(idx):
    def res(p): 
        F=lambda s:1+p[0]*s**2/(1+p[1]*s**4); return (Ex(F)/tar-1)[idx]
    best=None
    for p0 in ([.3,.05],[.5,.2],[1,1],[.2,.01]):
        sol=least_squares(res,p0,bounds=([0,1e-4],[50,100]))
        if best is None or sol.cost<best.cost: best=sol
    return best.x
for lab,idx in (('H+He exact',[0,8]),('H,He,Ne,Ar,Be,Li 1s',[0,8,9,10,11,12]),('H,He,Ne,Ar + H2+R<=2,HeH2+',[0,1,2,8,9,10,6,7])):
    p=joint(idx); F=lambda s,p=p:1+p[0]*s**2/(1+p[1]*s**4); print(lab,np.round(p,4),'maxF',round(max(F(np.linspace(0,30,3000))),3)); rep(lab,F)
