import numpy as np, pickle
from scipy.optimize import least_squares
allr=[r for r in pickle.load(open('lib.pkl','rb')) if r['k']==0]
names=['H','Li','H2+1','H2+2','H2+3','H2+4.5','H2+7','LiH2.5','LiH4','HeH2','HeH3.5']
sel=[0,2,3,4,9,10]
D=[]
for i in sel:
    r=allr[i]
    n=np.maximum(r['n'],1e-30); n2=2*n; sig2=4*r['sig']
    kf=(3*np.pi**2*n2)**(1/3); s=np.sqrt(sig2)/(2*kf*n2)
    exl=-0.75*(3/np.pi)**(1/3)*n2**(4/3)*0.5
    D.append((r['w'],exl,s,-r['J']))
tar=np.array([d[3] for d in D])
Ex=lambda F:np.array([np.sum(w*exl*F(s)) for w,exl,s,_ in D])
kap,mu=0.804,0.2195149727645171
pbe=lambda s:1+kap-kap/(1+mu*s**2/kap)
print('PBE',np.round(100*(Ex(pbe)/tar-1),1))
def mk(p): a,b,c=p; return lambda s:1+(a*s**2+c*s**3)/(1+b*s**3.5)
for lo in ([0,1e-3,0],):
    best=None
    for p0 in ([.3,.3,.1],[1,1,1],[.5,2,.2],[.3,.1,.0]):
        sol=least_squares(lambda p:Ex(mk(p))/tar-1,p0,bounds=([0,1e-3,0],[50,500,50]))
        if best is None or sol.cost<best.cost: best=sol
    p=best.x; F=mk(p); print(p,np.round(100*(Ex(F)/tar-1),2)); s=np.linspace(0,6,13); print(np.round(F(s),3),max(F(np.linspace(0,100,10000))))
# sonnet3 F1 (CTAIL form)
CT=4/(5*kap/mu)**1.25
F3=lambda s:1+mu*s**2/(1+CT*s**2.5)
print('sonnet3 orb',np.round(100*(Ex(F3)/tar-1),1))
