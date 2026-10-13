import numpy as np, pickle
from scipy.optimize import least_squares
rows=[r for r in pickle.load(open('lib.pkl','rb')) if r['k']==0]
rows=[r for i,r in enumerate(rows) if i in (0,2,3,4,5,6,9,10)]   # H, H2+ x5, HeH2+ x2 (drop duplicate hydrogenic bare-Li, LiH collapsed)
lams=[1.0,0.75,0.5,0.25]
D=[]
for r in rows:
    n=np.maximum(r['n'],1e-30); n2=2*n; sig2=4*r['sig']
    kf=(3*np.pi**2*n2)**(1/3); s=np.sqrt(sig2)/(2*kf*n2)
    exl=-0.75*(3/np.pi)**(1/3)*n2**(4/3)*0.5
    D.append((r['w'],exl,s,-r['J']))
def Ex(F,lam=1.0):
    return np.array([lam**(4/3)*np.sum(w*exl*F(s*lam**(-1/3))) for w,exl,s,_ in D])
tar=np.array([d[3] for d in D])
kap,mu=0.804,0.2195149727645171
pbe=lambda s:1+kap-kap/(1+mu*s**2/kap)
def report(name,F,lams=lams):
    print(name)
    for lam in lams:
        e=Ex(F,lam)/(lam**2*tar)-1
        print('  lam',lam,np.round(100*e,1))
def mk(p,tail):
    a,b,c=p
    return lambda s:1+a*s**2/(1+b*s**2+c*s**(2+tail))
report('PBE',pbe)
# family: tail s^-2: 1 + a s^2/(1 + b s^4)
def ff(p): a,b=p; return lambda s:1+a*s**2/(1+b*s**4)
def res(p,lams=lams):
    return np.concatenate([Ex(ff(p),l)/(l**2*tar)-1 for l in lams])
best=None
for p0 in ([1,1],[0.3,0.3],[3,10],[0.5,3]):
    sol=least_squares(res,p0,bounds=([0,1e-5],[200,2000]))
    if best is None or sol.cost<best.cost: best=sol
print(best.x,best.cost); report('s^-2 tail 2par',ff(best.x))
print('maxF',max(ff(best.x)(np.linspace(0,60,6000))))
def ff3(p): a,b,c=p; return lambda s:1+(a*s**2+c*s**4)/(1+b*s**6)    # tail s^-2? c/b s^-2
def res3(p): return np.concatenate([Ex(ff3(p),l)/(l**2*tar)-1 for l in lams])
best=None
for p0 in ([1,1,1],[.3,.3,.3],[3,10,10]):
    sol=least_squares(res3,p0,bounds=([0,1e-5,0],[200,2000,2000]))
    if best is None or sol.cost<best.cost: best=sol
print(best.x,best.cost); report('3par',ff3(best.x)); print('maxF',max(ff3(best.x)(np.linspace(0,60,6000))))
print('=========== free F(0)')
def g2(p): A,B=p; return lambda s:A/(1+B*s**2)
def g3(p): A,C,B=p; return lambda s:(A+C*s**2)/(1+B*s**4)
for name,g,p0s,bd in (('A/(1+Bs2)',g2,([1,1],[2,.5],[.5,.1]),([0,0],[50,1e3])),('(A+Cs2)/(1+Bs4)',g3,([1,1,1],[1,.1,.5],[2,.3,.1]),([0,0,0],[50,100,1e3]))):
    best=None
    for p0 in p0s:
        sol=least_squares(lambda p:np.concatenate([Ex(g(p),l)/(l**2*tar)-1 for l in lams]),p0,bounds=bd)
        if best is None or sol.cost<best.cost: best=sol
    print(name,best.x,best.cost); report(name,g(best.x)); s=np.linspace(0,6,13); print(np.round(g(best.x)(s),3),'max',max(g(best.x)(np.linspace(0,60,6000))))
