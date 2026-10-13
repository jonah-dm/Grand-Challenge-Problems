import numpy as np, pickle
from scipy.optimize import least_squares
rows=pickle.load(open('lib.pkl','rb'))
data=[]
for r in rows:
    n=np.maximum(r['n'],1e-30); n2=2*n; sig2=4*r['sig']
    kf=(3*np.pi**2*n2)**(1/3); s=np.sqrt(sig2)/(2*kf*n2)
    exl=-0.75*(3/np.pi)**(1/3)*n2**(4/3)*0.5
    data.append((r['w'],exl,s,-r['J']))
def Ex(F):
    return np.array([np.sum(w*exl*F(s)) for w,exl,s,_ in data])
tar=np.array([d[3] for d in data])
kap,mu=0.804,0.2195149727645171
pbe=lambda s:1+kap-kap/(1+mu*s**2/kap)
def rep(name,F):
    e=Ex(F)/tar-1
    print(name,'mean|err|%.2f%% max %.2f%% mean %.2f%%'%(100*np.mean(abs(e)),100*np.max(abs(e)),100*np.mean(e)))
    return e
rep('LDA',lambda s:1+0*s); e=rep('PBE',pbe)
print(np.round(100*e,1))
# orbital s-range
print('s max', [round(float(np.max(d[2][d[0]*d[1]<-1e-7])),1) for d in data][:6])
# form: F=1+ (a s^2+c s^3)/(1+b s^3.5)
def form(p):
    a,b,c=p
    return lambda s:1+(a*s**2+c*s**3)/(1+b*s**3.5)
def res(p):
    return Ex(form(p))/tar-1
best=None
for p0 in ([0.3,0.1,0.1],[1,1,1],[0.5,0.3,0.5],[0.2,0.5,0.0]):
    sol=least_squares(res,p0,bounds=([0,1e-4,0],[20,50,50]))
    if best is None or sol.cost<best.cost: best=sol
p=best.x; print('fit',p); e=rep('F1',form(p)); print(np.round(100*e,1))
s=np.linspace(0,5,26); print(np.round(form(p)(s),3)); print('max F',max(form(p)(np.linspace(0,50,5000))))
