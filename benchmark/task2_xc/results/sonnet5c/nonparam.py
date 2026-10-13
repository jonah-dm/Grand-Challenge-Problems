"""Non-parametric no-go: Ex[F] = sum_k F_k * A_{ik} is LINEAR in the values F_k of F on a log-s grid
(piecewise-linear F). Solve the minimax LP: min t s.t. |A F / (lam^2 J) - 1| <= t for the one-orbital library at
several amplitudes lam. This is the best any F(s) can do, parametrisation-free."""
import numpy as np, pickle
from scipy.optimize import linprog
allr=pickle.load(open('lib2.pkl','rb'))
sel=[0,2,3,4,5,6]; names=['H','H2+R1','R2','R3','R4.5','R7']
sg=np.exp(np.linspace(np.log(1e-2),np.log(60),70)); K=len(sg)
def hat(s):
    # piecewise-linear interpolation weights on log grid
    x=np.log(np.clip(s,sg[0],sg[-1])); idx=np.clip(np.searchsorted(np.log(sg),x)-1,0,K-2)
    t=(x-np.log(sg)[idx])/(np.log(sg)[idx+1]-np.log(sg)[idx])
    return idx,t
def build(lams,sel_):
    rows=[];tar=[]
    for lam in lams:
        for i in sel_:
            r=allr[i]; m=r['n']>1e-12
            n=r['n'][m]; n2=2*n; sig2=4*r['sig'][m]
            kf=(3*np.pi**2*n2)**(1/3); s=np.sqrt(sig2)/(2*kf*n2)*lam**(-1/3)
            exl=-0.75*(3/np.pi)**(1/3)*n2**(4/3)*0.5*lam**(4/3)
            idx,t=hat(s); a=np.zeros(K)
            np.add.at(a,idx,r['w'][m]*exl*(1-t)); np.add.at(a,idx+1,r['w'][m]*exl*t)
            rows.append(a/(lam**2*(-r['J']))*(-1)*(-1)); tar.append(1.0)
    return np.array(rows)
def solve(lams,sel_,Fmax=1.804,Fmin=0.0):
    A=build(lams,sel_); n=len(A)
    # variables F (K), t ; constraints  A F -1 <= t,  -(A F -1) <= t
    c=np.zeros(K+1); c[-1]=1
    Aub=np.vstack([np.hstack([A,-np.ones((n,1))]),np.hstack([-A,-np.ones((n,1))])]); bub=np.concatenate([np.ones(n),-np.ones(n)])
    bounds=[(Fmin,Fmax)]*K+[(0,None)]
    r=linprog(c,A_ub=Aub,b_ub=bub,bounds=bounds,method='highs')
    return r.x[-1],r.x[:K],A
for lams,label in (([1.0],'lam=1 only (shape family H,H2+ R<=7)'),([1.0,0.75],'lam in {1,.75}'),([1.0,0.75,0.5],'lam in {1,.75,.5}'),([1,.75,.5,.25],'lam in {1,.75,.5,.25}')):
    t,F,A=solve(lams,sel); print('%-40s minimax max-rel-error = %.1f%%'%(label,100*t))
for lams,label in (([1.0],'lam=1, bonded only (H, R<=3)'),):
    t,F,A=solve(lams,[0,2,3,4]); print('%-40s minimax = %.2f%%'%(label,100*t))
t,F,A=solve([1.0,0.75,0.5],[0,2,3,4]); print('bonded shapes + lam {1,.75,.5}: %.1f%%'%(100*t))
print('--- smooth F exact at lam=1 for H, H2+ R=1..7 (|err|<=tol)')
def smooth(sel_,tol,lams=[1.0],Fmax=1.804):
    A=build(lams,sel_); n=len(A)
    # vars: F(K), u(K-2) >= |second diff|
    nv=K+K-2
    c=np.zeros(nv); c[K:]=1
    rowsU=[];bU=[]
    for i in range(n):
        r=np.zeros(nv); r[:K]=A[i]; rowsU.append(r); bU.append(1+tol)
        r=np.zeros(nv); r[:K]=-A[i]; rowsU.append(r); bU.append(-(1-tol))
    for k in range(1,K-1):
        for sgn in (1,-1):
            r=np.zeros(nv); r[k-1]=sgn;r[k]=-2*sgn;r[k+1]=sgn; r[K+k-1]=-1; rowsU.append(r); bU.append(0)
    res=linprog(c,A_ub=np.array(rowsU),b_ub=np.array(bU),bounds=[(0,Fmax)]*K+[(0,None)]*(K-2),method='highs')
    return res
for tol in (0.01,0.03):
    res=smooth(sel,tol)
    if res.status==0:
        F=res.x[:K]; print('tol',tol,'ok; total variation',round(res.fun,3)); print(np.round(np.interp(np.array([0.05,0.1,0.2,0.3,0.5,0.7,1,1.5,2,3,5,8,15,30]),sg,F),3))
    else: print('tol',tol,'infeasible')
