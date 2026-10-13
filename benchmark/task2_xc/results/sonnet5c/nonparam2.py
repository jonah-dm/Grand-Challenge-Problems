import numpy as np, pickle
exec(open('nonparam.py').read().split("for lams,label in")[0])
from scipy.optimize import brentq
mu=0.2195149727645171
sel=[0,2,3,4,5,6,9,10]
A=build([1.0],sel)
# reference parametric H-exact form
bH=0.0377296
Fref=1+mu*sg**2/(1+bH*sg**4)
w=np.gradient(np.log(sg))  # log-grid weights
def solve_dev(tol,Fmax=1.804):
    n=len(A); nv=2*K
    c=np.concatenate([np.zeros(K),w])
    U=[];b=[]
    for i in range(n):
        r=np.zeros(nv); r[:K]=A[i]; U.append(r); b.append(1+tol)
        r=np.zeros(nv); r[:K]=-A[i]; U.append(r); b.append(-(1-tol))
    for k in range(K):
        r=np.zeros(nv); r[k]=1; r[K+k]=-1; U.append(r); b.append(Fref[k])
        r=np.zeros(nv); r[k]=-1; r[K+k]=-1; U.append(r); b.append(-Fref[k])
    res=linprog(c,A_ub=np.array(U),b_ub=np.array(b),bounds=[(0,Fmax)]*K+[(0,None)]*K,method='highs')
    return res
print('F_ref errors',np.round(100*(A@Fref-1),1))
S=np.array([0.1,0.3,0.5,0.7,1,1.5,2,3,5,8,15])
print('s',S); print('Fref',np.round(np.interp(S,sg,Fref),3))
for tol in (0.005,0.01,0.02):
    res=solve_dev(tol)
    F=res.x[:K]; print(tol,'dev',round(res.fun,3),np.round(np.interp(S,sg,F),3),'errs',np.round(100*(A@F-1),1))
np.save('F_lp_tol1.npy',solve_dev(0.01).x[:K]); np.save('sg.npy',sg)
