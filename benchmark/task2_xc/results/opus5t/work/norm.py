# Solve the H-atom norm for F_q(s) = F_LO/(1+F_LO s^2/c_q) (spin-scaled s), and check T2 constants.
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
cx=0.75*(3/np.pi)**(1/3); FLO=1.804
def Ex(F, lam=1.0, shape='H'):
    def integrand(r):
        if shape=='H':
            n=lam*np.exp(-2*r)/np.pi; g=2*n
        else:  # gaussian, a=1
            a=1.0; n=lam*(2*a/np.pi)**1.5*np.exp(-2*a*r*r); g=4*a*r*n
        n2=2*n; kf=(3*np.pi**2*n2)**(1/3); s=2*g/(2*kf*n2)
        return 0.5*(-cx*n2**(4/3))*F(s)*4*np.pi*r*r
    return quad(integrand,0,60,limit=400,points=[0.5,1,2,5])[0]
c=0.3125/(-Ex(lambda s: 1/s**2)); print('pure c',c)
print('gauss ratio', Ex(lambda s:c/s**2,shape='G')/(-1/np.sqrt(np.pi)))
f=lambda cq: Ex(lambda s: FLO/(1+FLO*s*s/cq))+0.3125
cq=brentq(f,0.05,50); print('c_q',repr(cq))
for lam in (0.25,0.5,0.75):
    print(lam,'f_q',Ex(lambda s: FLO/(1+FLO*s*s/cq),lam)/-0.3125,'lda',lam**(4/3),'exact',lam**2)
KAPPA,MU=0.804,0.2195149727645171; CT=4.0/(5.0*KAPPA/MU)**1.25
Fs3=lambda s: 1+MU*s*s/(1+CT*s**2.5)
Fpbe=lambda s: 1+KAPPA-KAPPA/(1+MU*s*s/KAPPA)
for nm,F in (('sonnet3-iso',Fs3),('PBE',Fpbe),('LDA',lambda s:1+0*s),('Fq',lambda s: FLO/(1+FLO*s*s/cq))):
    e1=Ex(F); print(nm,'Ex[H]',e1,'f(1/2)',Ex(F,0.5)/e1, 'f(1/4)',Ex(F,0.25)/e1)
def ExG(F):
    return quad(lambda r: 0.5*(-cx*(2*(2/np.pi)**1.5*np.exp(-2*r*r))**(4/3))*F(4*r/(2*(3*np.pi**2*2*(2/np.pi)**1.5*np.exp(-2*r*r))**(1/3)))*4*np.pi*r*r,1e-6,12,limit=400)[0]
print('gauss: pure c ratio',ExG(lambda s:c/s**2)/(-1/np.sqrt(np.pi)),' Fq ratio',ExG(lambda s: FLO/(1+FLO*s*s/cq))/(-1/np.sqrt(np.pi)),' sonnet3',ExG(Fs3)/(-1/np.sqrt(np.pi)))
