import numpy as np
from scipy.special import erfcx
from scipy.integrate import quad
def J4(b):
    E = erfcx(np.sqrt(b))
    J1 = np.pi/2*E
    J2 = np.pi/4*(1-2*b)*E + 0.5*np.sqrt(np.pi*b)
    J3 = ((3-2*b)*J2 + 2*b*J1)/4
    return ((5-2*b)*J3 + 2*b*J2)/6
def R(b): return 1-32/(5*np.pi)*J4(b)
for b in [1e-3,0.1,1,6.25,11.1,50,300,1000]:
    q=quad(lambda u: np.exp(-b*u*u)/(1+u*u)**4,0,np.inf,limit=200)[0]
    print(b, R(b), 1-32/(5*np.pi)*q)
from scipy.optimize import brentq
c=5/9*(6*np.pi)**(2/3)
w=brentq(lambda w: c*R(1/w**2)-1.804,0.1,2); print('c',c,'omega*',repr(w))
print(R(1/0.3**2),R(1/0.4**2))
z=4.0+1e-20j; print((R(z).imag/1e-20), (R(4+1e-6)-R(4-1e-6))/2e-6)
