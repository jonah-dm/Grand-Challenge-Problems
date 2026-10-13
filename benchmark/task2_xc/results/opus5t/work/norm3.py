import numpy as np
from scipy.optimize import brentq
exec(open('norm.py').read().split("c=0.3125")[0])
MU=0.2195149727645171
Fq=lambda s,c: (1+MU*s*s)/(1+MU*s**4/c)
c=brentq(lambda c: Ex(lambda s: Fq(s,c))+0.3125, 5, 200); print('c_q2',repr(c))
ss=np.linspace(0,10,100001); print('max F',Fq(ss,c).max(),'at s',ss[np.argmax(Fq(ss,c))])
e1=Ex(lambda s:Fq(s,c)); print('f(1/2)',Ex(lambda s:Fq(s,c),0.5)/e1,'f(1/4)',Ex(lambda s:Fq(s,c),0.25)/e1)
