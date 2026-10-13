import numpy as np, importlib.util, os
from scipy.integrate import quad
from scipy.special import erf
from scipy.optimize import brentq
spec = importlib.util.spec_from_file_location('o6', '../../opus6c/submission.py'); M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M)
n0 = 1/np.pi
eLDA = -0.75*(6/np.pi)**(1/3)*n0**(4/3)   # fully polarized LDA at nucleus
def vlr(w): return quad(lambda r: 4*r*np.exp(-2*r)*erf(w*r), 0, 60)[0]
def Ftot(w, a):
    eLR = -0.5*n0*vlr(w); eSR = -0.5*n0*(1-vlr(w))
    Fbr = 5/9*(6*np.pi**2*n0)**(2/3)/w**2*float(M.Qfun(np.array([1/w**2]))[0])  # branch enhancement vs polarized LDA
    return (eLR + a*eSR)/eLDA + (1-a)*Fbr, eLR/eLDA, eSR/eLDA, Fbr
for w in [0.3,0.4,0.5,0.5767,0.7,1.0,2.0]:
    print(w, [round(x,4) for x in Ftot(w,0)], round(Ftot(w,0.25)[0],4))
for a in [0, 0.25]:
    try: print('a',a,'w*',brentq(lambda w: Ftot(w,a)[0]-1.174, 0.2, 5))
    except Exception as e: print(a, e)
