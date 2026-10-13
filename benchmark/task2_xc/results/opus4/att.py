import numpy as np
from scipy.special import erf
def Fd(a):
    return 1 - 8/3*a*(np.sqrt(np.pi)*erf(1/(2*a)) + (2*a-4*a**3)*np.exp(-1/(4*a**2)) - 3*a + 4*a**3)
def Fs(a):
    x=1/a**2
    return x/36 - x**2/960 + x**3/26880 - x**4/829440
import mpmath as mp
mp.mp.dps=50
def Fm(a):
    a=mp.mpf(a); return 1-mp.mpf(8)/3*a*(mp.sqrt(mp.pi)*mp.erf(1/(2*a))+(2*a-4*a**3)*mp.exp(-1/(4*a**2))-3*a+4*a**3)
for a in [0.5,1,2,3,5,10,50]:
    print(a, float(Fm(a)), Fd(a), Fs(a))
