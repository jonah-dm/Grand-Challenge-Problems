import numpy as np, importlib.util
from scipy.optimize import brentq
from scipy.interpolate import CubicSpline
from hbranch import shoot, S0SQ
import atomnorm as A   # loads HF densities (prints suppressed below)
FMAX = 1.804
def Fcap(c, rmax=7.0):
    F0 = brentq(lambda f: shoot(f, c, rmax).y[0,-1], -5, 10, xtol=1e-13); sol = shoot(F0, c, rmax)
    rr = np.linspace(1e-5, 4.5, 3000); FF = sol.sol(rr)[0]
    k = np.argmax(FF >= FMAX) if np.any(FF >= FMAX) else len(rr)-1
    t = 4*rr[:k+1]/3; FF = FF[:k+1]; FF[-1] = FMAX
    cs = CubicSpline(t, FF, bc_type=((1,0.0),'natural'))
    f = lambda s2: np.where(np.log(s2/S0SQ)<=0, FF[0], np.where(np.log(s2/S0SQ)>=t[-1], FMAX, cs(np.clip(np.log(s2/S0SQ),0,t[-1]))))
    return f, t, FF, F0
c = brentq(lambda c: A.ex("Ne", Fcap(c)[0]) - A.data["Ne"][0], 0.2, 0.26, xtol=1e-7)
f, t, FF, F0 = Fcap(c)
print("c(Ne, capped) =", c, "F(s0)=", F0, "cap at s=", np.sqrt(S0SQ*np.exp(t[-1])), {at: round(A.ex(at,f)-A.data[at][0],5) for at in A.data})
np.save('fh_cap.npy', np.array([t, FF])); open('c_value.txt','w').write(repr(c))
