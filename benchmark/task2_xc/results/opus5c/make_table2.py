import numpy as np
from scipy.optimize import brentq
from hbranch import shoot, S0SQ
c = 0.2254754
F0 = brentq(lambda f: shoot(f, c, 6.0).y[0,-1], -5, 10, xtol=1e-13); sol = shoot(F0, c, 6.0)
tt = np.linspace(0, 8.0, 241); rr = np.maximum(3*tt/4, 1e-5)
Ft = sol.sol(rr)[0]
open('fh_table_literal.txt','w').write("_FH_T = np.array(" + repr([float('%.10g'%x) for x in tt]) + ")\n_FH_F = np.array(" + repr([float('%.10g'%x) for x in Ft]) + ")\n")
print("F(s0)", F0, "F at s=", [(round(float(np.sqrt(S0SQ*np.exp(t))),2), round(float(f),3)) for t,f in zip(tt[::24], Ft[::24])])
