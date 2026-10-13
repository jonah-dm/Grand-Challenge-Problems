import numpy as np
from hbranch import build, S0SQ
c = 0.18949707
F0, t, FF = build(c, 6.0, npts=4000)
# resample on 161 points uniform in t for embedding
tt = np.linspace(t[0], t[-1], 161)
from scipy.interpolate import interp1d
Ft = interp1d(t, FF, kind="cubic")(tt); tt[0] = 0.0
Ft[-1] = max(Ft[-1], 0.0)
np.save('fh_table.npy', np.array([tt, Ft]))
print("F(s0)=", F0, "t_end=", t[-1], "s_end=", np.sqrt(S0SQ*np.exp(t[-1])), "F_end=", Ft[-1])
print("tt=", repr(np.round(tt, 8).tolist())[:200])
s = np.sqrt(S0SQ*np.exp(tt)); print(np.round(np.c_[s[::16], Ft[::16]], 4))
