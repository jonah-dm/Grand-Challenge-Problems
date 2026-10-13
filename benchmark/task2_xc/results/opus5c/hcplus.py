import numpy as np
from hbranch import build, S0SQ
def ok(c):
    F0, t, FF = build(c, 6.0)
    return FF.min() > 0, FF[-1], np.sqrt(S0SQ*np.exp(t[-1]))
lo, hi = 0.18, 0.20
for _ in range(14):
    mid = 0.5*(lo+hi); good, fm, sm = ok(mid)
    if good: hi = mid
    else: lo = mid
print("c+ =", hi, ok(hi))
for c in [hi, 0.19, 0.195, 0.2]:
    print(c, ok(c))
