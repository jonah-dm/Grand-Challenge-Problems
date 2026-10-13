import numpy as np, sys
sys.path.insert(0, '..'); sys.path.insert(0, '.')
import submission as S
from scipy.optimize import brentq
r = np.exp(np.linspace(np.log(1e-7), np.log(60), 20000)); wr = 4*np.pi*r**2*np.gradient(r)
n = np.exp(-2*r)/np.pi; g = (2*n)**2
def ex(cq):
    S.CQ = cq
    return np.sum(S._ex_ch(n, g, S.F_q)*wr)
cq = brentq(lambda c: ex(c) + 5/16, 0.5, 10)
print('c_q (H norm, exact density) =', cq)
S.CQ = cq
# full functional on exact H, mode C: z=1 exactly (tau=tau_W)
z = np.zeros_like(n); tw = g/(8*n)
for m in 'ABC':
    e = S.edens(n, z, g, z, z, tw, z, mode=m)
    print(m, 'semilocal E_xc[H] =', np.sum(e*wr).real, ' + HF part', -0.25*5/16, ' total', np.sum(e*wr).real - 0.25*5/16)
# fractional H, C1 check and fractional spin (constancy) for each mode, incl. HF part
for m in 'ABC':
    out = []
    for fa, fb in [(0.5, 0), (0.5, 0.5), (0.25, 0.25), (0.75, 0.25)]:
        na, nb = fa*n, fb*n
        e = S.edens(na, nb, fa*fa*g, fa*fb*g, fb*fb*g, fa*tw, fb*tw, mode=m)
        exc = np.sum(e*wr).real - 0.25*(fa*fa+fb*fb)*5/16
        out.append(f"({fa},{fb}) Exc/J={exc/(5/16):+.3f} [exact {-(fa+fb)**2:+.3f}]")
    print(m, ' '.join(out))
b1 = brentq(lambda b: np.sum(S._ex_ch(n, g, lambda s2: 1 + S.MU*s2/(1 + b*s2*s2))*wr) + 5/16, 1e-4, 10)
print('B1 (H norm of F_1) =', b1)
