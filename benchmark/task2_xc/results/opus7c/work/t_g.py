from hfeval import *
from diag import parts, M
from exl import get
C = load()
def fx_g(ra, rb, saa, sab, sbb, ta, tb, mode='g'):
    A = (ra, rb, saa, sab, sbb, ta, tb)
    P = parts(*A)
    n = ra+rb; zeta = (ra-rb)/n
    ssa = M._ec_pol(ra, saa); ssb = M._ec_pol(rb, sbb); os_ = P['cPBE'] - ssa - ssb
    za = M._zsat(saa/(8*ra*ta)); zb = M._zsat(sbb/(8*rb*tb))
    e = 0
    for l, ss, zs, zo in (('a', ssa, za, zb), ('b', ssb, zb, za)):
        rem = zs**2*(1-zo**2*(1-zeta**2))*ss
        g = rem/(ss + 0.5*os_ - 1e-30)
        w = P['zs_'+l]**2*g if mode == 'g' else P['w_'+l]
        e = e + (1-w)*P['base_'+l] + w*P['br_'+l]
    return e
def fc(*A):
    P = parts(*A); return P['cPBE'] + P['cgate']
for a, w, cap in [(0, 0.5767436061911404, 1.174), (0.25, 0.43515205200281454, 1.174/0.75), (0.25, 0.43515205200281454, 1.174)]:
    M.OMEGA = w; M.F_LO = cap
    for nm in C: C[nm]['Exl'] = get(nm, w)
    for mode in ['zeta2', 'g']:
        score(species_terms(C, lambda *A: fx_g(*A, mode=mode), fc, a=a), tag=f'{mode} a={a} w={w:.4f} cap={cap:.3f}')
