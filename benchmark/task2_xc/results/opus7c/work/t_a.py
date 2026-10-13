from hfeval import *
from diag import parts, M
from exl import get
C = load()
def fx(*A):
    P = parts(*A); return sum((1-P['w_'+l])*P['base_'+l] + P['w_'+l]*P['br_'+l] for l in 'ab')
def fc(*A):
    P = parts(*A); return P['cPBE'] + P['cgate']
for a, w, cap in [(0, 0.5767436061911404, 1.174), (0.25, 0.43515205200281454, 1.174/0.75), (0.25, 0.43515205200281454, 1.174), (0.2, 0.4660213072225416, 1.174/0.8)]:
    M.OMEGA = w; M.F_LO = cap
    for nm in C: C[nm]['Exl'] = get(nm, w)
    score(species_terms(C, fx, fc, a=a), tag=f'a={a} w={w:.4f} cap={cap:.3f}')
