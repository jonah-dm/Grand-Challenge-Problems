from hfeval import *
from diag import parts, M
from exl import get
from t_zs import fx_ch
C = load()
def fx_tot(*A):
    P = parts(*A); return sum((1-P['w_'+l])*P['base_'+l] + P['w_'+l]*P['br_'+l] for l in 'ab')
def fc(*A):
    P = parts(*A); return P['cPBE'] + P['cgate']
cands = [('S0 K2b', fx_tot, 0, 0.5767436061911404, 1.174),
         ('S1', fx_ch, 0, 0.5767436061911404, 1.174),
         ('S2', fx_ch, 0.25, 0.5767436061911404, 1.174),
         ('S3 R2', fx_ch, 0.25, 0.43515205200281454, 1.174/0.75),
         ('S4 R3', fx_ch, 0.25, 0.48103212833008924, 1.422658968713173)]
for tag, f, a, w, cap in cands:
    M.OMEGA = w; M.F_LO = cap
    for nm in C: C[nm]['Exl'] = get(nm, w)
    score(species_terms(C, f, fc, a=a), tag=tag)
