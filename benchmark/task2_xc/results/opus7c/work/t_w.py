from hfeval import *
from diag import parts, M
import diag
C = load()
for nm in C:
    zw = np.load(os.path.join(D, 'cache_cc-pvtz', nm + '_w.npz')); C[nm]['ws'] = zw['ws']; C[nm]['ExlW'] = zw['Exl']
def fx(*A, br=True):
    P = parts(*A)
    if br: return sum((1-P['w_'+l])*P['base_'+l] + P['w_'+l]*P['br_'+l] for l in 'ab')
    return P['base_a']+P['base_b']
def fc(*A):
    P = parts(*A); return P['cPBE'] + P['cgate']
import sys
for a in [0.0, 0.25]:
  for i, w in enumerate(C['H_atom']['ws']):
    M.OMEGA = w
    for nm in C: C[nm]['Exl'] = C[nm]['ExlW'][i]
    for br in [True]:
        score(species_terms(C, lambda *A: fx(*A, br=br), fc, a=a), tag=f'a={a} w={w} br={br}')
