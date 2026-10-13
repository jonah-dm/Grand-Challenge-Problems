import importlib.util, os, numpy as np
from hfeval import *
spec = importlib.util.spec_from_file_location('o6', os.path.join(D, '../../opus6c/submission.py'))
M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M)
from diag import parts
C = load()
def mk(xkind, ckind):
    def fx(*A):
        P = parts(*A)
        if xkind == 'o6': return sum((1-P['w_'+l])*P['base_'+l] + P['w_'+l]*P['br_'+l] for l in 'ab')
        if xkind == 'base': return P['base_a'] + P['base_b']
    def fc(*A):
        P = parts(*A)
        return P['cPBE'] + (P['cgate'] if ckind == 'o4' else 0)
    return fx, fc
for xk in ['o6', 'base']:
    for ck in ['o4', 'pbe']:
        fx, fc = mk(xk, ck)
        score(species_terms(C, fx, fc), tag=f'{xk}/{ck}')
