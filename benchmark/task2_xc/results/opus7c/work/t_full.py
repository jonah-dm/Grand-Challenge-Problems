from hfeval import *
from diag import parts, M
from t_hjs import lx, C, fc_o4, fc_pbe
for nm in C: C[nm]['Exl'] = 0*C[nm]['Exl0']
def fx_s3(*A):
    M.OMEGA = 1e-8
    P = parts(*A); return P['base_a'] + P['base_b']
for code in ['GGA_X_PBE', 'GGA_X_B88', 'MGGA_X_R2SCAN', 'LDA_X']:
    score(species_terms(C, lambda *A: lx(code, A), fc_pbe), tag=code+'/PBEc')
score(species_terms(C, fx_s3, fc_o4), tag='sonnet3 full-range / o4 (=opus4)')
for nm in C: C[nm]['Exl'] = C[nm]['Exl0']
score(species_terms(C, fx_s3, fc_o4, a=1.0), tag='HF x / o4c')
