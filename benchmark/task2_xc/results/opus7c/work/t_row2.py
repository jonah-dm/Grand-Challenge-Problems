import glob, os, numpy as np
from hfeval import args_of, H2K
from diag import parts, M
from t_zs import fx_ch, fc
def fx_tot(*A):
    P = parts(*A); return sum((1-P['w_'+l])*P['base_'+l] + P['w_'+l]*P['br_'+l] for l in 'ab')
C = {os.path.basename(f)[:-4]: dict(np.load(f)) for f in glob.glob('cache_row2/*.npz')}
items = {'AE:HCl': {'H_atom': 1, 'Cl_atom': 1, 'HCl': -1}, 'AE:H2S': {'H_atom': 2, 'S_atom': 1, 'H2S': -1},
         'AE:PH3': {'H_atom': 3, 'P_atom': 1, 'PH3': -1}, 'AE:SiH4': {'H_atom': 4, 'Si_atom': 1, 'SiH4': -1},
         'AE:Cl2': {'Cl_atom': 2, 'Cl2': -1}, 'IP:Cl': {'Cl_cation': 1, 'Cl_atom': -1}, 'IP:S': {'S_cation': 1, 'S_atom': -1}}
for tag, f, wi, cap in [('S0', fx_tot, 0, 1.174), ('S1', fx_ch, 0, 1.174), ('S5', fx_ch, 1, 1.299259029906918)]:
    M.OMEGA = [0.5767436061911404, 0.5256070335369912][wi]; M.F_LO = cap
    d = {}
    for nm, z in C.items():
        A = args_of(z); m = (A[0]+A[1]) > 1e-11; A = [x[m] for x in A]
        with np.errstate(all='ignore'):
            exsl = np.sum(z['w'][m]*np.nan_to_num(f(*A))); ec = np.sum(z['w'][m]*np.nan_to_num(fc(*A)))
        d[nm] = (z['ExlW'][wi].sum() + exsl - z['Ex'].sum(), ec - float(z['ec']))
    out = []
    for k, ex in items.items():
        if all(s in d for s in ex):
            x = sum(c*d[s][0] for s, c in ex.items())*H2K; c = sum(c*d[s][1] for s, c in ex.items())*H2K
            out.append(f'{k}:{x+c:+.1f}(x{x:+.1f},c{c:+.1f})')
    print(tag, ' '.join(out))
