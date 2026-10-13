from hfeval import *
from diag import parts, M
C = load()
for nm in ['H_atom','N2','N_atom']:
    zw = np.load(os.path.join(D, 'cache_cc-pvtz', nm + '_w.npz'))
    z = C[nm]; A = args_of(z); m = (A[0]+A[1]) > 1e-11; A = [x[m] for x in A]; w = z['w'][m]
    print(nm, 'Ex', z['Ex'].sum(), 'Exl(0.5767)', z['Exl'].sum())
    for i, om in enumerate(zw['ws']):
        M.OMEGA = om
        with np.errstate(all='ignore'):
            P = parts(*A)
        sl = float(np.sum(w*np.nan_to_num(P['base_a']+P['base_b'])))
        print('  w=%.2f ExLR=%.5f ExSR_exact=%.5f  SRsl=%.5f  diff=%.2f kcal' % (om, zw['Exl'][i].sum(), z['Ex'].sum()-zw['Exl'][i].sum(), sl, (sl-(z['Ex'].sum()-zw['Exl'][i].sum()))*H2K))
