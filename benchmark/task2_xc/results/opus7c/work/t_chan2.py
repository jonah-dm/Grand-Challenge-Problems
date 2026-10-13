from hfeval import *
from diag import M
C = load()
def chan(r, s, t, mode):
    n2 = 2*r; kf = (3*np.pi**2*n2)**(1/3); s2 = 4*s/(4*kf**2*n2**2)
    z = {'pbe': 0*r, 'orb': 1+0*r, 'chan': M._zsat(s/(8*r*t))}[mode]
    return 0.5*M._ex_base_sr(n2, s2, z)
print('species   chan  exactSR   err(pbe)  err(orb)  err(chan)  [kcal]')
for nm in ['O_atom','F_atom','F_cation','HF','HF_cation','B_atom','B_cation','Li_atom']:
    z = C[nm]; A = args_of(z); m = (A[0]+A[1]) > 1e-11; A = [x[m] for x in A]; w = z['w'][m]
    for i, (r, s, t) in enumerate(((A[0], A[2], A[5]), (A[1], A[4], A[6]))):
        ex = z['Ex'][i]-z['Exl'][i]
        with np.errstate(all='ignore'):
            e = [float(np.sum(w*np.nan_to_num(chan(r, s, t, md)))) for md in ['pbe','orb','chan']]
        print(f"{nm:10s} {'ab'[i]} {ex:9.4f} " + ' '.join(f'{(x-ex)*H2K:8.1f}' for x in e))
