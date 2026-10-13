from hfeval import *
from diag import parts, M
from pyscf.dft import libxc
C = load()
for nm in C:
    zw = np.load(os.path.join(D, 'cache_cc-pvtz', nm + '_w.npz')); C[nm]['ws'] = zw['ws']; C[nm]['ExlW'] = zw['Exl']; C[nm]['Exl0'] = C[nm]['Exl']
def lx(code, A, om=None):
    ra, rb, saa, sab, sbb, ta, tb = A
    # build pyscf spin rho (2,4,N): need gradients; we only have sigmas -> fake gradient along x
    N = len(ra)
    rho = np.zeros((2, 6, N))
    rho[0, 0] = ra; rho[1, 0] = rb; rho[0, 1] = np.sqrt(saa); rho[1, 1] = np.sqrt(sbb)
    rho[0, 5] = ta; rho[1, 5] = tb
    exc = libxc.eval_xc(code, rho[:, [0,1,2,3,5]] if 'MGGA' in code or 'SCAN' in code else rho[:, :4], spin=1, deriv=0, omega=om)[0]
    return exc*(ra+rb)
def branchmix(A, base):
    P = parts(*A)
    return sum(P['w_'+l]*(P['br_'+l]) for l in 'ab'), sum(P['w_'+l] for l in 'ab')
import itertools
def fc_o4(*A):
    P = parts(*A); return P['cPBE'] + P['cgate']
def fc_pbe(*A): return lx('GGA_C_PBE', A)
for wi, om in [(None, 0.5767436061911404), (2, 0.4), (1, 0.3)]:
    M.OMEGA = om
    for nm in C: C[nm]['Exl'] = C[nm]['Exl0'] if wi is None else C[nm]['ExlW'][wi]
    def fx_hjs(*A):
        # per-channel: HJS SR exchange of each spin channel (libxc handles spin scaling)
        return lx('GGA_X_HJS_PBE', A, om)
    def fx_hjs_br(*A):
        ra, rb, saa, sab, sbb, ta, tb = A
        z0 = 1e-14*np.ones_like(ra)
        ea = lx('GGA_X_HJS_PBE', [ra, z0, saa, 0*ra, 0*ra, ta, 1e-20+0*ra], om)
        eb = lx('GGA_X_HJS_PBE', [z0, rb, 0*ra, 0*ra, sbb, 1e-20+0*ra, tb], om)
        P = parts(*A)
        return (1-P['w_a'])*ea + P['w_a']*P['br_a'] + (1-P['w_b'])*eb + P['w_b']*P['br_b']
    for a in [0.0, 0.25]:
        score(species_terms(C, fx_hjs, fc_pbe, a=a), tag=f'w={om:.3f} a={a} HJS/PBEc (LC-wPBE-like)')
        score(species_terms(C, fx_hjs_br, fc_o4, a=a), tag=f'w={om:.3f} a={a} HJS+branch/o4')
