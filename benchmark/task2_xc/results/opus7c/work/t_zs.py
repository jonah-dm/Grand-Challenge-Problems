from hfeval import *
from diag import parts, M
C = load()
def fx_ch(*A, zmode='chan', branch=True):
    ra, rb, saa, sab, sbb, ta, tb = A
    n = ra+rb; zeta = (ra-rb)/n
    e = 0
    for r, s, t in ((ra, saa, ta), (rb, sbb, tb)):
        n2 = 2*r; kf = (3*np.pi**2*n2)**(1/3); s2 = 4*s/(4*kf**2*n2**2)
        zc = M._zsat(s/(8*r*t))
        eb = 0.5*M._ex_base_sr(n2, s2, zc)
        if branch:
            tu_ = 0.3*(6*np.pi**2)**(2/3)*r**(5/3)*M.ZREG
            zs_ = M._zsat((s/(8*r)+tu_)/(t+tu_)); w = zeta**2*zs_**2
            eb = (1-w)*eb + w*0.5*M._ex_branch(r, s)
        e = e + eb
    return e
def fc(*A):
    P = parts(*A); return P['cPBE'] + P['cgate']
def fcp(*A):
    P = parts(*A); return P['cPBE']
score(species_terms(C, fx_ch, fc), tag='chan-z + branch / o4')
score(species_terms(C, lambda *A: fx_ch(*A, branch=False), fc), tag='chan-z nobranch / o4')
score(species_terms(C, fx_ch, fcp), tag='chan-z + branch / pbe')
