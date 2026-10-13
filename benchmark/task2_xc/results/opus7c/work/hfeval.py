"""Non-SCF evaluation on cached UHF densities: E[n_HF] = E_HF - Ex_HF + Ex_LR + a Ex_SR + (1-a) Ex_sl + Ec."""
import sys, os, json, glob
import numpy as np
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../kit')))
import dev_set as S
H2K = 627.5094740631
REF = json.load(open(os.path.join(os.path.dirname(__file__), '../../../kit/dev_reference.json')))
D = os.path.dirname(os.path.abspath(__file__))

def load(basis='cc-pvtz'):
    C = {}
    for f in glob.glob(os.path.join(D, 'cache_' + basis, '*.npz')):
        if f.endswith('_w.npz'): continue
        z = np.load(f); C[os.path.basename(f)[:-4]] = {k: z[k] for k in z.files}
    return C

def args_of(z):
    R = z['rho']
    ra = np.maximum(R[0][0], 1e-14); rb = np.maximum(R[1][0], 1e-14)
    ga, gb = R[0][1:4], R[1][1:4]
    return [ra, rb, (ga**2).sum(0), (ga*gb).sum(0), (gb**2).sum(0), np.maximum(R[0][4], 1e-20), np.maximum(R[1][4], 1e-20)]

def species_terms(C, fx, fc, a=0.0):
    """fx(args)->SR semilocal exchange energy density (per volume), fc(args)->corr energy density"""
    out = {}
    for nm, z in C.items():
        A = args_of(z); m = (A[0]+A[1]) > 1e-11
        A = [x[m] for x in A]; w = z['w'][m]
        with np.errstate(all='ignore'):
            exsl = float(np.sum(w*np.nan_to_num(fx(*A)))); ec = float(np.sum(w*np.nan_to_num(fc(*A))))
        Ex, Exl = z['Ex'].sum(), z['Exl'].sum()
        dx = Exl + a*(Ex-Exl) + (1-a)*exsl - Ex
        out[nm] = dict(E=REF['species'][nm]['hf_qz'] + dx + ec, dx=dx, ec=ec,
                       ecref=REF['species'][nm]['E_tot']-REF['species'][nm]['hf_qz'] if nm in REF['species'] else 0)
    return out

def score(out, verbose=True, tag=''):
    errs = {}; cls = {}
    for k, (c, expr) in S.ITEMS.items():
        if not all(s in out for s in expr): continue
        val = sum(co*out[s]['E'] for s, co in expr.items())*H2K
        # use reference HF energies at QZ for the HF part? keep own: compare against ref items
        refv = REF['items'][k]['kcal']
        hfpart = sum(co*(REF['species'][s]['hf_qz']) for s, co in expr.items())*H2K
        x = sum(co*out[s]['dx'] for s, co in expr.items())*H2K
        cerr = sum(co*(out[s]['ec']-out[s]['ecref']) for s, co in expr.items())*H2K
        errs[k] = (val-refv, x, cerr)
        cls.setdefault(c, []).append(abs(val-refv))
    summ = {c: np.mean(v) for c, v in cls.items()}
    meanref = {}
    for k, (c, expr) in S.ITEMS.items(): meanref.setdefault(c, []).append(abs(REF['items'][k]['kcal']))
    n = {c: len(v) for c, v in cls.items()}
    wt = sum(summ[c]/np.mean(meanref[c])*100*n[c] for c in summ)/sum(n.values())
    if verbose:
        print(tag, 'WTMAD~%.3f' % wt, ' '.join('%s=%.2f' % (c, v) for c, v in summ.items()))
        print('   ', ' '.join('%s:%+.1f(x%+.1f,c%+.1f)' % (k.split(':')[1], *v) for k, v in errs.items()))
    return wt, summ, errs
