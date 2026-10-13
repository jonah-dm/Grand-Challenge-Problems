"""Diagnosis of opus6c IP errors at UHF/cc-pVQZ densities: exchange vs correlation, branch vs base, by species."""
import sys, json, os, importlib.util
import numpy as np
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../kit')))
from pyscf import gto, scf, dft
from geom import to_str
import dev_set as S
H2K = 627.5094740631
spec = importlib.util.spec_from_file_location('o6', os.path.join(os.path.dirname(__file__), '../../opus6c/submission.py'))
M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M)
ref = json.load(open(os.path.join(os.path.dirname(__file__), '../../../kit/dev_reference.json')))['species']
W = M.OMEGA

def parts(ra, rb, saa, sab, sbb, ta, tb):
    n = ra + rb; sig = saa + 2*sab + sbb
    z = M._zsat(sig/(8*n*(ta+tb))); zeta = (ra-rb)/n
    out = {}
    for lab, r, s, t in (('a', ra, saa, ta), ('b', rb, sbb, tb)):
        n2 = 2*r; kf = (3*np.pi**2*n2)**(1/3); s2 = 4*s/(4*kf**2*n2**2)
        eb = 0.5*M._ex_base_sr(n2, s2, z)
        tu_ = 0.3*(6*np.pi**2)**(2/3)*r**(5/3)*M.ZREG
        zs_ = M._zsat((s/(8*r)+tu_)/(t+tu_))
        w = zeta**2*zs_**2
        e1 = 0.5*M._ex_branch(r, s)
        out['base_'+lab] = eb; out['br_'+lab] = e1; out['w_'+lab] = w
        out['zs_'+lab] = zs_
    out['zeta'] = zeta
    # correlation
    rs = (3/(4*np.pi*n))**(1/3); ec = M._pw92(rs, zeta)
    phi = 0.5*((1+zeta)**(2/3)+(1-zeta)**(2/3)); kf = (3*np.pi**2*n)**(1/3); ks2 = 4*kf/np.pi
    t2 = sig/(4*phi**2*ks2*n**2); phi3 = phi**3
    A = M.BETA/M.GAMMA/(np.exp(-ec/(M.GAMMA*phi3))-1)
    H = M.GAMMA*phi3*np.log(1+M.BETA/M.GAMMA*t2*(1+A*t2)/(1+A*t2+A*A*t2*t2))
    out['cPBE'] = n*(ec+H)
    za_ = M._zsat(saa/(8*ra*ta)); zb_ = M._zsat(sbb/(8*rb*tb))
    gate = 0
    for r, sg_s, t_s, zo in ((ra, saa, ta, zb_), (rb, sbb, tb, za_)):
        zs = M._zsat(sg_s/(8*r*t_s))
        gate = gate - zs**2*(1-zo**M.PW*(1-zeta**2))*M._ec_pol(r, sg_s)
    out['cgate'] = gate
    return out

def run(name, basis='cc-pvqz'):
    atoms, charge, spin = S.SPECIES[name]
    mol = gto.M(atom=to_str(atoms), basis=basis, charge=charge, spin=spin, verbose=0)
    mf = scf.UHF(mol); mf.conv_tol = 1e-10; mf.max_cycle = 200; mf.run()
    dm = mf.make_rdm1()
    if dm.ndim == 2: dm = np.array([dm/2, dm/2])
    Kf = mf.get_k(mol, dm)
    Kl = mf.get_k(mol, dm, omega=W)
    Ex = -0.5*sum(np.einsum('ij,ji', dm[s], Kf[s]) for s in range(2))
    Exl = -0.5*sum(np.einsum('ij,ji', dm[s], Kl[s]) for s in range(2))
    g = dft.gen_grid.Grids(mol); g.level = 4; g.build()
    ni = dft.numint.NumInt()
    ao = ni.eval_ao(mol, g.coords, deriv=1)
    rho = [ni.eval_rho(mol, ao, dm[s], xctype='MGGA', with_lapl=False) for s in range(2)]
    ra = np.maximum(rho[0][0], 1e-14); rb = np.maximum(rho[1][0], 1e-14)
    ga, gb = rho[0][1:4], rho[1][1:4]
    P = parts(ra, rb, (ga**2).sum(0), (ga*gb).sum(0), (gb**2).sum(0), np.maximum(rho[0][4], 1e-20), np.maximum(rho[1][4], 1e-20))
    mask = (ra+rb) > 1e-11
    I = lambda f: float(np.sum((g.weights*np.nan_to_num(f))[mask]))
    exsl = sum((1-P['w_'+l])*P['base_'+l] + P['w_'+l]*P['br_'+l] for l in 'ab')
    res = dict(E_hf=mf.e_tot, Ex_hf=Ex, Ex_lr=Exl, Ex_sr=Ex-Exl,
               Exsl=I(exsl), Exsl_base=I(P['base_a']+P['base_b']),
               dbranch=I(sum(P['w_'+l]*(P['br_'+l]-P['base_'+l]) for l in 'ab')),
               EcPBE=I(P['cPBE']), Ecgate=I(P['cgate']),
               Ec_ref=ref[name]['E_tot']-ref[name]['hf_qz'], hf_ref=ref[name]['hf_qz'])
    res['Ec'] = res['EcPBE'] + res['Ecgate']
    res['dEx'] = res['Exsl'] - res['Ex_sr']          # SR semilocal error vs exact SR exchange
    res['dEc'] = res['Ec'] - res['Ec_ref']
    # radial bins (atoms) of branch delta and exchange
    r = np.linalg.norm(g.coords - mol.atom_coords()[0], axis=1)
    bins = [0, 0.5, 1, 1.5, 2, 3, 4, 6, 100]
    res['branch_by_r'] = [I(np.where((r >= bins[i]) & (r < bins[i+1]), sum(P['w_'+l]*(P['br_'+l]-P['base_'+l]) for l in 'ab'), 0)) for i in range(len(bins)-1)]
    res['gate_by_r'] = [I(np.where((r >= bins[i]) & (r < bins[i+1]), P['cgate'], 0)) for i in range(len(bins)-1)]
    return res

if __name__ == '__main__':
    names = sys.argv[1:] or ['B_atom', 'B_cation', 'F_atom', 'F_cation', 'HF', 'HF_cation']
    out = {}
    for nm in names:
        out[nm] = run(nm)
        d = out[nm]
        print(nm, ' '.join(f'{k}={v:.6f}' for k, v in d.items() if not isinstance(v, list)), flush=True)
        print('   branch_by_r', np.round(np.array(d['branch_by_r'])*H2K, 2), ' gate_by_r', np.round(np.array(d['gate_by_r'])*H2K, 2))
    json.dump(out, open(os.path.join(os.path.dirname(__file__), 'diag_%s.json' % '_'.join(names)[:60]), 'w'), indent=1)
