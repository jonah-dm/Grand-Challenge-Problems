from hfeval import *
from diag import parts, M
C = load()
for nm in ['B_atom','B_cation','F_atom','F_cation','HF','HF_cation','O_atom','N_atom','C_atom','Li_atom','H_atom','OH','H2O']:
    z = C[nm]; A = args_of(z); m = (A[0]+A[1]) > 1e-11; A = [x[m] for x in A]; w = z['w'][m]
    with np.errstate(all='ignore'):
        P = parts(*A)
        I = lambda f: float(np.sum(w*np.nan_to_num(f)))
        # spin-resolved PBE correlation (Stoll): ss = Ec[n_s,0], os = rest
        ea = M._ec_pol(A[0], A[2]); eb = M._ec_pol(A[1], A[4])
    exl = z['Exl']; ex = z['Ex']
    row = []
    for i, l in enumerate('ab'):
        sl = I((1-P['w_'+l])*P['base_'+l] + P['w_'+l]*P['br_'+l])
        row.append((ex[i]-exl[i], sl, sl-(ex[i]-exl[i])))
    print(f"{nm:10s} SRx err a={row[0][2]*H2K:+7.1f} b={row[1][2]*H2K:+7.1f} | Ec_ss a={I(ea)*H2K:7.1f} b={I(eb)*H2K:7.1f} os={(I(P['cPBE'])-I(ea)-I(eb))*H2K:7.1f} gate={I(P['cgate'])*H2K:6.1f} tot={I(P['cPBE']+P['cgate'])*H2K:7.1f} ref(fc)={(REF['species'][nm]['E_tot']-REF['species'][nm]['hf_qz'])*H2K:7.1f}")
