from t_cand import *
for tag, f, a, w, cap in [('S5 cap1.299', fx_ch, 0, 0.5256070335369912, 1.299259029906918), ('S5t total-z', fx_tot, 0, 0.5256070335369912, 1.299259029906918)]:
    M.OMEGA = w; M.F_LO = cap
    for nm in C: C[nm]['Exl'] = get(nm, w)
    score(species_terms(C, f, fc, a=a), tag=tag)
