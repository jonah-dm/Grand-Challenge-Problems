from t_cand import *
for tag, f, a, w, cap in [('S6 cap1.299 a=1/4', fx_ch, 0.25, 0.5256070335369912, 1.299259029906918)]:
    M.OMEGA = w; M.F_LO = cap
    for nm in C: C[nm]['Exl'] = get(nm, w)
    score(species_terms(C, f, fc, a=a), tag=tag)
