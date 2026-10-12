# haiku submission notes

## Result
- Submission: `results/haiku/submission.py` (pure GGA, HYB = 0.0, META = False): Becke-88 exchange + LYP correlation (BLYP-type), via `pyscf.dft.libxc.eval_xc`.
- Dev-set WTMAD (def2-qzvp, scored with score.py): **1.68** (AE MAE 3.84, IP MAE 5.12). `score.json` and `score.log` are in this directory.
- This is better than PBE (2.95) but worse than B3LYP (0.97), r2SCAN (1.28) and wB97M-V (0.56). It does NOT beat the best baseline, so the goal was not met.

## What was tried
1. PBE0-family hybrids, `(1-a)*PBE_x + PBE_c` with HYB = a (a = 0.20, 0.25, 0.30). Dev WTMAD 8.66 (a=0.20) and 11.22 (a=0.25); a=0.30 was killed before finishing. These numbers are NOT valid hybrid results. `define_xc_` with a callable description ignored the `hyb` argument in this PySCF build (H2 energy identical for hyb=0 and hyb=0.2 in a direct check), so exact exchange was never added and only the semilocal exchange was scaled down. Scan files are in `results/haiku/scan/`.
2. BLYP-type pure GGA (B88 + LYP): WTMAD 1.68. This is the submitted candidate.

## Why it was not better
- Hybrid functionals (the source of B3LYP's and wB97M-V's gains) could not be evaluated through score.py's `define_xc_` path, so that direction was blocked for this run.
- Full dev-set runs take roughly 10 minutes each at def2-qzvp with NTHREADS=1 under parallel load, which left room for only one clean pure-GGA candidate.

## What I would try with more time
- Get the exact-exchange fraction honored: check how this PySCF build handles `hyb` for callable xc (e.g. set `mf.xc` to a libxc string such as a PBE0 or B3LYP mix, or inspect `is_hybrid_xc`), then rerun the PBE0-family scan with a in {0.20, 0.25}.
- Consider a B3LYP-like mix (B88 + LYP + VWN, 20% HF) once hybrids are wired correctly, since B3LYP already scores 0.97 on this set.
- Scan the semilocal exchange/correlation mixing coefficients on the dev set, while keeping the change physically motivated and checking the hidden-set transfer.
