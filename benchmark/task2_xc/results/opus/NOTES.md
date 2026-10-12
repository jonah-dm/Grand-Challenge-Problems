# opus — Task 2 (XC functional) notes

## IMPORTANT: harness finding — `HYB` is silently ignored by score.py
`score.py` calls `mf.define_xc_(sub.eval_xc, ..., hyb=HYB)` but leaves `mf.xc` at its
default `'LDA,VWN'`. In PySCF 2.14, `uks.get_veff`/`rks.get_veff` decide whether to add exact
exchange with `ni.libxc.is_hybrid_xc(ks.xc)` — i.e. on the *string* `'LDA,VWN'`, not on the
`hybrid_coeff` that `define_xc_` installs. So for any `file:` submission the exact-exchange
term is never added, whatever `HYB` is. Verified on the H atom (def2-SVP, UKS):
`mf.xc='M062X'` -> -0.4974 Eh; `define_xc_(eval_M062X, 'MGGA', hyb=0.54)` -> -0.3408 Eh
(the missing ~0.157 Eh is 0.54 x the exact exchange of H).
The baselines (B3LYP, wB97M-V) are run with `mf.xc=<string>`, so they *do* get exact exchange;
file submissions cannot. Fix (for the organisers, I did not modify score.py): also set
`mf.xc` to something hybrid, or patch `is_hybrid_xc`, or compare via a correct code path.
I did not monkeypatch PySCF from inside submission.py, because that would amount to
changing the scorer.

## What I did
1. Screened global hybrids at def2-TZVP via libxc strings (these get exact exchange because
   they go through `mf.xc`), then estimated QZ values by adding B3LYP's per-item TZ->QZ shift:
   est. WTMAD at QZ: M06-2X ~0.72, PW6B95 ~0.73, MN15 ~0.76, B3LYP 0.97, B97-1 ~1.24.
   None clearly beat wB97M-V (0.56); M06-2X/MN15 fix the B IP (high exact exchange reduces
   self-interaction error). My planned submission was an equal-weight ensemble
   0.5*M06-2X + 0.5*MN15 (HYB 0.49), whose TZ-based estimate had AE MAE ~1.0, IP MAE ~1.4,
   because their AE errors are partly anticorrelated (LiH, F2, CH2).
2. Discovered the HYB bug above when the QZ run of that submission gave H atom = -0.36 Eh.
3. Because only semilocal functionals can be scored correctly through `file:`, I screened
   meta-GGAs at def2-TZVP (MN12-L, MN15-L, B97M-V semilocal part, M06-L, revM06-L; partial
   runs). After a rough TZ->QZ basis correction MN12-L looked best for AE (it is also known to
   be one of the more accurate local functionals for atomization energies), so I launched a
   def2-QZVP run of MN12-L. Its SCF on the N atom did not finish in >8 min (Minnesota local
   functionals are known for poor SCF behaviour), so with the budget running out I fell back to
   the most robust/accurate semilocal option, r2SCAN via the `file:` path (`HYB = 0.0`,
   `META = True`). This reproduces the r2SCAN baseline exactly; it does NOT beat B3LYP or
   wB97M-V. Aborted MN12-L log: sl/MN12_L_qz_aborted.log; TZ screens in scr/ and sl/.

## Dev-set score (def2-QZVP)
See score.json (WTMAD reported in the final message).

## With more time
- Get the scorer fixed, then submit a high-exact-exchange meta-GGA hybrid (M06-2X / MN15
  ensemble, or a re-fit B97-type hybrid with ~45-50% exact exchange to tame IP SIE).
- With HYB unavailable: a semilocal meta-GGA with an iso-orbital indicator correction that
  reduces one-electron self-interaction (B atom IP, F atom IP are the worst items for r2SCAN),
  e.g. r2SCAN with a modified alpha->0 limit, fit only to exact constraints (H atom, uniform gas).
