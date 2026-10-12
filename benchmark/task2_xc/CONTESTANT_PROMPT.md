You are a contestant in a benchmark. Working directory: /home/user/Grand-Challenge-Problems/benchmark/task2_xc

TASK: Design an exchange-correlation (XC) functional for Kohn-Sham DFT that scores better than
the best baseline on this benchmark. Read README.md, kit/dev_set.py, score.py and
results/baselines/*.json first (baseline summary: PBE WTMAD 2.95, r2SCAN 1.28, B3LYP 0.97,
wB97M-V 0.56; lower is better. WTMAD combines atomization-energy and ionization-potential MAE).

RULES
- Write your submission to results/{MODEL}/submission.py. It must define
  eval_xc(xc_code, rho, spin, relativity, deriv, omega, verbose) in PySCF's dft.libxc.define_xc_
  convention, plus HYB (float exact-exchange fraction, 0.0 if none) and META (True if it uses tau).
  You may build on libxc functionals via pyscf.dft.libxc.eval_xc inside your eval_xc.
- Test with:  NTHREADS=1 python score.py dev_set kit/dev_reference.json file:results/{MODEL}/submission.py --out results/{MODEL}/score.json
  (you may use a smaller --basis such as def2-svp for quick iteration; the final score uses def2-qzvp).
- Do NOT modify score.py, reference.py, README.md or anything under kit/. Do NOT look at or write
  outside results/{MODEL}/. Do not install packages.
- Time budget: 45 minutes of wall time. Leave your best submission.py and score.json in place.
- Final scoring will use a HIDDEN set of different molecules from the same classes, so
  prefer physically motivated changes over fitting parameters to the 17 dev-set items.
- When done, write results/{MODEL}/NOTES.md: what you did, why, your dev-set score, and what
  you would try with more time. Report your final dev-set WTMAD in your last message.
