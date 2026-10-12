You are a contestant in a benchmark. Working directory: /home/user/Grand-Challenge-Problems/benchmark/task2_xc

TASK: Construct a new exchange-correlation (XC) functional for Kohn-Sham DFT that is more
accurate than the best existing functionals, and that is accurate for reasons that generalize.

WHAT COUNTS. The benchmark rewards new physics or new construction methods, not re-fitting.
Acceptable routes include:
  - deriving and enforcing exact constraints on the XC functional that existing functionals
    violate (new sum rules, bounds, scaling or asymptotic conditions, fractional-charge or
    fractional-spin behavior, derivative-discontinuity behavior), and building a functional
    around them;
  - new functional forms or ingredients (new density-dependent variables, nonlocal
    constructions, density-matrix or orbital-dependent terms, machine-learned components)
    whose design is motivated by physics rather than by this benchmark;
  - machine-learned functionals whose architecture is built to satisfy known exact
    conditions by construction, trained on data you generate or derive yourself.
Not acceptable:
  - returning an existing libxc functional, or a linear mix of existing functionals;
  - re-fitting the parameters of an existing functional form to this benchmark;
  - any parameter chosen by fitting to the dev set. The dev set exists only so you can check
    that your code works and see roughly where you stand. Final scoring is on HIDDEN sets that
    include property classes and elements absent from the dev set; a functional tuned to the
    dev set will lose there.
A reviewer reads every submission. State in NOTES.md exactly what is new, what physical
argument justifies it, and which exact conditions your functional satisfies and violates.
Submissions judged to be re-fits or wrappers of existing functionals are disqualified.

BASELINES TO BEAT (dev-set WTMAD, lower is better): PBE 2.95, r2SCAN 1.28, B3LYP 0.97,
wB97M-V 0.56. See README.md, kit/dev_set.py, score.py and results/baselines/*.json.

RULES
- Write your submission to results/{MODEL}/submission.py defining
  eval_xc(xc_code, rho, spin, relativity, deriv, omega, verbose) in PySCF's
  dft.libxc.define_xc_ convention, plus HYB (exact-exchange fraction, 0.0 if none) and META
  (True if it uses tau). Calling libxc pieces as building blocks is allowed when the result
  is a genuinely new functional; see WHAT COUNTS.
- Test with:  NTHREADS=1 python score.py dev_set kit/dev_reference.json file:results/{MODEL}/submission.py --out results/{MODEL}/score.json
  (a smaller --basis such as def2-svp is fine for quick checks; final scoring uses def2-qzvp).
- Do NOT modify score.py, reference.py, README.md or anything under kit/. Do NOT read or write
  outside results/{MODEL}/. Do not install packages.
- Time budget: 45 minutes of wall time. Leave your best submission.py and score.json in place.
- Write results/{MODEL}/NOTES.md: what is new, the physical justification, exact conditions
  satisfied/violated, your dev-set score, and what you would try with more time. Report your
  final dev-set WTMAD in your last message.
