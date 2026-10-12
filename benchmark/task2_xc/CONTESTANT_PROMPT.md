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
    conditions by construction, trained on data you generate or derive yourself (e.g. exact
    exchange energies of densities you compute, model systems, or exact-condition samples).
Suggested, not required: start from a specific failure of the best baseline on the dev set
(e.g. one-electron self-interaction in atomic IPs; X-H bond underbinding), identify the exact
condition whose violation explains it, and build the functional around satisfying it.
Not acceptable:
  - returning an existing libxc functional, or a linear mix of existing functionals;
  - re-fitting the parameters of an existing functional form to this benchmark;
  - any continuous parameter chosen by fitting to the dev set. Choosing between a handful of
    physically motivated candidate forms using dev-set scores is allowed, but you must list
    every candidate you compared and its score in NOTES.md; undisclosed selection is
    disqualifying. The dev set exists only so you can check
    that your code works and see roughly where you stand. Final scoring is on HIDDEN sets that
    include property classes and elements absent from the dev set; a functional tuned to the
    dev set will lose there.
An automated audit (audit.py) is run on every submission: it fits your functional's energy density
as a combination of ~20 common libxc functionals on random densities and checks the uniform-gas
limit, the Lieb-Oxford bound and negativity. WRAPPER and MIX are disqualified. NEAR-EXISTING
is not automatically disqualified: a functional derived independently from exact constraints can
legitimately land close to an existing one (constraints alone nearly rebuild r2SCAN). It goes
to the reviewer, who checks that the derivation is documented and the constants are derived,
not copied. DISTINCT is scored directly. Run it yourself:
  OMP_NUM_THREADS=1 python audit.py file:results/{MODEL}/submission.py
A reviewer also reads every submission.
Before any SCF run, verify your derivatives (wrong spin-polarized vxc is the most common failure):
  OMP_NUM_THREADS=1 python check_derivs.py file:results/{MODEL}/submission.py State in NOTES.md exactly what is new, what physical
argument justifies it, and which exact conditions your functional satisfies and violates.
Submissions judged to be re-fits or wrappers of existing functionals are disqualified.

PRIOR ART. Earlier rounds of this benchmark produced the attempts in results/prior_art/*.md (read
them first; *3.md are the latest). None beat B3LYP (0.97) yet. Frontier so far:
  - sonnet3 (dev WTMAD 1.00): 25% hybrid meta-GGA; correlation = E_c^PBE[n_a,n_b] - sum_s z_s^2
    E_c^PBE[n_s,0] (Stoll partition: no same-spin correlation in a channel holding one orbital).
  - opus3 (1.59): joint one-electron condition on (1-a)E_x^SL + E_c; H atom exact by construction;
    showed no PBE-form iso-orbital exchange can meet both H and He norms.
Known failure modes: one-orbital exchange norms that underbind X-H bonds; gates that kill genuine
opposite-spin correlation in open-shell atoms (N2, C2H2 overbind); wrong spin-polarized
derivatives; constraint-only exchange rebuilding r2SCAN. Build on or go beyond these; do not
repeat them. Ideas named but untried: a per-spin same-spin removal paired with an opposite-spin
correlation built for it; an iso-orbital exchange meeting both H and He norms (s^-1/2 tail);
range separation for IPs/barriers; nonlocal correlation for noncovalent interactions.

HARNESS. Besides HYB and META, submission.py may define RSH = (omega, alpha, beta) for a
range-separated hybrid (alpha*HF + beta*short-range HF with erfc(omega r)/r; eval_xc then gets
omega and must return the matching semilocal part) and NLC = (b, C) to have the harness add a
VV10 nonlocal correlation term. score.py reports WTMAD and MCM (plain mean of class MAEs).

HIDDEN TEST CLASSES. Final ranking uses hidden sets with: first-row and SECOND-ROW atomization
energies, ionization potentials, a reaction BARRIER HEIGHT, a NONCOVALENT interaction and PROTON
AFFINITIES. The dev set has only first-row AEs and IPs. Design for all of these.

BASELINES TO BEAT (dev-set WTMAD, lower is better): PBE 2.95, r2SCAN 1.28, B3LYP 0.97,
wB97M-V 0.56. See README.md, kit/dev_set.py, score.py and results/baselines/*.json.

RULES
- Write your submission to results/{MODEL}/submission.py defining
  eval_xc(xc_code, rho, spin, relativity, deriv, omega, verbose) in PySCF's
  dft.libxc.define_xc_ convention, plus HYB (exact-exchange fraction, 0.0 if none) and META
  (True if it uses tau). Calling libxc pieces as building blocks is allowed when the result
  is a genuinely new functional; see WHAT COUNTS.
- Test with:  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NTHREADS=1 python score.py dev_set kit/dev_reference.json file:results/{MODEL}/submission.py --out results/{MODEL}/score.json
  (a smaller --basis such as def2-svp is fine for quick checks; final scoring uses def2-qzvp).
- Do NOT modify score.py, reference.py, README.md or anything under kit/. Do NOT read or write
  outside results/{MODEL}/. Do not install packages.
- Time budget: 90 minutes of wall time. A full dev-set run at def2-qzvp takes ~10 min on one
  thread; use def2-svp for iteration and reserve one QZ run for your final submission. Leave your best submission.py and score.json in place.
- Write results/{MODEL}/NOTES.md: what is new, the physical justification, exact conditions
  satisfied/violated, your dev-set score, and what you would try with more time. Report your
  final dev-set WTMAD in your last message.
