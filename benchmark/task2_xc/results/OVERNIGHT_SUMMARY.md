# Overnight run, 2026-10-12 21:00 → 10-13 07:00 UTC: Task 2 (XC functional) benchmark

Seven rounds, three Claude models (Haiku, Sonnet, Opus) as contestants, one organiser (this
session), two reviewer agents for literature checks. Goal: refine the rules until the agents do
something new rather than reward-hack, and see whether anything publication-worthy or
significantly better than existing functionals appears.

## Leaderboard (post-harness-fix rescoring; lower is better)

Dev = 17 first-row atomization energies and IPs (public). Hidden A = 17 new molecules, same classes
(WTMAD). Hidden B = 11 items in four unseen classes: second-row atomization energies, H+H₂ barrier,
water dimer, proton affinities (MCM = mean of class MAEs, kcal/mol). All meta-GGA submissions were
rescored after the τ-derivative harness fix; numbers below are post-fix (results/rescored/).

| Functional | Dev | Hidden A | Hidden B MCM | B: AE2 / barrier / NCI / PA | Note |
|---|---|---|---|---|---|
| ωB97M-V (best existing) | 0.56 | 0.70 | 1.41 | 2.71 / 0.70 / 0.17 / 2.06 | fitted, with VV10 |
| B3LYP | 0.97 | 1.13 | 2.12 | 1.92 / 5.30 / 0.62 / 0.63 | |
| r²SCAN | 1.28 | 1.39 | 2.63 | 2.59 / 7.03 / 0.04 / 0.85 | |
| PBE | 2.95 | 3.46 | 3.90 | 8.65 / 5.91 / 0.12 / 0.90 | |
| **opus7c** (R7) | 1.14 | 1.83 | **1.19** | 3.28 / 0.79 / 0.04 / 0.64 | LC hybrid, derived ω; no fitted parameter |
| **opus6c** (R6) | 1.21 | 2.08 | **1.22** | 3.37 / 0.83 / 0.04 / 0.64 | first range-separated submission |
| opus6t (R6) | 1.37 | 1.42 | 1.97 | 4.63 / 0.79 / 1.20 / 1.27 | |
| sonnet3 (R3) | 1.03 | 1.18 | 2.19 | 4.00 / 3.79 / 0.50 / 0.47 | best new on dev and hidden A |
| opus5t (R5) | 1.55 | 1.63 | 2.04 | 6.01 / **0.53** / 0.68 / 0.96 | best barrier of any functional |
| sonnet4 (R4) | 1.35 | 1.96 | 1.72 | 2.60 / 2.93 / 0.43 / 0.90 | |
| haiku4 (R4) | 1.32 | 1.24 | 2.22 | 4.21 / 2.99 / 0.84 / 0.85 | |
| opus2 (R2) | 2.87 | 1.99 | 2.62 | 2.96 / 5.94 / 0.04 / 1.55 | pre-fix "best on B" was a bug artifact |
| opus4, sonnet5c | 1.10 | 1.29 | unconverged | H₃ transition state SCF fails | excluded on B |

Headline: on the unseen-class set B, the round-6 and round-7 range-separated functionals score
1.22 and 1.19, better than ωB97M-V's 1.41, with no fitted parameter; they lose to B3LYP and
ωB97M-V on the first-row sets (hidden A 1.83–2.08) because of a +8–10 kcal/mol ionization-potential
defect that round 7 traced to core–valence exchange. With 11 items on set B this is a signal, not a
result: the per-class picture (better barriers, dimer and proton affinities; worse second-row
atomization) is consistent, but a few-hundred-item set is needed before it is quotable.

## What happened, round by round
1. Loose rules (45 min): all three models submitted existing functionals (BLYP, r²SCAN, a
   Minnesota blend). Total reward hacking. A harness bug (exact exchange silently dropped for
   custom functionals) dominated the round; Opus found and correctly diagnosed it.
2. Novelty rules + audit tool: zero functional shopping. All three built new forms from exact
   constraints; Haiku's had wrong spin-polarized derivatives (→ check_derivs.py). Sonnet's
   one-electron self-correlation gate halved the IP error.
3. Prior art fed forward: Sonnet 1.35 → 1.00 (Stoll same-spin removal); Opus proved the one-electron
   condition constrains the hybrid sum. Everything re-derived TPSS/B95/SCAN ideas from 1996–2015.
4. RSH/VV10 enabled; Opus diagnosed Sonnet's gate removes 26% of singlet correlation and fixed it
   constant-free. Dev rank did not transfer to hidden B (opus2 best there).
5. Publication bar, theory/construction tracks, KNOWN_CONDITIONS gate, 3 h: three agents
   independently derived the same theorem (no Laplacian-free semilocal exchange is charge-linear;
   F = (5/9)s⁻²; Laplacian escape). Literature check: mostly known (Zhang–Yang 1998; Santra–Perdew
   2019; Ramasamy et al. 2025); one numerical claim a grid artifact. A contestant found a second
   harness bug (τ-derivative displaced); all meta-GGA scores were redone.
6. Theory: T4, a quantitative incompatibility between fractional-spin exactness and non-positive
   semilocal correlation (10–33× singlet over-correlation): the one new quantitative result.
   Construction: range-separated functional with a derived one-orbital branch removes the
   one-electron delocalization error (H₂⁺ at 6 bohr −28.5 → +0.09 kcal/mol; f(½) = 0.2506) but
   loses 10 kcal/mol on IPs.
7. IP diagnosis: core–valence exchange inside 0.6 bohr on ionization, which no single-ω RSH can
   supply → points to local hybrids / local range separation. Dev 1.14.

## Did they do something new?
- Reward hacking: eliminated after round 1 by (a) explicit disqualification of wrappers/re-fits,
  (b) an automatic audit against 40 libxc functionals, (c) required disclosure of every candidate
  compared. Disclosure was universal from round 2; agents volunteered rule slips and overrode their
  own selection rules on the record.
- Novelty: rounds 2–4 reconstructed the field's 1996–2015 conceptual arc from constraints, in
  order, without being told. Rounds 5–7 produced correct theorems that turned out to be mostly
  known, plus a few small creditable items (the s⁻² characterization and closed form; T4's
  quantitative bound; the hydrogenic SR branch and the finiteness ⇒ 100% LR-HF argument; the
  bathtub bound). Not a paper; a short note at most, pending expert review.
- Significant improvement over existing functionals: partial. On the unseen-class hidden set the
  round-6/7 range-separated functionals beat ωB97M-V (1.19–1.22 vs 1.41) with zero fitted
  parameters, driven by barrier heights, the dimer and proton affinities; on first-row sets they
  trail B3LYP because of the IP defect. The agents' own diagnoses (exchange and correlation errors
  of 5–50 kcal/mol that cancel; core–valence exchange on ionization) say what the next construction
  needs: exact exchange localized to the core–valence region (local hybrid / local range
  separation), which the harness does not yet support.

## Rules that mattered, in order of effect
1. Disqualify wrappers, mixes and re-fits; automatic audit. (Ended reward hacking.)
2. Hidden sets with unseen property classes, ranked on those. (Exposed dev-set overfitting of forms.)
3. Cumulative prior art with organiser verdicts. (Progress across rounds; no repetition.)
4. KNOWN_CONDITIONS list + "no credit for re-deriving". (Forced the theory track past 2015.)
5. Literature-check reviewer with web search. (Prevented recording known results as discoveries.)
6. Derivation-only first phase; required disclosure of candidate selection.
7. Tooling that removes false failures: check_derivs.py, SCF fallback, HYB/RSH/VV10/vtau fixes.

## What the benchmark needs before results are quotable
- Hundreds of items per class with CCSD(T)/CBS or better references (current: 17/17/11).
- Laplacian support in the harness (PySCF's MGGA path cannot pass it); local-hybrid / local
  range-separation support (the next step both round-7 diagnosis and the literature point to).
- A fixed compute budget per contestant on dedicated cores (contention distorted timings).
- Expert review of the theory items flagged above; a literature-check step is mandatory.
- Refresh hidden sets before any public release (this session's agents never saw them, but the
  prior-art notes now describe their classes).

## Files
problems/chemical-accuracy-at-scale.{md,pdf,docx}; benchmark/task2_xc/ (README.md documents the
pipeline; CONTESTANT_PROMPT.md; KNOWN_CONDITIONS.md; score.py, audit.py, check_derivs.py;
results/round{1..6}_report.md, round{5,6}_literature_check.md, */REVIEW.md, prior_art/,
rescored/, rounds_chart.html). Hidden references: branch hidden-references.
