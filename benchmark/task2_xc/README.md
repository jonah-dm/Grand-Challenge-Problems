# Task 2 benchmark: design a better XC functional

## Task given to the model

See `CONTESTANT_PROMPT.md`. In short: construct a *new* XC functional, more accurate than the
best existing ones, for reasons that generalize. Rewarded: new exact constraints and
functionals built to satisfy them; new functional forms or ingredients motivated by physics;
machine-learned functionals that satisfy exact conditions by construction. Disqualified:
wrapping or mixing existing functionals, re-fitting existing forms, or any parameter fit to
the dev set. Every submission is read by a reviewer.

## Reference data
CCSD(T)/CBS (HF/QZ + T,Q-extrapolated correlation), frozen core, nonrelativistic, fixed
geometries, no ZPE. Scored at the def2-QZVP basis. `kit/dev_reference.json`.

## Scoring
1. Hidden set A: same classes as the dev set, different molecules.
2. Hidden set B: classes absent from the dev set (second-row atomization energies, barrier
   heights, noncovalent interactions, proton affinities).
3. Reviewer classification: new constraint / new form / ML-with-constraints / re-fit / wrapper.
   Re-fits and wrappers are disqualified.
Metric per set: MAE per class and WTMAD (each class scaled by its mean |reference|), lower is
better. A functional is ranked by its hidden-set B score and must not be worse than the best
baseline on any class by more than a stated tolerance. Must converge on every species.

## Files
- `kit/dev_set.py`, `kit/geom.py` — public development set
- `reference.py` — regenerate reference energies (`python reference.py dev_set out.json`)
- `score.py` — score a functional: `python score.py dev_set kit/dev_reference.json pbe`
- `results/baselines/` — baseline functionals
- `results/<model>/` — model submissions and scores

## Evaluation pipeline (as run overnight 2026-10-12/13)
1. Contestant agent in a sandbox with CONTESTANT_PROMPT.md, KNOWN_CONDITIONS.md, the public dev
   set, score.py, audit.py, check_derivs.py and prior_art/ (notes from earlier rounds, with
   organiser verdicts). Fixed time budget, one CPU thread, derivation-only first phase.
2. Automatic gates: check_derivs.py (vxc vs finite differences, both spin paths); audit.py
   (wrapper/mix detection against ~40 libxc functionals; UEG, Lieb–Oxford, negativity; a
   spin-polarized one-orbital probe). WRAPPER/MIX disqualify; NEAR-EXISTING goes to review.
3. Scoring: dev set, hidden set A (same classes, new molecules), hidden set B (unseen classes).
   Ranking on hidden set B by MCM. The hidden sets live on the `hidden-references` branch.
4. Organiser review (REVIEW.md per submission): numerical re-verification of the central claim,
   check of disclosed candidate selection, classification against KNOWN_CONDITIONS.md.
5. Literature check by a separate reviewer agent with web search for any claimed theory result
   (results/round5_literature_check.md is the template). Round 5 showed this step is essential:
   three agents' "new theorem" was mostly known, and one quantitative claim was an artifact.
6. Rounds are cumulative: each round's notes and verdicts feed the next round's prior art.
Harness bugs found and fixed along the way: HYB ignored for custom functionals (round 1);
Newton fallback needing fxc (round 2); vtau displaced by a zeros vlapl (round 5, found by a
contestant).
