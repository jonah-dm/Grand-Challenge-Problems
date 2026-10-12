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
