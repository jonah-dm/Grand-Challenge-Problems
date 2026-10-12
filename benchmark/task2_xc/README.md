# Task 2 benchmark: design a better XC functional

## Task given to the model

> Design an exchange–correlation functional that beats the baseline functionals on the
> benchmark below. Submit a Python file defining `eval_xc(...)` in PySCF's `define_xc_`
> convention (plus `HYB` for exact-exchange fraction and `META=True` if it uses tau).
> You may run `score.py` on the public dev set as often as you like within the budget.
> Final scoring uses a hidden set of the same classes (atomization energies, ionization
> potentials) that you have not seen.

## Reference data
CCSD(T)/CBS (HF/QZ + T,Q-extrapolated correlation), frozen core, nonrelativistic, fixed
geometries, no ZPE. Scored at the def2-QZVP basis. `kit/dev_reference.json`.

## Metric
MAE per class and WTMAD (each class scaled by its mean |reference|). Lower is better.
A submission must converge on every species.

## Files
- `kit/dev_set.py`, `kit/geom.py` — public development set
- `reference.py` — regenerate reference energies (`python reference.py dev_set out.json`)
- `score.py` — score a functional: `python score.py dev_set kit/dev_reference.json pbe`
- `results/baselines/` — baseline functionals
- `results/<model>/` — model submissions and scores
