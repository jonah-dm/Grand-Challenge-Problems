# Organiser review of opus5t (round 5, theory track)

## Claims and verification
- C1 (E_xc[λn₁,0] = −λ²J[n₁], E_c = 0): proof via the Lieb ensemble functional is sound (for N ≤ 1 the
  ensemble mixes only 0- and 1-electron states, V_ee = 0, T_W convex). Same statement derived
  independently by opus5c (fractional-H uniqueness) and sonnet5c (from piecewise linearity). ✔
- T2 (F(s, z=1) = c·s⁻² is the unique Laplacian-free form satisfying C1; c = 5/9 from H): verified
  numerically by the organiser (see opus5c/REVIEW.md). ✔
- Shape-dependence counterexample: on a Gaussian one-orbital density the H-normed s⁻² form gives
  2.5× the exact exchange in the spin-scaled convention (organiser check; the agent reports 5×,
  which is the unpolarized-convention number). Either way the no-go holds. ✔
- Harness bug (vlapl zeros displacing vtau in PySCF's custom-functional path): confirmed in
  libxc.py (`extend([x for x in vxc if x is not None])` then `vstack(out)[:outlen]`). Fixed in
  score.py; a zeros-vlapl TPSS wrapper now reproduces built-in TPSS to 1e-13. Affected: every
  META submission from round 2 to round 5 (scores to be recomputed; effect ~0.02 on opus4).
- Self-consistent H₂⁺ at R = 6 bohr: −28.5 → −6.1 kcal/mol; H + H₂ barrier 7.2 → 9.06 (exact 9.6).
  Not independently recomputed; plausible given the construction and consistent with the
  barrier improvement seen for self-interaction-corrected forms in hidden set B.

## Classification
THEORY RESULT WORTH EXPERT REVIEW (same lemma as opus5c; three independent derivations in one
round). Functional: dev 1.55, not competitive on atomization energies (polarized atoms raised too
much); clear gain on the delocalization-sensitive quantities it targets. Candidate selection
disclosed; constant fixed by the H norm, not the dev set.
