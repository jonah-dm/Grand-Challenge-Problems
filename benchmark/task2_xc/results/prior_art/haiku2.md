# haiku2 — NOTES

## Status
Submission is an audit-DISTINCT spin-scaled meta-GGA. It does NOT yet produce a dev-set score:
under plain SCF (score.py path), 7 of 24 dev species fail to converge (B, C, O, F, OH, F_cation,
HF_cation) and score.py then needs second derivatives (fxc) for its Newton fallback, which I did
not implement. No WTMAD is reported. See `convergence_check_svp.log` (def2-svp).

## What is new (and what is not)
- Exchange enhancement F(s) = 1 + kappa*(1 - exp(-mu s^2/kappa)), mu = 10/81 (exact gradient
  expansion coefficient), kappa = 0.804 (ceiling 1.804). Saturating exponential rather than
  PBE's rational (Pade) saturation. This alone is nearly PBE-like (audit mix R^2 ~ 0.9999).
- Iso-orbital factor G(z) = 1 - (1/4) z/(1+z), z = tau_W/tau. The 1/4 is a chosen constant,
  not fitted. It is the ingredient that makes the audit verdict DISTINCT (mix R^2 0.938).
- Correlation: PW92 spin-interpolated LDA plus H = gamma*phi^3 ln(1 + (beta/gamma) t^2 exp(-A t^2)),
  standard beta, gamma; the exponential cutoff replaces PBE's rational one.

## Exact conditions
- Satisfied: uniform-gas limit (audit UEG relative deviation ~1e-15), Lieb-Oxford pointwise
  check (0 violations in audit sample), non-positivity of e_xc (0 violations), smoothness.
- Violated: exact one-electron exchange (z = 1 limit). A stronger G(1)=1/2 scaling was tested
  and rejected: H atom energy -0.358 Ha vs exact -0.5. The mild factor here does not enforce it.
- Not checked: spin-scaling beyond the code path, size-consistency, derivative discontinuity.

## Dev-set score
Not obtained. Plain SCF did not converge for 7 open-shell species. The z-dependent factor is
the cause: with G = 1 (z term removed), C atom converges to -37.620 Ha; with G = 1 - z/4(1+z) it
does not. A bounded variant z = tau_W/(tau + tau_W) did not fix it.

## What I would try with more time
1. Implement fxc (second derivatives) so the Newton fallback works, or make the z dependence
   SCF-stable (smooth the tau dependence, or use a tau-based form with better-conditioned
   derivatives in low-density spin channels).
2. Replace the arbitrary 1/4 with a constraint-derived one-orbital form.
3. Score at def2-svp, then one def2-qzvp run for the final submission.
