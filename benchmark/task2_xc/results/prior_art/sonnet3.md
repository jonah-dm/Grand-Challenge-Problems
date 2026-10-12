# sonnet3: spin-resolved one-orbital correlation gate (Stoll partition) in a 25% hybrid meta-GGA

Final dev-set WTMAD (def2-qzvp): **1.00** (AE MAE 2.22, IP MAE 3.56 kcal/mol; all species converge). Baselines: r2SCAN 1.28, B3LYP 0.97, sonnet2 1.35.
Audit: DISTINCT (mix R2 ~0.998, closest single PBE0 R2 0.970). check_derivs: OK (complex-step derivatives).

## What is new
Exact condition used: a spin channel occupied by a single orbital (z_s = tau_W,s/tau_s = 1) has no same-spin correlation
(one electron per spin orbital). Semilocal PBE correlation violates this (self-correlation of H, He, H2 same-spin part).
Construction: Stoll partition E_c = E_c^opp + sum_s E_c^ss, with E_c^ss = E_c^PBE[n_s,0], and the same-spin part is gated:
  E_c = E_c^PBE[na,nb] - sum_s z_s^2 * E_c^PBE[n_s, 0]
- Exactly zero for H atom and any fully polarised one-orbital density (as in sonnet2's gate), AND the same-spin part vanishes for
  two-electron singlets (He, H2, bonding regions), which sonnet2's zeta^2 z^2 gate (polarised only) cannot do.
- Uniform gas limit untouched (z_s -> 0). Not a B95 copy: B95 uses (tau_s - tau_W,s)/tau_unif factors in a fitted form; here it is a gate on PBE's own Stoll component, no new constants.
- Exchange (inherited idea from sonnet2, re-derived): z-blend (weight z^2, z total) of a PBE-form GGA branch (z->0) and a one-orbital branch with
  F_x ~ s^(-1/2) tail and max F_x = 1+kappa (Lieb-Oxford saturated), constant c analytic. 25% exact exchange (1/n, n=4 adiabatic connection), semilocal x scaled 0.75.
- Constants: kappa = 0.804 (LO), beta = 0.06672, gamma = (1-ln2)/pi^2, mu = beta*pi^2/3 = 0.2195 (xc-level second-order gradient / linear-response cancellation of PBE, derived
  from beta, not fitted). All derived; nothing fitted to data.

## Exact conditions
Satisfied: uniform gas; spin scaling of exchange; Lieb-Oxford (exchange F_x <= 1.804, audit 0 violations); negativity; E_c = 0 for one-orbital polarised densities;
same-spin E_c = 0 for single-orbital spin channels (new); s^-1/2 exchange tail in one-orbital regions; high-density correlation limit (PBE H).
Violated: exchange exact-GE2 coefficient 10/81 (uses mu_PBE; xc-level condition instead); one-electron exchange not exact pointwise; no derivative discontinuity or
fractional-charge condition beyond what 25% EXX gives; no dispersion (NCI will be weak); opposite-spin part is PBE's (not gated); PBE correlation low-density limit inherited.

## Candidates compared (all def2-svp dev WTMAD; for reference svp r2SCAN 2.11, PBE0 2.12, sonnet2 2.21)
| exchange GGA-branch mu | same-spin gate | WTMAD svp |
|---|---|---|
| 10/81 (GE2) | none | 3.81 |
| 10/81 (GE2) | Stoll z^2 gate | 2.53 |
| 0.2195 (PBE) | none | 2.74 |
| 0.2195 (PBE) | Stoll z^2 gate (SUBMITTED) | 1.99 |
Selection among these four physically motivated forms used dev scores (disclosed). No continuous parameter was tuned. Only the final form was run at QZ (1.00).
The gate improves both mu choices (-1.3, -0.75) so it is not a dev artifact of one choice. Env switches S3_MU/S3_GATE in submission.py exist only for this ablation (defaults = submitted form).

## Weaknesses seen
H2 (-6) and LiH (-5.8) underbound, F2 -3.9, HF -5.0; B IP +5.4. Second-row, barriers, NCI not tested.

## With more time
Gate the opposite-spin part with an exact two-electron condition (E_c of He/H2 exact relation); derive an exchange hole-based one-orbital exchange to replace the 25% constant;
range separation for derivative discontinuity; non-local dispersion term; test z_s^p exponents from a He/H2 model density rather than the plain z^2.
