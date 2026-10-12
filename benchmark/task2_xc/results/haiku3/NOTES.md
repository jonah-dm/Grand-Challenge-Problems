# haiku3 NOTES

Submission: `results/haiku3/submission.py` (meta-GGA, HYB = 0, META = True).
Final dev-set score (def2-qzvp, all species converged): **WTMAD 1.94** (AE MAE 4.52, IP MAE 5.42 kcal/mol).
Audit (`audit.py`): UEG rel. dev 3.5e-7, LO violations 0, negativity 0, mix R2 0.99993, closest single
functional PBE (R2 0.9996), verdict **NEAR-EXISTING** (goes to reviewer). check_derivs: OK.

## What is new
1. **One-orbital-interpolated exchange.** Spin-scaled exchange, per spin channel with n = 2 rho_s:
   e_x = -C_x n^{4/3} F_x(s, z), C_x = (3/4)(3/pi)^{1/3},
   F_x = F0(s) + z (F1(s) - F0(s)), z = tau_W/tau (per spin, so z = 1 exactly for one-orbital densities).
   F_k(s) = 1 + kappa_k - kappa_k / (1 + mu_k s^2 / kappa_k).
   - z -> 0 (many-orbital, slowly varying): F0 with kappa0 = 0.804 and mu0 = 0.2195 (PBE value; see candidates).
   - z = 1 (one-orbital): F1 with kappa1 = 0.804 and mu1 = 0.27583, solved from the **exact hydrogen-atom
     exchange energy, E_x(H) = -5/16 Ha** (one constant, one datum; solved at import with brentq, no dev-set data).
   - Linear-in-z blend is the only interpolation choice; it is not fitted.
2. Correlation: PBE correlation from libxc (`,PBE`) as a building block, unchanged.

Physical justifications:
- kappa = 0.804 saturates the Lieb-Oxford-type bound F_x <= 1.804 at large s (audit: 0 violations).
- The one-orbital form cannot use the gradient-expansion coefficient: at z = 1 the density is far from the GE
  regime, and the exact H-atom exchange requires a larger mu1 (mu = 10/81 with kappa = 0.804 gives E_x(H) = -0.293 Ha
  at best, so the exact -0.3125 Ha is unreachable with the GE value). mu1 = 0.2758 is therefore a derived constant,
  not a GE value.
- The exchange derivatives are exact (complex-step, no finite differences).

## Exact conditions
Satisfied: uniform-gas limit (F -> 1 as s -> 0, UEG rel. dev 3.5e-7); Lieb-Oxford bound on exchange (F_x <= 1.804);
non-positivity; exact spin-scaling of exchange; exact one-orbital exchange for the H atom (by construction, z = 1);
smoothness (derivatives OK).
Violated / not satisfied:
- Second-order gradient expansion of exchange for z -> 0 (mu0 = 0.2195, not 10/81). This is deliberate; see candidates.
- One-electron correlation: PBE correlation is retained, so the H atom total energy is -0.5065 Ha vs exact -0.5
  (self-correlation ~ -6.5 mEh). Candidate C addressed this but scored worse (see below).
- Derivative discontinuity, fractional charge/spin, strict-confinement, and noncovalent (dispersion) physics: none.
- Correlation GE and low-density limit are inherited from PBE.

## Candidates compared (dev-set, def2-svp; selection made by dev score among these three; all disclosed)
| id | z->0 mu0 | one-electron correlation gate | svp WTMAD | svp AE MAE | svp IP MAE | status |
|---|---|---|---|---|---|---|
| A | 10/81 (GE) | no | 3.50 | 8.17 | 9.61 | not selected |
| B | 0.2195 (PBE) | no | **3.08** | 7.33 | 7.50 | **selected** |
| C | 0.2195 (PBE) | yes, E_c x (1 - zeta^2 z^2) | 4.43 | 11.16 | 6.71 | not selected |
Logs: `score_svp_A_mu10_81_nogate.log`, `score_svp_B_mu0p2195_nogate.log`, `score_svp_C_mu0p2195_gate.log`.
Sources: `cand_A_mu10_81_nogate.py`, `cand_B_mu0p2195_nogate.py`, `cand_C_mu0p2195_gate.py`. The submitted file
equals candidate B (its gate code path is present but `GATE = False`).
The choice of mu0 (A vs B) is a choice between two named constants made by dev-set score; the reviewer should weigh that.
The gate (C) is physically motivated but worsened atomization energies more than it helped IPs.

Final QZ result for B: WTMAD 1.94 (AE 4.52, IP 5.42). Baselines: PBE 2.95, r2SCAN 1.28, B3LYP 0.97, wB97M-V 0.56.
So this submission beats PBE but not r2SCAN; it is not a new best.

## Known weaknesses (from the QZ per-item errors)
- X-H-containing molecules underbind: H2O -5.9, CH4 -8.0, NH3 -2.9, HF -4.1, H2 -4.6 kcal/mol.
- Second-row/heavier: N2 +11.4, F2 +9.0, HCN +6.9 (overbinding) — PBE-like behaviour.
- IP:B +11.1 kcal/mol (overestimated). IP:HF +0.2, IP:F +5.0.

## What I would try with more time
1. Replace the one-orbital form with one that matches the exact one-electron exchange hole pointwise (not just the H-atom
   energy), so the X-H underbinding from the one-orbital region is reduced without a dev-set fit.
2. Make the one-electron correlation gate weaker or smoothly switched so AEs are not penalized (candidate C showed
   the trade-off), or use a self-interaction-free correlation derived from the spin-resolved exact condition.
3. Add a derivation of mu0 from the exact gradient expansion of the spin-scaled exchange (10/81) together with an
   explicit correction for the low-density tail, and check whether it removes the PBE-like behaviour (candidate A).
4. Check second-row and barrier classes once a hidden-class proxy is available; no such data is used here.
