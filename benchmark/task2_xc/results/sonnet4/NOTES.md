# sonnet4: two-electron-singlet-consistent Stoll correlation in a 25% hybrid meta-GGA

submission.py: HYB = 0.25, META = True (tau used only through z_s = tau_W,s/tau_s). No RSH, no NLC (see candidates).

    E_xc = 1/4 E_x^HF + 3/4 E_x^PBE[n_a,n_b] + E_c
    E_c  = E_c^PBE[n_a,n_b] - sum_s z_s^2 E_c^PBE[n_s,0]                      (Stoll same-spin removal)
           + (c-1) (z_a z_b)^4 { E_c^PBE[n_a,n_b] - sum_s E_c^PBE[n_s,0] }      (NEW: opposite-spin completion)

## What is new
Stoll partition E_c = E_c^opp + sum_s E_c^ss. A spin channel holding one orbital (z_s = 1) has no same-spin correlation, so the
same-spin part is gated by z_s^2 (this first term is the frontier idea from sonnet3; I re-derived it, it is not new here).
Gating alone leaves PBE's *opposite-spin* part, which for PBE is right in He only by cancellation with its (wrong) same-spin part.
Removing the same-spin part therefore leaves He (and H2, and every two-electron singlet region) under-correlated.
New: the "opposite-spin correlation built for it". In regions where both spin channels are one-orbital (z_a z_b -> 1: a
two-electron singlet, e.g. He, H2 bonds, core 1s^2 pairs) the opposite-spin component is rescaled by c so that the **exact He
correlation energy** (-0.04204 Ha) is reproduced on the HF He density (aug-cc-pV5Z, computed by me, cnorm.py):
PBE opposite-spin part on that density = -0.030895 Ha, so c = 0.04204/0.030895 = 1.3607 (derived, not fitted to dev data).
Where either channel is not one-orbital (atoms, one-electron densities, uniform gas: z -> 0) the rescale switches off, so
- H atom: opposite part = E_PBE[n,0] - E_ss = 0 exactly, same-spin gated to 0, E_c = 0 exactly (H atom correlation-free);
- He: E_c exact by construction (2e singlet condition); H2 (dev) improves from -8.7 to -1.7 kcal/mol error in AE;
- uniform gas: z_s -> 0, E_c -> PBE = LDA/PW92 limit.
The switch (z_a z_b)^p uses p = 4 (candidate list below); p only controls how sharply the 2e-singlet regime is isolated.
Cap on z_s: C^1 cubic Hermite cap at 1 (identity below 0.8) so numerical z>1 does not make SCF noisy (a hard clip made H atom SCF not converge).
Exchange is plain PBE-form (kappa = 0.804 Lieb-Oxford, mu = beta pi^2/3 from the PBE gradient cancellation), scaled by 3/4.

## Exact conditions
Satisfied: uniform gas limit (audit dev 4e-7); spin scaling of exchange; Lieb-Oxford (audit: 0 violations); negativity (0);
E_c = 0 for every fully polarised one-orbital density; same-spin E_c = 0 for single-orbital channels;
exact He correlation energy on the HF-density (two-electron singlet normalization); high-density limit of PBE correlation.
Violated: exchange one-electron condition (PBE-form exchange, H -0.3 mEh/ He 1.2% error in Ex; the 25% HF only partly cancels);
10/81 gradient coefficient in exchange alone; no derivative discontinuity beyond 25% HF; no dispersion; opposite-spin correlation
low-density limit as PBE; the (c-1) rescale is calibrated on He only (He is the sole two-electron reference used).

## Audit
Verdict NEAR-EXISTING (closest PBE0, R2 0.99996; NEG 0, LO 0, UEG ok). Expected: audit samples unpolarized, mostly
one-or-many-orbital smooth densities where the gate/rescale is small; the new physics acts for z_s -> 1 regions (H, He, H2, cores).
Constants: kappa (LO), mu (PBE cancellation), 0.25 (1/n adiabatic connection), c (exact He Ec), none from the dev set.

## Candidates compared (all disclosed; dev WTMAD, def2-svp unless stated)
| candidate | AE MAE | IP MAE | WTMAD |
|---|---|---|---|
| PBE x + z^2 Stoll gate only (c = 1) | 4.99 | 7.15 | 2.21 |
| + opp-spin rescale, p=1 | 4.95 | 6.90 | 2.18 |
| + p=2 | 4.41 | 7.31 | 2.01 |
| + p=4 (SUBMITTED) | 4.34 | 7.23 | 1.98 |
| + p=8 | 4.47 | 6.89 | 2.00 |
| no gate, no rescale, two-norm exchange (below) | 18.09 | 13.76 | 7.34 |
| two-norm iso-orbital exchange + gate (z^2 weight) | 20.45 | 13.02 | 8.16 |
| two-norm exchange + gate, weight z^4 / z^8 / z^16 | 21.9 / 20.9 / 19.8 | 10.6 / 8.5 / 6.4 | 8.55 / 8.06 / 7.56 |
| submitted + VV10 NLC (b=6.0, C=0.0093, published constants) | 5.02 | 8.01 | 2.26 |
Selection among p in {1,2,4,8}, the exchange variants and NLC on/off used dev scores. No continuous parameter was fitted.
Two-norm exchange (rejected): F1 = (1+A s^2)/(1+B s^2.5), s^-1/2 tail, A=5.690, B=5.0207 solved from exact H-atom (-5/16) and
He (-1.02573 Ha) exchange simultaneously (norms.py); max F1 = 1.283 (LO ok). It meets both norms (the thing opus3 showed is
impossible for PBE-form) but ruins atomization energies (+20 kcal/mol overbinding): exact one-electron-region exchange
without matching correlation breaks the exchange-correlation error cancellation that molecules rely on. NLC rejected: it overbinds
molecules relative to atoms in the dev set, and no NCI data to justify it.

## Final dev score (def2-qzvp)
def2-qzvp: AE MAE 3.02, IP MAE 3.70 kcal/mol, WTMAD 1.30, MCM 3.36, all species converged (score.json). Baselines: r2SCAN 1.28, B3LYP 0.97, PBE0 ~1.3. Does NOT beat B3LYP or sonnet3 (1.00).

## Weaknesses / with more time
HF -11.6, H2O -5.9 (X-H, O/F underbinding: exchange SIE), CH4 +8.4, C2H2 +9.5, H2CO +9.7 (C overbinding), B IP +7.5. Not tested: second row,
barrier, NCI, PA. Next: a two-electron/one-electron exchange branch re-normed jointly with the correlation rescale; calibrate c from H2
exact (rather than He) as an independent check; range separation for barriers/IPs; a dispersion term with b derived rather than taken from literature.
Files: norms.py, cnorm.py (data generation for the derived constants), svp_*.log/json (ablations).
