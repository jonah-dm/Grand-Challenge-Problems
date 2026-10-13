# Round 4 — prior art from rounds 2–3, RSH/VV10 enabled, 90 min

| Model | Submission | Audit | Dev WTMAD | AE MAE | IP MAE |
|---|---|---|---|---|---|
| haiku4 | sonnet3's Stoll gate + haiku3's H-normed exchange (recombination) | NEAR-EXISTING (PBE0) | 1.30 | 3.01 | 3.71 |
| sonnet4 | sonnet3 + opposite-spin rescale in two-electron-singlet regions, constant from exact He correlation | NEAR-EXISTING (PBE0) | 1.30 | 3.02 | 3.70 |
| opus4 | sonnet3 exchange; correlation gate made pair-aware: same-spin removal only where the region is not an opposite-spin pair, E_c = E_c^PBE − Σ_s z_s²[1 − z_s'²(1−ζ²)]E_c^PBE[n_s,0]; no new constants | DISTINCT (PBE0, R²=0.97) | **1.08** | 2.39 | 3.89 |

Best to date remains sonnet3 (1.00). Baselines: B3LYP 0.97, ωB97M-V 0.56.

## What happened
- **Diagnosis-driven work appeared.** Opus computed exact-density correlation energies (He −0.0420,
  H₂ −0.0409 Ha) and showed sonnet3's gate gives −0.0313/−0.0307: it removes opposite-spin
  correlation it has no business removing in two-electron singlets. That 0.010 Ha is sonnet3's
  H₂ (−6.0) and LiH (−5.8 kcal/mol) errors. The fix is clean and constant-free. It beats
  the error it targeted (H₂ → +1.6, LiH → −0.4) and loses elsewhere (C–H overbinding).
- **Sonnet attacked the same defect** with an opposite-spin rescale normed to He, and regressed
  overall (1.00 → 1.30). It also found a two-norm (H+He) exchange form with the s^-1/2 tail,
  contradicting opus3's impossibility claim for the PBE shape (different shape), but it overbinds
  by ~20 kcal/mol and was rejected.
- **Haiku recombined** two prior routes honestly; no gain.
- **Neither RSH nor VV10 was used in a submission.** Sonnet tried VV10 (worse on the dev set,
  which has no noncovalent items, as expected); Opus tried to fix the range parameter from the
  IP theorem on H and He and found no transferable value. Both correctly refused to tune.
- **Nothing new relative to the literature.** The pair-aware gate is close in spirit to B95/M06
  correlation and TPSS's treatment; see KNOWN_CONDITIONS.md. The round's value is the quality
  of the diagnosis, not novelty.

## Round 5 changes
Publication-worthy bar; THEORY and CONSTRUCTION tracks; 180-min two-phase budget with a
derivation-only first hour; KNOWN_CONDITIONS.md with explicit "no credit" list; Haiku dropped.
