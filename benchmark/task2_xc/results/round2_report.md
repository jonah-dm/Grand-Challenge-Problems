# Round 2 — novelty rules

Rules: new constraints / new forms / ML-with-constraints; wrappers, mixes and re-fits disqualified;
audit.py given to contestants; 60 min, 1 thread; HYB bug fixed.

| Model | Submission | Audit | Dev WTMAD (QZ) | AE MAE | IP MAE |
|---|---|---|---|---|---|
| haiku2 | meta-GGA: PBE-type exchange × one-orbital factor G(z)=1−z/(4(1+z)), PW92/PBE correlation | DISTINCT (closest SCAN0, R²=0.88) | — (7/24 species unconverged) | — | — |
| sonnet2 | "OSC-PBE0": exchange blend z²-weighted between GE2-exact PBE form (μ=10/81) and an s^-1/2-tail form with κ from Lieb–Oxford; correlation × [1−ζ²z²]; 25% HF | DISTINCT (closest PBE0, R²=0.97) | **1.35** | 3.21 | **3.32** |
| opus2 | "BIG" exchange: linear interpolation in β=α/(1+α) between iso-orbital and UEG limits, GE2 enforced on the integrated energy (μ₁ derived), κ₁ from LO on the α→∞ branch; r²SCAN correlation; no HF | DISTINCT (R²=0.99) | 3.23 | 7.51 | 9.12 |

Baselines: PBE 2.95 (AE 7.3 / IP 5.2), r²SCAN 1.28 (2.6 / 6.3), B3LYP 0.97 (1.7 / 6.4), ωB97M-V 0.56 (1.3 / 1.3).
Round-1 Sonnet mix (MN15/M06-2X/PW6B95), rescored with the fixed harness: 1.62 (4.2 / 1.5).

## What changed from round 1
- **Zero functional shopping.** All three built new forms from exact constraints and derived their
  constants (gradient expansion 10/81, Lieb–Oxford 1.804, two-electron bound 1.174, H-atom and He norms).
  The rule change plus the audit gate was sufficient; no agent tried to game the audit.
- **Sonnet's self-correlation gate works as designed.** Making E_c vanish for one-electron densities
  cuts the IP error from ~6 (B3LYP, r²SCAN) to 3.3 kcal/mol. The cost is in atomization energies.
- **Opus's constraint-only exchange underbinds X–H bonds by 11–17 kcal/mol.** It reports that the
  variant closest to satisfying all constraints was flagged NEAR-EXISTING to r²SCAN and withheld it.
  Lesson: the constraints alone nearly rebuild r²SCAN; the audit should not punish that.
- **Haiku failed on numerics, not novelty.** Closed-shell derivatives correct; spin-polarized ones
  wrong by orders of magnitude. check_derivs.py (new) catches this in seconds.
- **Disclosure behaviour was good.** Sonnet and Opus both volunteered that they selected between
  structural variants after seeing dev-set energies, without being asked.

## Rule / harness changes for round 3
1. NEAR-EXISTING goes to reviewer judgment rather than automatic disqualification.
2. Selection among a few disclosed candidate forms is allowed; undisclosed selection disqualifies.
3. check_derivs.py required before any SCF run.
4. score.py: damping/level-shift fallback before Newton; Newton failure no longer aborts.
5. Round-2 NOTES.md files are given to round-3 contestants as prior art, so rounds accumulate.
6. Hidden-set B class list (not molecules) disclosed so contestants design for generality.
7. Budget 90 min.
