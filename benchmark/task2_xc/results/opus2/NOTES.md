# opus2 submission: "BIG" exchange (bounded-indicator, integrated-gradient-expansion meta-GGA)

`submission.py`: pure meta-GGA (HYB = 0, META = True).

    E_xc = E_x^BIG[rho, sigma, tau] + E_c^r2SCAN      (r2SCAN correlation taken from libxc as a building block)

## What is new (exchange only; correlation is unmodified r2SCAN)

    Fx(p, alpha) = [ F0 + (F1(p) - F0) * W ] * g(s)
    F0    = 1 + kappa0,  kappa0 = 0.174                 (alpha = 0, the iso-orbital limit)
    F1(p) = 1 + kappa1 - kappa1/(1 + mu1 p/kappa1)
    W     = 2*beta,   beta = ah/(1+ah) in [0, 1),   ah = alpha^2 (1+eps)/(alpha+eps)
    alpha = (tau - tau_W)/(tau_unif + eta*tau_W)
    g(s)  = 1 - exp(-a1/sqrt(s))

1. **The interpolation is linear in a bounded indicator.** SCAN-type functionals switch between
   the iso-orbital limit (alpha = 0) and the slowly varying limit (alpha = 1) with a fitted
   switching function, plus a separately fitted alpha -> infinity branch (c1x, c2x, dx in SCAN).
   Here Fx is *linear* in beta = alpha/(1+alpha) (smoothed near 0), which maps alpha in [0, inf)
   to [0, 1). Fixing the two physical end points (beta=0: iso-orbital; beta=1/2: uniform gas)
   fixes the entire alpha dependence, including the alpha -> infinity limit
   Fx -> (2 F1 - F0) g, with no extra parameter. Fx decreases monotonically with alpha, as
   required for the correct ordering of single-orbital, metallic and weak-bond regions.
2. **The gradient expansion is enforced at the energy level, not pointwise.** Because W has a
   nonzero slope at alpha = 1, Fx contains a term linear in (alpha - 1). That term is second order
   in gradients (alpha - 1 = -(40/27)p + (20/9)q). After integrating by parts
   (int rho^{4/3} q = (1/3) int rho^{4/3} p), it shifts the effective coefficient of p by
   (20/27 + 5 eta/3) kappa0 dW/dalpha|_1. mu1 is then *derived* so that the total second-order
   coefficient of the exchange energy is exactly 10/81. This gives mu1 = 0.0408, which is not a
   chosen value. Existing meta-GGAs instead make the switch flat at alpha = 1 so that the
   pointwise expansion is untouched. Using the energy-level form of the exact condition is the
   construction-level novelty.
3. **kappa1 is derived from Lieb-Oxford through the linear extrapolation.** The alpha -> infinity
   limit 2 F1 - F0 must obey Fx <= 1.804 for all s. That gives kappa1 <= (0.804 + 0.174)/2 = 0.489,
   and the bound is saturated. PBE's 0.804 and SCAN's fitted 0.065 are both replaced by this
   derived value.

## Parameters and where they come from (none fitted to the dev set or to any reference data)
| constant | value | origin |
|---|---|---|
| mu (GE) | 10/81 | second-order gradient expansion of exchange |
| kappa0 | 0.174 | tight bound Fx <= 1.174 for two-electron / iso-orbital densities |
| kappa1 | 0.489 | Lieb-Oxford applied to the alpha -> infinity limit (item 3) |
| mu1 | 0.0408 | integrated GE = 10/81 (item 2) |
| a1 | 4.9487 | exact exchange energy of the H atom, -5/16 Ha (solved numerically, `_hydrogen_norm`) |
| eta | 0.1 | regularization of alpha in density tails (as in r2SCAN, eta=1e-3 there); value set by SCF convergence on Li/LiH, not by accuracy; its effect on the GE is included in mu1 |
| eps | 0.05 | removes the kink of alpha-clipping at alpha=0 (SCF stability on H atom); its effect on the slope is included in mu1 |

Validation of the H norm on a second system (data I generated myself: He HF/aug-cc-pV5Z
density, `he_density.py`): the alpha = 0 limit gives E_x(He) within 4.8 mHa (0.5%) of exact
HF exchange (-1.0257 Ha). No parameter was fitted to it.

## Exact conditions
Satisfied (exchange): uniform-gas limit (Fx(0,1) = 1; audit UEG deviation 3.5e-7); correct spin
scaling (exact); uniform coordinate scaling (Fx depends only on s and alpha); second-order GE
coefficient 10/81 at the integrated level; Lieb-Oxford Fx <= 1.804 everywhere (checked
analytically, including the alpha -> infinity branch); tight two-electron bound Fx <= 1.174 at alpha = 0;
Fx ~ s^{-1/2} as s -> infinity (exchange energy per particle vanishes under nonuniform scaling); exact
hydrogen-atom exchange energy; negativity of e_x. Correlation (r2SCAN) inherits r2SCAN's
conditions: one-electron self-correlation free, UEG, GE, and so on.
Audit: NEG 0, LO violations 0, UEG ok, verdict **DISTINCT** (mix R2 = 0.9981).

Violated / not satisfied: the *pointwise* second-order GE is not reproduced, only its integral.
The fourth-order GE is not imposed. Exchange is not exact for every one-electron density, only
normed to H, so one-electron self-interaction remains. The functional is semilocal, so it has no
derivative discontinuity, wrong -1/r asymptotics of v_x, and fractional-charge delocalization
error (which hurts IPs). There is no dispersion. Exchange and correlation were not co-designed:
r2SCAN correlation was built to pair with SCAN exchange.

## Dev-set results
SVP (quick checks; SVP baselines I ran: B3LYP 1.46, PBE0 2.12, r2SCAN 2.11, r2SCAN0 4.18):
final functional WTMAD 3.81 (AE MAE 8.5, IP MAE 13.1) at def2-SVP.
def2-QZVP (final, `score.json`): see bottom of this file.

The functional does **not** beat the baselines on the dev set. Its systematic error is
underbinding of X-H bonds (NH3, H2O, HF) and too-small IPs for F and HF. My diagnosis: W rises
almost linearly from alpha = 0, so bonded H regions (nearly single-orbital, alpha ~ 0.2-0.5)
quickly lose the H-normed iso-orbital enhancement.

## Design history (honest record; all on def2-SVP)
- v3: log-symmetric W(alpha) = W(1/alpha) with 25% HF: WTMAD 4.10, overbinding hydrocarbons. Rejected
  on physics: it violates the monotone decrease of Fx with alpha.
- v5/v8: beta-linear W with 25% HF: 5.3. Hybridizing SCAN-type exchange underbinds, as r2SCAN0 does,
  so HYB = 0 was adopted.
- v10/v11: cubic W = 1 + t^3, flat at alpha = 1, pointwise GE: 4.7-4.9.
- v12 (submitted): beta-linear + integrated GE + LO-derived kappa1: 3.81.
- Tried: pushing the whole GE into the alpha slope (mu1 = 0, W = 2 alpha^n/(1+alpha^n),
  n = 1.564 derived). The audit flags it NEAR-EXISTING (closest r2SCAN, R2 = 0.9989). The
  constraints alone essentially rebuild r2SCAN, so it was not submitted.
The dev set was used to catch SCF failures and to reject forms. That is a form of selection
among a handful of physically motivated forms, which I am disclosing; no continuous parameter
was tuned.

## What I would do with more time
- Derive the alpha-shape near 0 from a model (for example the exact exchange-energy density of
  bonded H2 / He2 at stretched geometries) instead of choosing it. Then add a second
  iso-orbital norm (He exchange density, not just the energy) to fix how quickly the H-like
  enhancement decays.
- Co-design the correlation: impose the same bounded-indicator linear interpolation on
  correlation, with the one-electron (alpha = 0, Ec = 0) and the alpha -> infinity limits derived together,
  so that exchange-correlation error cancellation is built in.
- Provide analytic fxc so PySCF's Newton fallback works, and reduce eta.
- Range-separated or local-hybrid exact exchange to fix IPs (delocalization error). That needs
  score-harness support for RSH, which the current harness lacks.

## Final def2-QZVP dev-set score (score.json)
**WTMAD = 3.23** (AE MAE 7.51 kcal/mol, IP MAE 9.12 kcal/mol). All 27 species converged.
For reference, PBE scores 2.95, r2SCAN 1.28, B3LYP 0.97 and wB97M-V 0.56, so this functional
does not beat any of the baselines on the dev set.
