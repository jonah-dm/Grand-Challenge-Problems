# opus3: "OEC-h": one-electron-consistent hybrid meta-GGA

`submission.py`: global hybrid meta-GGA, HYB = 0.25, META = True.

    E_xc = 1/4 E_x^HF + 3/4 E_x^SL[n_a, n_b, tau_a, tau_b] + E_c^G

## What is new
Prior rounds imposed the one-electron (self-interaction) condition on only one of the two
semilocal pieces. opus2 normed the *exchange* to the H atom and left correlation alone. sonnet2
gated the *correlation* to zero for one-orbital densities and left exchange un-normed, and its
atomization energies got worse. In a global hybrid the one-electron condition is a condition on
the **sum**: for any one-electron density, (1-a) E_x^SL + E_c must equal (1-a)(-J).
If you fix one piece without the other, you break the error cancellation PBE0 relies on: PBE0's
H atom is right only because the -5 mEh exchange error cancels the +6 mEh self-correlation.
Here both conditions are imposed together, and every constant is derived:

1. **Iso-orbital exchange with a derived H-atom norm.** Fx(s,z) = 1 + kappa - kappa/(1 + mu(z) s^2/kappa),
   with z = tau_W/tau evaluated per spin channel (so exchange obeys exact spin scaling),
   mu(z) = mu0 + (mu1 - mu0) z^2.
   - mu0 = beta*pi^2/3 = 0.21951. In the slowly varying (z -> 0) limit this cancels the
     second-order gradient term of the PBE correlation that it is paired with. That is the
     exchange-correlation gradient-cancellation argument, not a fitted value.
   - mu1 = 0.275832 is solved numerically at import time (`_hydrogen_mu1`) so that
     E_x^SL[H atom] = -5/16 Ha exactly. Because z = 1 exactly for any one-orbital density,
     the norm holds for the H atom whatever w(z) is.
   - kappa = 0.804 (Lieb-Oxford, as in PBE). I did not raise kappa: normalizing through kappa
     would need kappa = 1.525, which breaks the spin-scaled LO bound.
2. **Correlation that vanishes for every fully spin-polarized one-orbital density:**
   E_c = E_c^PBE x G, with G = 1 - zeta^4 z^2 (z = tau_W/tau of the total density). For any
   one-electron density zeta = 1 and z = 1, so G = 0. For spin-unpolarized densities and the
   uniform gas G = 1.
   Together with item 1 and a = 1/4, the H atom's xc energy is exact: SCF E(H) = -0.49994 Ha,
   compared with -0.50104 for PBE0 at the same basis.
3. **Validation on data I generated myself (not used for any parameter):** the He HF density
   (aug-cc-pV5Z, `he_check.py`) gives E_x = -1.0362 Ha with the H-normed z = 1 form, against
   -1.0257 exact and -1.0136 for PBE. The error is +1.0%, compared with -1.2% for PBE. I tried to
   meet both the H and He norms with the two constants (kappa1, mu1) of the PBE form
   (`two_norm.py`). No solution exists (kappa1 -> infinity), so a single PBE-shaped
   iso-orbital enhancement cannot be exact for both one- and two-electron densities. That is a
   real limitation of the form, recorded under "Violated".

Derivatives: exchange and gate derivatives use complex-step differentiation, which is exact to
machine precision. PBE correlation and its derivatives come from libxc and are combined with
the gate by the product rule. `check_derivs.py` passes for both the unpolarized and the
spin-polarized paths.

## Exact conditions
Satisfied:
- uniform-gas limit (audit deviation 4e-7);
- exact spin scaling of exchange;
- uniform coordinate scaling of exchange;
- Lieb-Oxford: Fx <= 1.804, and the audit finds no pointwise LO violations;
- e_xc <= 0;
- exact H-atom exchange (semilocal part) and zero H-atom correlation, so the hybrid xc
  energy is exact for the H atom;
- E_c = 0 for every fully polarized one-orbital density;
- the high-density correlation limit and the GE cancellation inherited from PBE.

Violated or not imposed:
- exchange is exact only in the integrated sense for the H atom, not pointwise, and not for
  general one-electron densities;
- the He exchange norm is missed by 1%;
- the second-order gradient expansion of exchange alone (mu0 is the PBE value, not 10/81);
- no derivative discontinuity or fractional-charge linearity beyond what 25% HF gives;
- no dispersion, so noncovalent interactions rely on PBE0-like behaviour;
- the z^2 and zeta^4 shapes are interpolation choices among the candidates listed below.

## Audit
Verdict NEAR-EXISTING (closest PBE0, R^2 = 0.99965; mix R^2 0.99996). This is expected. The
audit samples only spin-unpolarized densities, where G = 1 and correlation is exactly PBE, and the
z-dependence of mu changes exchange only mildly. The derivation above sets each constant from
an exact condition: mu0 from GE cancellation, mu1 from the H norm, kappa from LO, a = 1/4 from
the PBE0 adiabatic-connection argument. The new physics, a jointly one-electron-exact
correlation and exchange, acts on spin-polarized species, which the audit does not sample.
NEG 0, LO 0, UEG ok.

## Candidates compared (all disclosed; dev-set WTMAD)
At def2-TZVP (AE MAE / IP MAE / WTMAD):
| candidate | exchange w(z) | gate | norm | TZ WTMAD | AE | IP |
|---|---|---|---|---|---|---|
| PBE0 (reference, libxc) | - | none | - | 1.61 | 3.68 | 4.93 |
| nogate | z^2 | none | H | 3.40 | 8.47 | 5.84 |
| w1 | z | 1-zeta^2 z^2 | H | 2.31 | 5.66 | 4.55 |
| w2 | z^2 | 1-zeta^2 z^2 | H | 1.98 | 4.83 | 4.11 |
| he | z^2 | 1-zeta^2 z^2 | He (mu1 = 0.2490) | 2.49 | 6.22 | 4.08 |
| g4 | z^2 | 1-zeta^4 z^2 | H | 1.67 | 3.79 | 5.31 |
| w4 | z^4 | 1-zeta^2 z^2 | H | 1.71 | 4.06 | 4.20 |
| w4g4 | z^4 | 1-zeta^4 z^2 | H | 1.82 | 4.13 | 5.76 |
At def2-QZVP (final basis):
| candidate | QZ WTMAD | AE MAE | IP MAE |
|---|---|---|---|
| PBE0 (reference) | 1.31 | 2.91 | 4.66 |
| g4 (**submitted**) | **1.59** | 3.63 | 4.87 |
| w4 | 1.64 | 3.94 | 3.76 |
| w4g4 | 1.58 | 3.53 | 5.31 |
g4 and w4g4 tie (difference 0.01). I submitted g4 because w = z^2 is the lower-order and
smoother choice, which is better for SCF stability.
No continuous parameter was adjusted. The only choices made after seeing dev scores were the
powers in w(z) (1, 2, 4) and in the gate (zeta^2 vs zeta^4), and H vs He norm, all listed above.

## Diagnosis (what the dev set taught)
- Each one-electron condition on its own shifts atoms a lot relative to molecules. The H norm on
  its own lowers the C atom by 37 kcal/mol and CH4 by only 41, so CH4 underbinds by 7.7 kcal/mol
  more. The correlation gate on its own raises H by 3.7 and N by 6.3 kcal/mol. Imposed together,
  the two largely cancel, which is the point of the construction. The residual is overbinding of
  C-H and underbinding of O/F-containing molecules.
- The gate with zeta^2 also suppresses correlation in partially polarized open-shell valence
  regions (C, N atoms), where genuine opposite-spin correlation exists. That overbinds N2 and
  C2H2. zeta^4 restricts the removal to near-full polarization. It fixes this but loses part
  of the IP gain, because the B 2p tail is only partly polarized.

## What I would try with more time
- An iso-orbital exchange shape that meets both H and He norms, for example with s-dependent
  saturation that decays as s^-1/2.
- Replace the zeta-power gate with a per-spin-channel Stoll-type self-correlation removal that
  touches only same-spin correlation and keeps opposite-spin correlation of open-shell atoms
  intact. This needs an opposite-spin correlation built for it, because PBE's He correlation is
  right only by cancellation.
- Range separation (the harness supports only global HYB) for barrier heights and IPs, and a
  nonlocal VV10-type term for noncovalent interactions.

## Final dev-set score (def2-QZVP, score.json)
**WTMAD = 1.59** (AE MAE 3.63 kcal/mol, IP MAE 4.87 kcal/mol). All species converged.
This does **not** beat the baselines: PBE 2.95 (beaten), r2SCAN 1.28, B3LYP 0.97, wB97M-V 0.56.
It is also worse than plain PBE0 at QZ (1.31). Satisfying the joint one-electron condition
exactly makes the H atom exact. It does not improve AEs on its own, because the residual
atom-versus-molecule imbalance (C-H overbinding, O/F underbinding) is set by the many-electron
behaviour of the iso-orbital exchange, which a one-electron norm does not constrain.
