# Known exact conditions and known constructions (no credit for re-deriving these)

A submission earns novelty credit only for a condition or construction NOT on this list. Using
items on this list is fine and expected; presenting them as new is not. State explicitly which
listed items you use and what you add beyond them.

## Exact conditions on E_xc already in the literature
- Uniform electron gas limit (LDA exact for constant density); spin-scaling of exchange
  E_x[n_a,n_b] = (E_x[2n_a] + E_x[2n_b])/2; uniform coordinate scaling of exchange (E_x ~ lambda);
  scaling inequalities for correlation (Levy–Perdew).
- Second-order gradient expansion of exchange (mu = 10/81) and correlation (beta = 0.0667);
  fourth-order exchange GE (Svendsen–von Barth).
- Lieb–Oxford bound: E_xc >= -C_LO ∫ n^{4/3}, with the tight-ish constant for exchange-only and
  the two-electron bound F_x <= 1.174 (used by SCAN; kappa_0 = 0.174).
- Exchange hole normalization (−1), correlation hole normalization (0), negativity of the exchange
  hole and of e_xc.
- One-electron limits: E_x[n_1] = −J[n_1], E_c[n_1] = 0 (self-interaction freedom; Perdew–Zunger
  1981). Fully spin-polarized one-orbital densities have no correlation; a spin channel holding one
  orbital has no same-spin correlation (Stoll partition; Becke B95, 1996; TPSS correlation, 2003).
- Iso-orbital indicators: z = tau_W/tau (TPSS), alpha = (tau − tau_W)/tau_UEG (SCAN), and their
  use to detect one-orbital, slowly varying and overlap regions.
- High-density (strong-correlation-free) limit of correlation; low-density / strictly correlated
  limit (Seidl, Perdew, Levy; SCAN's beta(rs)); strongly interacting limit of the adiabatic
  connection.
- Asymptotic −1/r behaviour of v_xc; e_x tail decaying as s^{−1/2} in one-orbital regions;
  exact exchange in the high-density limit.
- Piecewise linearity in fractional electron number and the derivative discontinuity (Perdew,
  Parr, Levy, Balduz 1982); flat-plane condition in fractional spin (Mori-Sánchez, Cohen, Yang
  2009); constancy condition for fractional spins.
- Koopmans-type conditions on orbital energies (IP theorem; generalized KS gap = fundamental gap
  for the exact GKS, Perdew et al. 2017).
- Appropriate norms: exact energies of the H atom, He atom, hydrogenic ions, the Hooke's atom,
  jellium surfaces, large-Z asymptotics of atoms (Elliott–Burke).
- Adiabatic-connection argument for the global hybrid fraction a = 1/n (PBE0's 1/4).
- Size consistency; invariance under unitary rotations of occupied orbitals; translational and
  rotational invariance; correct behaviour under density scaling of the kinetic energy (von
  Weizsäcker bound tau >= tau_W).

## Known constructions
- GGA forms: B88, PBE/PBEsol/revPBE/RPBE, PW91, B97 power series.
- Meta-GGA forms: TPSS, revTPSS, SCAN/rSCAN/r2SCAN (alpha interpolation + norms), M06-L/MN15-L
  (fitted), B97M-V (combinatorial fit), TM, MVS.
- Correlation: PW92, LYP, PBE, B95 (Stoll partition + tau), TPSS one-electron SIE-free gate.
- Hybrids: PBE0, B3LYP, SCAN0; range-separated: HSE, CAM-B3LYP, wB97X/wB97M; local hybrids;
  double hybrids (B2PLYP, XYG3).
- Nonlocal correlation: vdW-DF, VV10; dispersion corrections D3/D4.
- Self-interaction corrections: Perdew–Zunger SIC, FLOSIC, LSIC; Koopmans-compliant functionals.
- ML functionals: NeuralXC, DM21 (fractional-charge/spin trained), Skala, Kohn–Sham regularizer.

## What would count as new (examples, not a checklist)
- A proven inequality or sum rule not above (e.g. a tighter bound for a density class; a relation
  between exchange and correlation hole moments; a constraint from the strictly-correlated limit
  on gradient-dependent terms; a bound relating tau, grad n and E_xc beyond the von Weizsäcker
  inequality).
- A construction that satisfies a known condition that NO existing functional satisfies (e.g.
  the flat-plane condition in a semilocal or local-hybrid form; exact one-electron exchange
  pointwise for all one-orbital densities, not just the H norm).
- A new ingredient (beyond n, grad n, tau, exact exchange) with a derived, not fitted, role.
