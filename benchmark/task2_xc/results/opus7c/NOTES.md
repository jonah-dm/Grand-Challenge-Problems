# opus7c (CONSTRUCTION track, round 7): diagnosing opus6c's IP error, and two repairs to its one-orbital LC construction

Started 05:23 UTC. Phase 1 (05:23-06:08): diagnosis and derivation. No candidate SCF runs. The only SCF runs were UHF
reference densities. All candidate numbers in Phase 1 are non-self-consistent: the functional evaluated on UHF/cc-pVTZ densities
(work/hfeval.py). I call this the "HF-density proxy". Phase 2 is appended below.

## 0. What I build on and what I claim
Inherited from opus6c and claimed by no one here (organiser literature check, results/round6_literature_check.md):
- An LC hybrid with 100% LR-HF. ITYH attenuation. sonnet3's z-blended base exchange. opus4's pair-aware Stoll correlation.
- The hydrogenic one-orbital SR branch e = -(pi/2w^2) n_s^2 5Q(b) with weight zeta^2 z_s^2. Its R function is opus5c's.
- The finiteness => alpha_LR = 1 argument (opus6c, new per the check).
- L4, "fix omega so the branch respects a bound on the H atom". Fixing omega on H is itself known (Modrzejewski et al. 2016).
- The large-omega on-top expansion (Gill-Adamson-Pople 1996; Toulouse-Colonna-Savin 2004).
From KNOWN_CONDITIONS I use: exchange spin scaling, the one-electron limits, the Lieb-Oxford family of bounds, and the adiabatic-connection 1/4.
I accept the organiser's qualification (E6): the branch is finite pointwise but grows like n^(2/3), so it needs a cap, and the cap breaks exact lambda^2 scaling on compact shapes.

## 1. Diagnosis of the IP error (ROUND-7 FOCUS)
Method: UHF/cc-pVQZ densities (work/diag.py) and UHF/cc-pVTZ densities for the whole dev set (work/cache.py, work/hfeval.py).
The reference energy of each species is split as E_HF/QZ + E_c^CCSD(T)/CBS. Both pieces are stored in kit/dev_reference.json.
Any functional's error at the HF density then splits exactly into an exchange part, x = Delta(E_x^DFA - E_x^HF), and a correlation part,
c = Delta(E_c^DFA - E_c^ref). The proxy reproduces opus6c's SCF QZ results: WTMAD 1.255 against 1.21, IP(B) +9.7 against +9.9,
IP(F) +9.1 against +10.4, IP(HF) +2.7 against +4.5. So the IP error is not density-driven. At most 1-2 kcal/mol is.

**Table D1. opus6c IP error at the HF density (kcal/mol), cation minus neutral**
| item | total | exchange x | of which branch | correlation c |
|---|---|---|---|---|
| IP(B)  | +9.7 | +3.9  | +0.9 | +5.8  |
| IP(F)  | +9.1 | +21.3 | +1.3 | -12.2 |
| IP(HF) | +2.7 | +20.6 | +0.2 | -17.9 |
(QZ densities give the same split: B x +3.8, c +5.7; F x +22.1, c -12.6; HF x +22.0, c -18.4.)

Findings:
1. **The branch is not the cause.** The one-orbital branch shifts the IPs by at most +1.3 kcal/mol. The organiser suggested the compact-shape
   overestimate in the valence of B and F might matter. It accounts for +0.9 (B) and +1.3 (F) at most. In F+ the branch is in fact
   *less* binding than the base it replaces.
2. **For F and HF the error is in exchange, and it sits in the cation-neutral difference of the inherited SR semilocal base.** The SR-semilocal
   exchange loses about 21 kcal/mol too much exchange on ionization. The overly positive correlation loss (-12/-18) cancels only part of it.
   Neither channel alone explains it. The occupied alpha channel of F, unchanged in occupation, contributes +8.4 through relaxation, and the beta channel contributes +13.
   It is generic to semilocal exchange. Full-range PBE has x = +14.1 and B88 +13.3. HJS-attenuated PBE at the same omega has +17. So LC
   raises it by about 7 kcal/mol. The long-range semilocal part, which LR-HF replaces, used to offset it.
3. **Spatial region (exact conventional SR-HF energy density vs the semilocal one, radial bins about the nucleus, work/t_radial.py; gauge caveat).**
   F+ - F: the IP exchange error is +5.9 (r < 0.3 bohr), +14.8 (0.3-0.6), +7.2 (0.6-1.0) and -7.5 (r > 1.5).
   The exact SR exchange becomes more negative inside 0.6 bohr on ionization (-21.7). The 2p shell contracts onto the 1s/2s shell, so the
   core-valence exchange grows. The semilocal SR exchange barely responds (-1.0). So the error is core-valence exchange relaxation in
   the inner shell-boundary region, not the valence tail where the ionized electron lives. For B it is small and spread out (+1 to +2 per bin inside 1.5 bohr).
4. **For B the larger part is correlation (+5.8).** The DFA adds 10.4 kcal/mol of correlation for the 2p electron. The reference adds only 4.6.
   B+ is a Be-like 2s^2 ion with strong 2s/2p near-degeneracy. The semilocal correlation under-correlates it by much more than it under-correlates B.
   Every non-fitted functional shows the same: PBE +9.9, B3LYP +11.5, r2SCAN +12.1, LC-wPBE-like +10.6. No semilocal gate change addresses it.
5. **Levers from the focus list, quantified (HF-density proxy, all disclosed in section 5):**
   - SR-HF fraction a = 1/4 scales the F/HF exchange error by 3/4: IP MAE 7.2 -> 4.7-5.1. But at any of the derived omegas
     (0.5767, 0.4352, 0.4810, 0.5256) the AEs collapse (AE MAE 4.0-6.3; N2 -9 to -17). With this correlation the AE needs the full
     semilocal SR exchange in bonding regions, and the IP needs less of it in the core-valence region. A global a cannot do both.
     The construction needs region-selective exact exchange: a second, shorter range for core-valence exchange. The harness's single-omega RSH cannot express it.
   - The z_s weight of the branch: at most +-1.3 kcal/mol on the IPs. One x-c consistent alternative weight
     ("branch fraction = fraction of the s-electron correlation hole removed by the gate") leaks into closed shells, because opus4's gate is nonzero for 0 < z < 1. Proxy 1.38. Rejected.
   - The correlation gate for the cation: it removes 1.1 kcal/mol in B+ and 7.5 in F+. Gate changes can only move about 1 kcal/mol.
   - The range parameter: at fixed a, IP(F) is flat for omega in 0.3-0.6 (+9 to +10). The AE MAE is steep in omega.
   Conclusion: with 100% LR-HF and opus4 correlation, the IPs of F and HF cannot reach the B3LYP level through any of the listed levers
   without breaking the AEs. IP(B) is a static-correlation error that no non-fitted functional here fixes. The only target left is the
   AE side of WTMAD together with the derivation's soundness. Those are the two repairs below.

## 2. Repair 1: the base exchange violates exact spin scaling (known condition; a fix, not a novelty)
opus6c (and sonnet3) blend the base enhancement with z = tau_W/tau of the TOTAL density. Each spin channel's exchange then depends on
the other channel's density, which violates E_x[n_a, n_b] = (E_x[2n_a] + E_x[2n_b])/2. opus6c's NOTES list spin scaling as "satisfied". It is not.
Fix: use the channel indicator z_s = tau_W,s/tau_s for each channel's base exchange. Proxy: WTMAD 1.255 -> 1.187 (AE 2.37 -> 2.24, IP 7.18 -> 6.80).
Closed shells are unchanged (z_s = z when n_a = n_b).

## 3. Repair 2 (the construction step): replace the pointwise use of the integrated 1.174 bound by a pointwise bound, which fixes omega
opus6c's L4 caps the branch at F <= 1.174 and sets omega* by "the cap is just inactive at the H nucleus". The organiser (E6d) notes that 1.174
is an INTEGRATED bound for two-electron/one-orbital densities, used pointwise without proof. Here is a pointwise statement that is provable.

**Lemma B (bathtub bound on the one-orbital exchange energy density).** Let a spin channel hold one orbital, n_s = |phi|^2, int n_s = 1, with
m = sup n_s. In the conventional gauge e_x(r) = -(1/2) n_s(r) v_H[n_s](r), because the exact hole is -n_s(r'). Over all densities with
0 <= n <= m and unit charge, int n(r') |r-r'|^-1 dr' is maximized by the uniform ball of density m centred at r (bathtub principle; the
kernel is radially decreasing). So v_H[n_s](r) <= (3/2)(4 pi m/3)^(1/3). At the density maximum, n_s(r) = m and
    |e_x(r)| <= (3/4)(4 pi/3)^(1/3) m^(4/3),   i.e.   F_x = e_x/e_x^LDA(spin-polarized) <= (4 pi^2/18)^(1/3) = 1.29926.
For a fractional channel lambda n_s the same argument gives F <= 1.29926 lambda^(2/3). Because erfc <= 1, the same bound holds for SR exchange.
(The erfc-attenuated bathtub gives a tighter, omega-dependent value. See section 4, which explains why it cannot be used for the branch.)
Check: on the exact H atom at the nucleus, the conventional exchange has F = 0.788 <= 1.299.
**Corollary (lambda-consistency).** The uncapped branch scales at fixed shape as F proportional to lambda^(2/3) (e proportional to n^2). That is exactly the lambda
scaling of Lemma B's bound, which is attained by the same maximizer at every lambda. So the comparison "branch vs bound" made on the H atom at
lambda = 1 holds unchanged for every fractional hydrogenic channel. (The smooth cap form F(1+(F/F_cap)^8)^(-1/8), kept from opus6c, still trims
the branch by up to 8% where F is close to F_cap, as in opus6c. This is disclosed. It does not break lambda^2 on the H family, because the trimming scales with F/F_cap, which is
lambda-independent there only if F_cap also carried lambda^(2/3). It does not, so the trimming is lambda-dependent near the nucleus. That is a small residual violation.)
**L4 with Lemma B:** omega* solves (5/9)(6 pi)^(2/3) R(1/omega^2) = (4 pi^2/18)^(1/3), giving **omega* = 0.525607**. The cap is F <= 1.29926.
Novelty status: the bathtub principle is textbook (Lieb-Loss). Using it as a pointwise one-orbital exchange bound, the lambda^(2/3)
matching to the amplitude-quadratic branch, and the omega it gives are, as far as I know, not on the known list. I claim only these, and only modestly.

## 4. A caveat found on the way (negative result; it applies to opus6c's L4 as well)
The branch is not the conventional-gauge energy density. It integrates exactly for hydrogenic channels, but it puts more energy at the
nucleus than the conventional density does: at the H nucleus with omega = 0.5767, F_branch = 1.174 against F_conv^SR = 0.380.
Applying the erfc-attenuated bathtub bound (rigorous for conventional SR energy densities, work/bound.py) would require omega* = 1.14.
That is almost HF. So every pointwise-bound rule for omega, opus6c's 1.174 included, mixes gauges. The cap and the omega rule are a
pointwise regularization of a gauge-free model, and they are justified only as that. I state this openly. The full-range Lemma B bound is the
weakest rigorous pointwise statement of the right form (density-only, attained at the maximum, where the branch attains its supremum).
The value used for the cap is a choice between two derived bounds: 1.174 (integrated, opus6c) and 1.29926 (pointwise, Lemma B).
I disclose it as a discrete candidate choice (section 5).

## 5. Candidates (pre-registered before any SCF run; HF-density proxy shown; SCF results appended in Phase 2)
| id | base z | a (SR-HF) | omega | cap | rule | proxy WTMAD (AE/IP) |
|---|---|---|---|---|---|---|
| S0 = opus6c K2b | total | 0 | 0.5767 | 1.174 | L4 + 1.174 | 1.255 (2.37/7.18) |
| S1 | channel | 0 | 0.5767 | 1.174 | Repair 1 | 1.187 (2.24/6.80) |
| S2 | channel | 1/4 | 0.5767 | 1.174 | + adiabatic-connection SR fraction | 2.543 (6.25/4.88) |
| S3 | channel | 1/4 | 0.4352 | 1.565 | bound on the (1-a) semilocal share | 1.928 (4.60/4.63) |
| S4 | channel | 1/4 | 0.4810 | 1.423 | bound on the hybrid-averaged total (conv. SR-HF at nucleus) | 2.046 (4.91/4.71) |
| S5 | channel | 0 | 0.5256 | 1.299 | Repair 1 + Lemma B | 0.954 (1.63/6.61) |
| S6 | channel | 1/4 | 0.5256 | 1.299 | Lemma B + 1/4 | 2.267 (5.51/4.79) |
| S5t | total | 0 | 0.5256 | 1.299 | Lemma B only | 1.074 (1.89/7.01) |
Diagnostic-only evaluations, not candidates, also disclosed: omega scans at a = 0 and 1/4 with cap 1.174 (work/t_w.py); HJS-PBE/PBE-c
(LC-wPBE-like) at omega = 0.3/0.4/0.577; full-range PBE/B88/r2SCAN exchange with PBE correlation; channel-z without branch;
the "g" weight (rejected); S2' (a = 1/4, omega = 0.4352, cap 1.174, an incoherent mix, proxy 1.054, not a candidate because its cap and omega come from different rules).
Selection rule (fixed now): among S1, S5 (a = 0, which keeps opus6c's delocalization fix intact), take the lower SVP dev WTMAD, subject to
all species converging, the SCF H atom at QZ within 1.5 kcal/mol, H2+ at 6 bohr within 1 kcal/mol, and the H+H2 barrier within 1.5 kcal/mol of 9.6.
The a = 1/4 candidates are excluded on the proxy AE collapse. They are listed because the focus asked for the lever.

### Phase-1 checks of Lemma B and of the implementation (no SCF)
- Conventional one-orbital exchange at the density maximum, F = e_x/e_x^LDA,pol: H atom 0.787; normalized Gaussian 1.075 (any exponent);
  uniform ball 1.29926 (the bound is attained). The bound holds only at the maximum. Elsewhere the exact F is unbounded (the s^-2 law).
  The branch attains its supremum at critical points (b -> 0), which is where the bound applies.
- Fractional one-electron tests at exact densities (work/frac.log):
  S5: E_xc[H] = -0.31234 (exact -0.3125), f(1/2) = 0.2508 (exact 0.25), f(1/4) = 0.0629, H2+ dissociation limit at the exact density -0.31 kcal/mol;
  Gaussian lambda = 1 ratio 1.093, Gaussian amplitude ratio 4E[n/2]/E[n] = 1.108 (cap, as disclosed). S0/S1 give 0.2507 and -0.28.
- The submission (O7_K = S5) reproduces the proxy model to machine precision. check_derivs: DERIVATIVES OK (both spin paths; check_derivs.log).
  Audit: DISTINCT (closest single CAM-B3LYP R^2 = 0.908; mix R^2 = 0.967; NEG 0; LO 0). The UEG flag (0.72) is the audit's global-hybrid
  assumption, as for opus6c: eval_xc returns only the SR semilocal part.
- Second-row transfer check. It was not used to set any constant. The reference is my own CCSD(T)/cc-pVTZ (fc), so every error carries
  the same positive basis bias, and only differences between candidates mean anything (work/row2.py, work/t_row2.py; HF-density proxy):
  S0: HCl +2.0, H2S +4.0, PH3 +5.3, SiH4 +7.0, Cl2 +4.0, IP(Cl) +11.0, IP(S) +10.6
  S1: HCl +1.9, H2S +3.8, PH3 +5.1, SiH4 +6.9, Cl2 +3.8, IP(Cl) +10.7, IP(S) +10.4
  S5: HCl +1.9, H2S +3.4, PH3 +4.0, SiH4 +4.0, Cl2 +5.3, IP(Cl) +10.9, IP(S) +10.7
  No candidate degrades second-row AEs. The IP(Cl), IP(S) pattern (x about +6, c about +5) is the same exchange/correlation split as in the first row.
## Phase 2 started 06:08 UTC
