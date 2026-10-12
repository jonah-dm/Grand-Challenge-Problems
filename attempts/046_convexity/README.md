# Problem 046 (Physics Open Problems catalogue): convexity in particle number of atomic ground-state energies

**Problem.** Fix a nuclear charge Z > 0 and let E(Z,N) be the ground-state energy (bottom of the
spectrum) of the nonrelativistic N-electron atomic Hamiltonian with a fixed nucleus,

    H(Z,N) = sum_i ( -Delta_i/2 - Z/|x_i| ) + sum_{i<j} 1/|x_i - x_j| ,

acting on antisymmetric functions of N space-spin variables. Prove that N -> E(Z,N) is convex,

    E(Z,N-1) + E(Z,N+1) >= 2 E(Z,N)   for every N >= 1,

or exhibit a pair (Z,N) violating it. (Catalogue entry 046; status "open", difficulty 5/5.
Partial-progress target named there: a machine-certified verification for one explicit atom.)

This directory is an honest attempt at that problem. It does **not** solve it. It proves the
N = 2 case of the inequality for every nuclear charge outside a small window, with a
computer-certified proof in exact rational arithmetic, and it explains quantitatively why the same
tools cannot reach the first physically interesting case (N = 3, lithium: ionization potential
>= electron affinity).

## Result

**Theorem.** E(Z,1) + E(Z,3) >= 2 E(Z,2) holds for every real Z in (0, 1] and every real Z >= 1.45.
In particular it holds for every integer nuclear charge Z >= 1.

Equivalently: the first ionization energy of the two-electron ion of charge Z is at least the
ionization energy of the three-electron ion, for all those Z. The N = 1 case of convexity,
E(Z,0) + E(Z,2) >= 2 E(Z,1), is trivial (H(Z,2) >= h_1 + h_2 >= -Z^2), so the theorem settles
convexity at N = 1 and N = 2 for all physical atoms. It also gives the "weakest consequence"
named in the catalogue at N = 3 for the same Z range: if the nucleus cannot bind two electrons
(E(Z,2) = E(Z,1)) then it cannot bind three.

What is proved here is elementary. I could not find this statement in the references listed in the
catalogue entry (Lieb 1983, Nam 2020, Burgess–Linscott–O'Regan 2023), but it may well be folklore;
it should be read as a certified base case, not as progress on the hard part of the conjecture.

The window Z in (1, 1.45) is not covered. Lieb's theorem covers Z <= 1; the lower-bound method used
here fails intrinsically below Z of about 1.35 (see "Limits"), so closing the window needs a
different idea.

## Proof

Write h = -Delta/2 - Z/|x| for the one-electron operator and E(Z,1) = -Z^2/2.

### Step 1 (Z <= 1): Lieb's bound on maximal negative ionization

Lieb (Phys. Rev. A 29, 3018 (1984)) proved that H(Z,N) has no bound state when N >= 2Z + 1.
For Z <= 1 and N = 3 this applies, so by the HVZ theorem E(Z,3) = E(Z,2). The inequality then
reads E(Z,1) >= E(Z,2), which holds because E(Z,N) is nonincreasing in N (put an electron at
infinity).

### Step 2: a spin-resolved Hall–Post bound for three electrons

Write the three-electron Hamiltonian as a sum over pairs,

    H(Z,3) = sum_{i<j} g_ij,       g_ij = (h_i + h_j)/2 + 1/|x_i - x_j| ,

(each one-body term appears in two pairs). The dilation x -> x/2 maps g_ij unitarily to
2 H(Z/2, 2): the pair operator is twice the two-electron atom of charge Z/2.

Let Psi be a normalized antisymmetric three-electron function and let P^s_ij, P^t_ij be the
spin-singlet and spin-triplet projectors of the pair (ij). Since g_ij commutes with spin,
<Psi, g_ij Psi> = <P^s Psi, g_ij P^s Psi> + <P^t Psi, g_ij P^t Psi>. Because Psi is antisymmetric
under the full exchange of i and j, P^s_ij Psi is spatially symmetric in (x_i, x_j) and
P^t_ij Psi is spatially antisymmetric. Hence

    <Psi, H Psi> >= lambda_s sum_{ij} ||P^s_ij Psi||^2 + lambda_t sum_{ij} ||P^t_ij Psi||^2 ,

where lambda_s (lambda_t) is the bottom of the spectrum of g over spatially symmetric
(antisymmetric) pair functions, i.e. lambda_s = 2 E_sym(Z/2, 2), lambda_t = 2 E_anti(Z/2, 2).
The spin identity sum_{i<j} P^t_ij = 9/8 + S^2/2 gives tau := <Psi, sum P^t Psi> in [3/2, 3] and
sum P^s = 3 - tau. The ground state of the two-electron problem over all functions is positive
(Perron–Frobenius), hence symmetric, so lambda_s <= lambda_t; therefore the right-hand side is
nondecreasing in tau and is minimized at tau = 3/2:

    E(Z,3) >= (3/2)(lambda_s + lambda_t) = 3 [ E_sym(Z/2, 2) + E_anti(Z/2, 2) ].        (HP)

(This is the Hall–Post subsystem bound, improved by the spin bookkeeping; without the spin
resolution one only gets 6 E_sym(Z/2, 2), which is too weak here.)

### Step 3: Bazley lower bounds for the two-electron atom

For the two-electron atom H(Z',2) = H_0 + V with H_0 = h_1 + h_2 (spectrum known exactly) and
V = 1/|x_1 - x_2| >= 0, Bazley's method of intermediate problems with the "special choice" gives,
for any finite family of H_0-eigenfunctions u_1..u_K (so that p_i = |x_1 - x_2| u_i satisfies
V p_i = u_i), the operator inequality

    V >= sum_{ij} |u_i> (B^{-1})_ij <u_j| ,       B_ij = <u_i, |x_1 - x_2| u_j>,

which is Bessel's inequality in the Hilbert space with inner product <f, V g>. Hence
H >= H_K := H_0 + sum_ij |u_i>(B^{-1})_ij<u_j|. H_K leaves span{u_i} invariant and equals H_0 on
its orthogonal complement, so, sector by sector,

    E_sector(Z',2) >= min( lowest eigenvalue of diag(eps) + B^{-1} S restricted to span{u},
                           lowest H_0 level in the sector not among the u_i ).

Here eps_i are the H_0 eigenvalues of the u_i and S = diag(<u_i,u_i>) (the u_i are taken
orthogonal but not normalized, which keeps every matrix entry rational). The u_i used are the pair
functions 1s(1) nl(2) +/- nl(1) 1s(2) with n <= K and l <= n-1 (and 1s 1s in the symmetric
sector), for both symmetry sectors; B is block diagonal in l and independent of m. The excluded
H_0 level is then 1s (K+1) l, at -(Z'^2/2)(1 + 1/(K+1)^2), and the continuum starts at -Z'^2/2.
All B and S entries are computed exactly by sympy in `bazley_blocks.py` using the multipole
expansion |x_1 - x_2| = sum_l c_l(r_<, r_>) P_l(cos theta_12) with
c_l = (r_<^l / r_>^{l+1}) ( r_<^2/(2l+3) - r_>^2/(2l-1) ) and hydrogenic radial functions, in the
scaled variable Z' r, so that S ~ Z'^{-6}, B ~ Z'^{-7}, eps ~ Z'^2 and the bound to certify is the
lowest generalized eigenvalue of Z'^2 S diag(e) + Z' S B^{-1} S with metric S.

With K = 1 (u = 1s 1s only; B = 35/(16 Z')) this is the classical bound
E_sym(Z',2) >= min(-Z'^2 + 16 Z'/35, -5 Z'^2/8), and E_anti(Z',2) >= -5 Z'^2/8.

### Step 4: variational upper bound for the two-electron atom

Any normalized trial function gives E(Z,2) <= <psi, H psi>. For the two-parameter Hylleraas
function psi = exp(-zeta s)(1 + c u), s = r_1 + r_2, u = r_12, the Rayleigh quotient is the
rational function (derived and checked in `hylleraas.py`)

    E_hyl(zeta,c;Z) = zeta [ c^2 (35 - 144 Z) + 8 c zeta (8c + 4 - 15 Z) + 16 zeta^3
                              + zeta^2 (50 c + 10 - 32 Z) ] / ( 2 (48 c^2 + 35 c zeta + 8 zeta^2) ).

It is affine and decreasing in Z for fixed (zeta, c). With c = 0 and zeta = Z - 5/16 it reduces to
E(Z,2) <= -(Z - 5/16)^2. (Check: Z = 2 gives -2.8911 at zeta = 1.85, c = 0.365, the textbook value.)

### Step 5: assembling, Z >= 2.054 (analytic)

Using K = 1 in Step 3 and c = 0 in Step 4, E(Z,1) + E(Z,3) - 2E(Z,2) >= min(P1(Z), P2(Z)) with

    P1(Z) = (9/32) Z^2 - (79/140) Z + 25/128,      P2(Z) = (9/16) Z^2 - (5/4) Z + 25/128 .

Both quadratics are positive beyond their larger roots, 1.5617 and 2.0531 respectively; so the
inequality holds for every Z >= 2.054 (`certify.py`, function `analytic_part`, verifies the roots
and signs exactly). For integer Z this already gives margins of 1.03 Ha at Z = 3 growing to 22.7 Ha
at Z = 10.

### Step 6: assembling, 1.45 <= Z <= 2.054 (certified cell by cell)

`certify.py` covers [1.45, 2.054] by 331 cells [Z_a, Z_b] of width 1/2000 to 1/200. On each cell:

* lower bound: for Z' = Z/2 in [a, b], Z'^2 S diag(e) >= b^2 S diag(e) (e < 0) and
  Z' S B^{-1} S >= a S B^{-1} S (B^{-1} > 0), so a single generalized eigenvalue bound mu valid on
  the whole cell is obtained; mu is found in floating point, rounded down to a rational, and the
  matrix inequality A - mu S >= 0 is then certified by exact rational Gaussian elimination. The
  excluded-level term is evaluated at b. Both sectors and all l-blocks for K = 6 are used.
* upper bound: Hylleraas parameters are optimized at Z_a and rationalized; since E_hyl is
  decreasing in Z, E(Z,2) <= E_hyl(zeta,c;Z_a) on the cell, and -E(Z,1) = Z^2/2 <= Z_b^2/2.
* the cell passes if 3[E_sym + E_anti](lower bounds) >= 2 E_hyl + Z_b^2/2, evaluated exactly.

The certificate (`data/certificate_K6.log`) reports the smallest margin, 0.0073 Ha on the first cell
[1.45, 1.4505], growing to 0.37 Ha at Z = 1.8. Together with Step 5 and Step 1 this proves the theorem.

Sanity checks of the lower bounds against known energies (all below the exact values, as required):
Bazley K = 6 gives E_sym(1,2) >= -0.5534 (H^-: -0.5277), E_sym(2,2) >= -3.0608 (He: -2.9037),
E_anti(2,2) >= -2.2151 (He 2^3S: -2.1752); the Hall–Post bound gives E(3,3) >= -8.257 (Li: -7.478).

## Limits of the method, and why N = 3 (lithium) is out of reach this way

The interesting instance is N = 3, i.e. for Li the statement IP(Li) >= EA(Li), with exact numbers
E(2) = -7.2799, E(3) = -7.4781, E(4) = -7.5008 Ha and a true margin of only 0.175 Ha. A
certified proof needs a rigorous lower bound on E(Li^-) within about 0.17 Ha of the truth
(given essentially exact bounds for the other two). I checked what subsystem decompositions can give:

* pair decomposition of H(Z,4) (each pair operator is 3 H(1,2), i.e. H^-): E(3,4) >= 18 E(1,2) = -9.50;
* the spin-resolved version: E(3,4) >= 3 [E_sym + E_anti](1,2) = -9.25;
* triple decomposition (each triple operator is (3/4) H(2,3), i.e. He^-): E(3,4) >= 3 E(2,3) >= -8.71
  even if one takes E(2,3) = E(He) as known.

All are far above the needed -7.676 Ha. At large Z the same happens for every N >= 3: the
fractional distribution of the nuclear attraction over subsystems loses a term of order Z^2
(for N = 3 the bound reads -1.5 Z^2 + O(Z) against a requirement of -1.25 Z^2 + O(Z)). The
Pauli principle for the 1s shell is exactly what these decompositions cannot see. Convexity at N = 2
is accessible only because the required bound, E(Z,3) >= 2E(Z,2) + Z^2/2, is loose by order Z.

For the N = 2 case itself the method also has an intrinsic floor: below Z of about 1.35 even the
exact pair energies make the Hall–Post right-hand side fail (for Z <= 1.82 the pair atom of charge
Z/2 does not bind, so (HP) can at best give -0.75 Z^2, while 2E(Z,2) + Z^2/2 is above that only
if E(Z,2) <= -0.625 Z^2, which fails for Z below roughly 1.35). The finite basis (K = 6) and the
two-parameter upper bound move the floor up to 1.45. Larger K and a longer Hylleraas expansion
would reach about 1.40 but not 1.0.

What a certified lithium check would actually need: a rigorous lower bound for Li^- within
0.17 Ha, which means either a Temple-type bound with an explicitly correlated (Hylleraas or ECG)
four-electron wave function whose energy variance is below about 0.004 Ha^2 (the Temple gap is
only the 0.023 Ha binding energy of Li^-), or a Bazley–Fox intermediate problem for four electrons
with the full interaction, whose matrix elements involve 1/(r_12 r_13)-type three-electron
integrals. Both are substantial projects; neither is close to a session's work, and this is the
genuine obstacle behind the catalogue's remark that no interacting case has ever been certified.

## Files

* `README.md` — this note (statement, proof, limits).
* `hylleraas.py` — derives the closed-form Rayleigh quotient of Step 4 with sympy and checks it.
* `bazley_blocks.py` — exact rational B, S, e blocks of Step 3 for a given K (`python bazley_blocks.py 6`
  writes `data/bazley_K6.json`; about 4 minutes).
* `certify.py` — Steps 5 and 6 in exact arithmetic (`python certify.py 6`; writes
  `data/certificate_K6.log`).
* `data/bazley_K{3,5,6}.json` — the exact blocks; `data/certificate_K6.log` — the certificate.

Dependencies: sympy, numpy, scipy (floating point is used only to guess values that are then
certified exactly).
