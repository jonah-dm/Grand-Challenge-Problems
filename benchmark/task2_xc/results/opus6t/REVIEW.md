# Organiser review of opus6t (round 6, theory track)

Functional B2: dev WTMAD 1.37 (post-vtau-fix harness), all converged, audit DISTINCT (SCAN0 R²=0.94).
No new construction element (opus5t exchange on sonnet5c base with opus4 correlation; constants
from the H norm). H+H₂ barrier 10.38 (exact 9.6); H₂⁺ at 6 bohr error −29.4 → −5.8 kcal/mol.
Disclosed: phase-1 verification scripts ran before the 60-min mark (no SCF, no functional code).

Theory:
- Occupation-square particle triangle: correctly not claimed (flat-plane / constancy condition).
- C2(ii) hole-triangle bounds on E_c for N ≥ 1 with particle–hole symmetric "tent": plausible;
  agent flags possible overlap with Gould–Dobson 2013 / Gould–Pittalis 2017 ensemble work. UNCHECKED.
- T3: flat-plane axis forces semilocal correlation to a unique s⁻² form on hydrogenic densities,
  wrong by 5× on Gaussians: the correlation analogue of round 5's T2. Only the form is claimed.
- T4 (main): nowhere-positive semilocal correlation + flat-plane axis on hydrogenic ions ⇒
  over-correlation of every two-electron singlet (He ≤ −0.876 Ha vs −0.042 exact). Verified by the
  agent on HF densities and against an exact 1D soft-Coulomb Lieb-functional calculation over the
  whole square (kink at N = 1 that no semilocal form reproduces). Organiser has not re-derived
  it. NEEDS LITERATURE CHECK against Mori-Sánchez–Cohen–Yang 2009 (discontinuous nature of the
  exact functional; semilocal forms cannot satisfy the flat plane) and Cohen–Mori-Sánchez–Yang
  2008 (fractional spins): the qualitative incompatibility is likely known; the quantitative
  bound and the 1D exact-square verification may be new.
- Negative result: candidate C (direct enforcement via whole-density indicator) fails because
  open-shell valence regions mimic fractional-spin one-orbital densities; withdrawn and explained.

Verdict: functional not new; theory T4 to literature check; good scientific conduct.
