"""Out-of-class hidden set: second-row atomization energies, barrier heights,
noncovalent interactions, proton affinities. Geometries approx. experimental (Angstrom)."""
import numpy as np
from geom import linear, diatomic, bent, pyramid, tetra, d2r
from collections import Counter

ATOM_SPIN = {"H": 1, "C": 2, "N": 3, "O": 2, "Si": 2, "P": 3, "S": 2, "Cl": 1}
SPECIES = {f"{a}_atom": ([(a, (0, 0, 0))], 0, s) for a, s in ATOM_SPIN.items()}

# water dimer (Cs, approx. CCSD(T) geometry, donor O at origin)
H2O_dimer = [("O", (-1.551007, -0.114520, 0.0)), ("H", (-1.934259, 0.762503, 0.0)),
             ("H", (-0.599677, 0.040712, 0.0)),
             ("O", (1.350625, 0.111469, 0.0)), ("H", (1.680398, -0.373741, -0.758561)),
             ("H", (1.680398, -0.373741, 0.758561))]
H2O = [("O", (0, 0, 0)), ("H", (0.957, 0, 0)), ("H", (-0.24, 0.927, 0))]  # not used; see dev
NH3 = pyramid("N", "H", 1.012, 106.7)
NH4 = tetra("N", "H", 1.024)
H3O = pyramid("O", "H", 0.976, 111.3)
H3_TS = linear(("H", 0), ("H", 0.930), ("H", 0.930))
SO2 = bent("S", "O", 1.431, 119.3)

SPECIES.update({
    "SiH4": (tetra("Si", "H", 1.480), 0, 0),
    "PH3": (pyramid("P", "H", 1.420, 93.5), 0, 0),
    "H2S": (bent("S", "H", 1.336, 92.1), 0, 0),
    "HCl": (diatomic("Cl", "H", 1.275), 0, 0),
    "Cl2": (diatomic("Cl", "Cl", 1.988), 0, 0),
    "SO2": (SO2, 0, 0),
    "CS": (diatomic("C", "S", 1.535), 0, 0),
    "H2": (diatomic("H", "H", 0.741), 0, 0),
    "H3_TS": (H3_TS, 0, 1),
    "H2O_dimer": (H2O_dimer, 0, 0),
    "H2O_a": (H2O_dimer[:3], 0, 0),
    "H2O_b": (H2O_dimer[3:], 0, 0),
    "NH3": (NH3, 0, 0),
    "NH4_cation": (NH4, 1, 0),
    "H2O": (bent("O", "H", 0.958, 104.5), 0, 0),
    "H3O_cation": (H3O, 1, 0),
})

def atomization(name):
    expr = {f"{s}_atom": c for s, c in Counter(a for a, _ in SPECIES[name][0]).items()}
    expr[name] = -1
    return expr

ITEMS = {f"AE2:{m}": ("AE2", atomization(m)) for m in ["SiH4", "PH3", "H2S", "HCl", "Cl2", "SO2", "CS"]}
ITEMS.update({
    "BH:H+H2": ("BH", {"H3_TS": 1, "H_atom": -1, "H2": -1}),
    "NC:H2O_dimer": ("NC", {"H2O_a": 1, "H2O_b": 1, "H2O_dimer": -1}),
    "PA:NH3": ("PA", {"NH3": 1, "NH4_cation": -1}),
    "PA:H2O": ("PA", {"H2O": 1, "H3O_cation": -1}),
})
