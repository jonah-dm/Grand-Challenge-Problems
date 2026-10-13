"""Hidden test set. Same classes as the dev set, different molecules. Geometries approx. experimental (Angstrom)."""
import numpy as np
from geom import linear, diatomic, bent, pyramid, tetra, d2r
from collections import Counter

ATOM_SPIN = {"H": 1, "Li": 1, "B": 1, "C": 2, "N": 3, "O": 2, "F": 1}
SPECIES = {f"{a}_atom": ([(a, (0, 0, 0))], 0, s) for a, s in ATOM_SPIN.items()}

def planar3(c, h, r):
    return [(c, (0, 0, 0))] + [(h, (r * np.cos(2 * np.pi * k / 3), r * np.sin(2 * np.pi * k / 3), 0)) for k in range(3)]

s, c = np.sin(121.3 * d2r / 2), np.cos(121.3 * d2r / 2)
C2H4 = [("C", (0, 0, 0.6695)), ("C", (0, 0, -0.6695))] + \
       [("H", (1.085 * np.sin(121.3 * d2r), 0, 0.6695 - 1.085 * np.cos(121.3 * d2r)))] + \
       [("H", (-1.085 * np.sin(121.3 * d2r), 0, 0.6695 - 1.085 * np.cos(121.3 * d2r)))] + \
       [("H", (1.085 * np.sin(121.3 * d2r), 0, -0.6695 + 1.085 * np.cos(121.3 * d2r)))] + \
       [("H", (-1.085 * np.sin(121.3 * d2r), 0, -0.6695 + 1.085 * np.cos(121.3 * d2r)))]
# H2O2: O-O 1.452, O-H 0.965, OOH 100.0, dihedral 118.5
oo, oh, a, dih = 1.452, 0.965, 100.0 * d2r, 118.5 * d2r
H2O2 = [("O", (0, 0, oo / 2)), ("O", (0, 0, -oo / 2)),
        ("H", (oh * np.sin(a), 0, oo / 2 - oh * np.cos(a))),
        ("H", (oh * np.sin(a) * np.cos(dih), oh * np.sin(a) * np.sin(dih), -oo / 2 + oh * np.cos(a)))]
# CH3OH: staggered, approximate
CH3OH = [("C", (0.0470, 0.6617, 0)), ("O", (0.0470, -0.7570, 0)), ("H", (-0.8950, -1.0620, 0)),
         ("H", (-0.4530, 1.0790, 0.8880)), ("H", (-0.4530, 1.0790, -0.8880)), ("H", (1.0860, 0.9820, 0))]
HNO = [("N", (0, 0, 0)), ("O", (0, 0, 1.212)), ("H", (1.063 * np.sin(108.6 * d2r), 0, -1.063 * np.cos(108.6 * d2r)))]

SPECIES.update({
    "C2H4": (C2H4, 0, 0),
    "CO2": (linear(("O", 0), ("C", 1.162), ("O", 1.162)), 0, 0),
    "H2O2": (H2O2, 0, 0),
    "O2_triplet": (diatomic("O", "O", 1.208), 0, 2),
    "NO": (diatomic("N", "O", 1.151), 0, 1),
    "CH3": (planar3("C", "H", 1.079), 0, 1),
    "NH2": (bent("N", "H", 1.024, 103.4), 0, 1),
    "CH3OH": (CH3OH, 0, 0),
    "HNO": (HNO, 0, 0),
    "LiF": (diatomic("Li", "F", 1.564), 0, 0),
    "BH3": (planar3("B", "H", 1.190), 0, 0),
    "BF": (diatomic("B", "F", 1.263), 0, 0),
    # vertical IPs at the neutral geometry
    "C_cation": ([("C", (0, 0, 0))], 1, 1),
    "N_cation": ([("N", (0, 0, 0))], 1, 2),
    "O_cation": ([("O", (0, 0, 0))], 1, 3),
    "H2O_cation": (bent("O", "H", 0.958, 104.5), 1, 1),
    "NH3_cation": (pyramid("N", "H", 1.012, 106.7), 1, 1),
    "H2O": (bent("O", "H", 0.958, 104.5), 0, 0),
    "NH3": (pyramid("N", "H", 1.012, 106.7), 0, 0),
})

def atomization(name):
    expr = {f"{s}_atom": c for s, c in Counter(a for a, _ in SPECIES[name][0]).items()}
    expr[name] = -1
    return expr

AE = ["C2H4", "CO2", "H2O2", "O2_triplet", "NO", "CH3", "NH2", "CH3OH", "HNO", "LiF", "BH3", "BF"]
ITEMS = {f"AE:{m}": ("AE", atomization(m)) for m in AE}
ITEMS.update({
    "IP:C": ("IP", {"C_cation": 1, "C_atom": -1}),
    "IP:N": ("IP", {"N_cation": 1, "N_atom": -1}),
    "IP:O": ("IP", {"O_cation": 1, "O_atom": -1}),
    "IP:H2O": ("IP", {"H2O_cation": 1, "H2O": -1}),
    "IP:NH3": ("IP", {"NH3_cation": 1, "NH3": -1}),
})
