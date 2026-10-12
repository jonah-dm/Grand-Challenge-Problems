"""Small helpers for building molecular geometries (Angstrom)."""
import numpy as np

d2r = np.pi / 180


def linear(*pairs):
    """linear(('H',0), ('C',1.066), ...) -> atoms along z at cumulative positions."""
    z, out = 0.0, []
    for sym, dz in pairs:
        z += dz
        out.append((sym, (0.0, 0.0, z)))
    return out


def diatomic(a, b, r):
    return [(a, (0, 0, 0)), (b, (0, 0, r))]


def bent(c, h, r, ang, h2=None, r2=None):
    """XY2 (or XYZ) bent molecule with central atom c."""
    h2 = h2 or h
    r2 = r2 or r
    a = ang * d2r / 2
    return [(c, (0, 0, 0)), (h, (r * np.sin(a), 0, r * np.cos(a))),
            (h2, (-r2 * np.sin(a), 0, r2 * np.cos(a)))]


def pyramid(c, h, r, ang):
    """XH3 with C3v symmetry and H-X-H angle `ang` (120 -> planar)."""
    ct = np.cos(ang * d2r)
    # cos(theta) = cos^2(b) - 0.5 sin^2(b), b = angle of X-H bond from the C3 axis
    s2 = (1 - ct) / 1.5
    b = np.arcsin(np.sqrt(min(s2, 1.0)))
    out = [(c, (0, 0, 0))]
    for k in range(3):
        p = 2 * np.pi * k / 3
        out.append((h, (r * np.sin(b) * np.cos(p), r * np.sin(b) * np.sin(p), -r * np.cos(b))))
    return out


def tetra(c, h, r):
    s = r / np.sqrt(3)
    return [(c, (0, 0, 0))] + [(h, (s * x, s * y, s * z)) for x, y, z in
                               [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]]


def to_str(atoms):
    return "; ".join(f"{s} {x:.6f} {y:.6f} {z:.6f}" for s, (x, y, z) in atoms)
