"""Tests for Bianchi identities and Riemann tensor symmetries.

Verified on Schwarzschild spacetime, which has nontrivial curvature.
"""

from atlas.simplify import simplify


def test_first_bianchi_identity(schw_riemann):
    """First Bianchi identity: R_{a[bcd]} = 0.

    R_{abcd} + R_{acdb} + R_{adbc} = 0 for all index combinations.
    """
    R = schw_riemann
    n = R.dim
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    r1 = R.fully_covariant(a, b, c, d)
                    r2 = R.fully_covariant(a, c, d, b)
                    r3 = R.fully_covariant(a, d, b, c)
                    total = simplify(r1 + r2 + r3)
                    assert str(total) == "0", (
                        f"Bianchi failed: R_{{{a}{b}{c}{d}}} + R_{{{a}{c}{d}{b}}} "
                        f"+ R_{{{a}{d}{b}{c}}} = {total}"
                    )


def test_riemann_pair_symmetry(schw_riemann):
    """Pair symmetry: R_{abcd} = R_{cdab}."""
    R = schw_riemann
    n = R.dim
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    r_abcd = R.fully_covariant(a, b, c, d)
                    r_cdab = R.fully_covariant(c, d, a, b)
                    diff = simplify(r_abcd - r_cdab)
                    assert str(diff) == "0", (
                        f"R_{{{a}{b}{c}{d}}} != R_{{{c}{d}{a}{b}}}: diff={diff}"
                    )


def test_riemann_first_pair_antisymmetry(schw_riemann):
    """First pair antisymmetry: R_{abcd} = -R_{bacd}."""
    R = schw_riemann
    n = R.dim
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    r_abcd = R.fully_covariant(a, b, c, d)
                    r_bacd = R.fully_covariant(b, a, c, d)
                    total = simplify(r_abcd + r_bacd)
                    assert str(total) == "0", (
                        f"R_{{{a}{b}{c}{d}}} + R_{{{b}{a}{c}{d}}} = {total} != 0"
                    )
