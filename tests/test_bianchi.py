"""Tests for Bianchi identities and Riemann tensor symmetries.

Verified on Schwarzschild spacetime, which has nontrivial curvature.
"""

from conftest import assert_zero


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
                    assert_zero(
                        r1 + r2 + r3,
                        f"Bianchi failed: R_{{{a}{b}{c}{d}}} + R_{{{a}{c}{d}{b}}} "
                        f"+ R_{{{a}{d}{b}{c}}} = {r1 + r2 + r3}",
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
                    assert_zero(
                        r_abcd - r_cdab,
                        f"R_{{{a}{b}{c}{d}}} != R_{{{c}{d}{a}{b}}}: diff={r_abcd - r_cdab}",
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
                    assert_zero(
                        r_abcd + r_bacd,
                        f"R_{{{a}{b}{c}{d}}} + R_{{{b}{a}{c}{d}}} = {r_abcd + r_bacd} != 0",
                    )
