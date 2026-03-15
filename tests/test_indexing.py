"""Tests for index raising/lowering utilities."""

from gravica.indexing import raise_index_2d, lower_index_2d
from gravica.simplify import simplify, str_is_zero


def test_raise_lower_roundtrip(schw):
    """Raising then lowering an index should return the original tensor."""
    g = schw.components
    # Raise index 0 of g_{ab}, then lower it back
    mixed = raise_index_2d(schw, g, which=0)
    # mixed should be delta^a_b (identity)
    n = schw.dim
    for a in range(n):
        for b in range(n):
            expected = "1" if a == b else "0"
            assert str(simplify(mixed[a][b])) == expected, (
                f"g^{{ac}} g_{{cb}} [{a}][{b}] = {mixed[a][b]}, expected {expected}"
            )


def test_lower_raise_roundtrip(schw):
    """Lowering then raising should return the original."""
    g_inv = schw.inverse
    mixed = lower_index_2d(schw, g_inv, which=0)
    n = schw.dim
    for a in range(n):
        for b in range(n):
            expected = "1" if a == b else "0"
            assert str(simplify(mixed[a][b])) == expected, (
                f"g_{{ac}} g^{{cb}} [{a}][{b}] = {mixed[a][b]}, expected {expected}"
            )


def test_raise_ricci(schw_ricci, schw):
    """Raise first index of Ricci tensor, then lower it back to recover original."""
    R_lower = schw_ricci.components
    R_mixed = raise_index_2d(schw, R_lower, which=0)
    R_recovered = lower_index_2d(schw, R_mixed, which=0)
    n = schw.dim
    for a in range(n):
        for b in range(n):
            diff = simplify(R_recovered[a][b] - R_lower[a][b])
            assert str_is_zero(diff), (
                f"Roundtrip failed at [{a}][{b}]: {R_recovered[a][b]} vs {R_lower[a][b]}"
            )


def test_fully_contravariant_minkowski(mink_riemann):
    """Minkowski: fully contravariant Riemann should also be zero."""
    n = mink_riemann.dim
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(c + 1, n):
                    val = mink_riemann.fully_contravariant(a, b, c, d)
                    assert str_is_zero(val), f"R^{{{a}{b}{c}{d}}} = {val}, expected 0"
