"""Tests for the Schouten tensor."""

from gravica.schouten import SchoutenTensor
from gravica.simplify import str_is_zero


def test_schouten_minkowski(mink_ricci):
    """Minkowski: all Schouten components vanish."""
    S = SchoutenTensor(mink_ricci)
    n = S.dim
    for a in range(n):
        for b in range(n):
            assert str_is_zero(S[a, b]), f"S[{a},{b}] = {S[a, b]}, expected 0"


def test_schouten_schwarzschild(schw_ricci):
    """Schwarzschild is vacuum (R_ab = 0, R = 0), so Schouten vanishes."""
    S = SchoutenTensor(schw_ricci)
    n = S.dim
    for a in range(n):
        for b in range(n):
            assert str_is_zero(S[a, b]), f"S[{a},{b}] = {S[a, b]}, expected 0"


def test_schouten_symmetric(ds_ricci):
    """Schouten tensor should be symmetric: S_ab = S_ba."""
    from gravica.simplify import simplify
    S = SchoutenTensor(ds_ricci)
    n = S.dim
    for a in range(n):
        for b in range(a + 1, n):
            diff = simplify(S[a, b] - S[b, a])
            assert str_is_zero(diff), f"S[{a},{b}] != S[{b},{a}]"
