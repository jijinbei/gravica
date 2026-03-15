"""Tests for the Kretschner scalar."""

from symbolica import Expression, S

from gravica.kretschner import kretschner_scalar
from gravica.simplify import simplify, str_is_zero


def test_kretschner_minkowski(mink_riemann):
    """Minkowski spacetime has zero Kretschner scalar."""
    K = kretschner_scalar(mink_riemann)
    assert str_is_zero(K), f"Expected 0, got {K}"


def test_kretschner_schwarzschild(schw_riemann):
    """Schwarzschild: K = 12 r_s^2 / r^6."""
    K = kretschner_scalar(schw_riemann)
    r = S("r")
    r_s = S("r_s")
    expected = Expression.num(12) * r_s**2 / r**6
    diff = simplify(K - expected)
    assert str_is_zero(diff), f"Expected 12*r_s^2/r^6, got {K}"
