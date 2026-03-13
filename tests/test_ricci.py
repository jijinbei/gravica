"""Tests for Ricci tensor and scalar."""

from atlas.ricci import ricci_scalar
from atlas.simplify import simplify


def test_minkowski_ricci_zero(mink_ricci):
    """All Ricci components vanish for flat spacetime."""
    for a in range(4):
        for b in range(4):
            assert str(mink_ricci[a, b]) == "0"


def test_minkowski_scalar_zero(mink_ricci):
    assert str(ricci_scalar(mink_ricci)) == "0"


def test_schwarzschild_vacuum(schw_ricci):
    """Schwarzschild is a vacuum solution: R_{ab} = 0."""
    for a in range(4):
        for b in range(4):
            val = schw_ricci[a, b]
            assert str(val) == "0", f"R_{{{a}{b}}} = {val}"


def test_schwarzschild_scalar_zero(schw_ricci):
    """Schwarzschild Ricci scalar = 0."""
    R = ricci_scalar(schw_ricci)
    assert str(R) == "0"
