"""Tests for Ricci tensor and scalar."""

from gravica.ricci import RicciTensor, ricci_scalar
from gravica.simplify import str_is_zero


def test_minkowski_ricci_zero(mink_ricci):
    """All Ricci components vanish for flat spacetime."""
    for a in range(4):
        for b in range(4):
            assert str_is_zero(mink_ricci[a, b])


def test_minkowski_scalar_zero(mink_ricci):
    assert str_is_zero(ricci_scalar(mink_ricci))


def test_schwarzschild_vacuum(schw_ricci):
    """Schwarzschild is a vacuum solution: R_{ab} = 0."""
    for a in range(4):
        for b in range(4):
            val = schw_ricci[a, b]
            assert str_is_zero(val), f"R_{{{a}{b}}} = {val}"


def test_schwarzschild_scalar_zero(schw_ricci):
    """Schwarzschild Ricci scalar = 0."""
    R = ricci_scalar(schw_ricci)
    assert str_is_zero(R)


def test_from_metric(mink):
    """RicciTensor.from_metric should produce the same results."""
    ricci = RicciTensor.from_metric(mink)
    for a in range(4):
        for b in range(4):
            assert str_is_zero(ricci[a, b])
