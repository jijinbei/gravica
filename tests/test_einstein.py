"""Tests for Einstein tensor."""

from gravica.simplify import str_is_zero


def test_minkowski_einstein_zero(mink_einstein):
    for a in range(4):
        for b in range(4):
            assert str_is_zero(mink_einstein[a, b])


def test_schwarzschild_vacuum(schw_einstein):
    """Schwarzschild: G_{ab} = 0 (vacuum)."""
    for a in range(4):
        for b in range(4):
            val = schw_einstein[a, b]
            assert str_is_zero(val), f"G_{{{a}{b}}} = {val}"
