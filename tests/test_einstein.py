"""Tests for Einstein tensor."""

from atlas.simplify import simplify


def test_minkowski_einstein_zero(mink_einstein):
    for a in range(4):
        for b in range(4):
            assert str(mink_einstein[a, b]) == "0"


def test_schwarzschild_vacuum(schw_einstein):
    """Schwarzschild: G_{ab} = 0 (vacuum)."""
    for a in range(4):
        for b in range(4):
            val = schw_einstein[a, b]
            assert str(val) == "0", f"G_{{{a}{b}}} = {val}"
