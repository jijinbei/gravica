"""Tests for Riemann curvature tensor."""

from atlas.simplify import simplify, str_is_zero
from conftest import assert_zero


def test_minkowski_all_zero(mink_riemann):
    """All Riemann components vanish for flat spacetime."""
    R = mink_riemann
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    assert str_is_zero(R[a, b, c, d])


def test_schwarzschild_antisymmetry_cd(schw_riemann):
    """R^a_{bcd} = -R^a_{bdc} (antisymmetric in last two indices)."""
    R = schw_riemann
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(c + 1, 4):
                    assert_zero(
                        R[a, b, c, d] + R[a, b, d, c],
                        f"R^{a}_{{{b}{c}{d}}} + R^{a}_{{{b}{d}{c}}} != 0",
                    )


def test_schwarzschild_nonzero_component(schw_riemann):
    """R^t_{rtr} should be nonzero for Schwarzschild."""
    R = schw_riemann
    assert not str_is_zero(R[0, 1, 0, 1])
