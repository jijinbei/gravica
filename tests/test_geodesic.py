"""Tests for geodesic equations."""

from gravica.geodesic import geodesic_equations
from gravica.simplify import str_is_zero


def test_geodesic_minkowski(mink_christoffel):
    """Minkowski: all Christoffel symbols vanish, so equations are just dd<coord>_ = 0."""
    result = geodesic_equations(mink_christoffel)

    # Each equation should be just the acceleration (free particle)
    for i, eq in enumerate(result.equations):
        diff = eq - result.accelerations[i]
        assert str_is_zero(diff), (
            f"Equation {i} should be just acceleration, got {eq}"
        )


def test_geodesic_returns_correct_symbols(mink_christoffel):
    """Check that velocity and acceleration symbols are generated correctly."""
    result = geodesic_equations(mink_christoffel)

    assert len(result.velocities) == 4
    assert len(result.accelerations) == 4
    assert len(result.equations) == 4

    # Check naming convention
    vel_names = [str(v) for v in result.velocities]
    acc_names = [str(a) for a in result.accelerations]
    assert vel_names == ["dt_", "dx_", "dy_", "dz_"]
    assert acc_names == ["ddt_", "ddx_", "ddy_", "ddz_"]


def test_geodesic_schwarzschild_nontrivial(schw_christoffel):
    """Schwarzschild geodesics should have non-trivial Christoffel contributions."""
    result = geodesic_equations(schw_christoffel)

    # At least some equations should have more than just the acceleration term
    has_nontrivial = False
    for i, eq in enumerate(result.equations):
        diff = eq - result.accelerations[i]
        if not str_is_zero(diff):
            has_nontrivial = True
            break
    assert has_nontrivial, "Schwarzschild should have non-trivial geodesic equations"
