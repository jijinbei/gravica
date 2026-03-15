"""Friedmann equations from the FLRW Einstein tensor.

The Einstein tensor for the FLRW metric encodes the Friedmann equations:

  First Friedmann equation  (G_00):
      3(ȧ² + k) / a² = 8πG ρ

  Acceleration equation  (spatial diagonal):
      (2ä·a + ȧ² + k) / a² = -8πG p

  Ricci scalar:
      R = -6(ä·a + ȧ² + k) / a²

where ȧ = da/dt and ä = d²a/dt².
"""

from symbolica import Expression, S

from gravica.ricci import ricci_scalar
from gravica.simplify import simplify, str_is_zero

ZERO = Expression.num(0)


def test_friedmann_first_equation(flrw_einstein):
    r"""G_{00} = 3(ȧ² + k) / a².

    Multiplying both sides by a² gives:
        G_{00} · a² - 3(ȧ² + k) = 0.
    """
    t = S('t')
    k = S('k')
    a = S('a')(t)
    da = a.derivative(t)       # ȧ = der(1, a(t))

    G_00 = flrw_einstein[0, 0]
    lhs = G_00 * a**2 - Expression.num(3) * (da**2 + k)
    assert str_is_zero(simplify(lhs)), f"First Friedmann eq failed: {simplify(lhs)}"


def test_friedmann_acceleration_equation(flrw_einstein):
    r"""G_{11} · (1 - kr²) = 2ä·a + ȧ² + k.

    This is the spatial component that, combined with G_{00},
    yields the acceleration equation ä/a = -4πG(ρ + 3p)/3.
    """
    t = S('t')
    r = S('r')
    k = S('k')
    a = S('a')(t)
    da = a.derivative(t)       # ȧ
    dda = da.derivative(t)     # ä = der(2, a(t))

    G_11 = flrw_einstein[1, 1]
    # G_11 = (k + 2ä·a + ȧ²) / (-1 + kr²)
    # so G_11 · (-1 + kr²) = k + 2ä·a + ȧ²
    lhs = G_11 * (Expression.num(-1) + k * r**2) - (
        k + Expression.num(2) * dda * a + da**2
    )
    assert str_is_zero(simplify(lhs)), f"Acceleration eq failed: {simplify(lhs)}"


def test_friedmann_ricci_scalar(flrw_ricci):
    r"""R = -6(ä·a + ȧ² + k) / a²."""
    t = S('t')
    k = S('k')
    a = S('a')(t)
    da = a.derivative(t)
    dda = da.derivative(t)

    R = ricci_scalar(flrw_ricci)
    expected = Expression.num(-6) * (dda * a + da**2 + k) / a**2
    diff = simplify(R - expected)
    assert str_is_zero(diff), f"Ricci scalar failed: got {R}, diff = {diff}"


def test_friedmann_off_diagonal_zero(flrw_einstein):
    """All off-diagonal G_ab = 0 (consistent with a perfect fluid)."""
    for a in range(4):
        for b in range(4):
            if a != b:
                val = flrw_einstein[a, b]
                assert str_is_zero(val), f"G[{a},{b}] = {val}, expected 0"


def test_friedmann_spatial_isotropy(flrw_einstein):
    r"""Spatial isotropy: G^i_i are all equal (same pressure in every direction).

    G^1_1 = G^2_2 = G^3_3 = (2ä·a + ȧ² + k) / a².
    We verify G_{ii}/g_{ii} is the same for i = 1, 2, 3.
    """
    t = S('t')
    k = S('k')
    a = S('a')(t)
    da = a.derivative(t)
    dda = da.derivative(t)

    expected_mixed = (
        Expression.num(2) * dda * a + da**2 + k
    ) / a**2

    metric = flrw_einstein.metric
    for i in range(1, 4):
        G_ii = flrw_einstein[i, i]
        g_ii = metric[i, i]
        # G^i_i = G_{ii} g^{ii} = G_{ii} / g_{ii}  (diagonal metric)
        mixed = simplify(G_ii / g_ii)
        diff = simplify(mixed - expected_mixed)
        assert str_is_zero(diff), (
            f"G^{i}_{i} = {mixed}, expected {expected_mixed}"
        )
