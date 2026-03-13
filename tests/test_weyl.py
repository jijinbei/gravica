"""Tests for Weyl tensor."""


def test_minkowski_weyl_zero(mink_weyl):
    """Minkowski is conformally flat: C_{abcd} = 0."""
    W = mink_weyl
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    assert str(W[a, b, c, d]) == "0"


def test_schwarzschild_weyl_nonzero(schw_weyl):
    """Schwarzschild Weyl tensor should be nonzero (not conformally flat)."""
    W = schw_weyl
    # For vacuum, Weyl = Riemann, so at least one component should be nonzero
    has_nonzero = any(
        str(W[a, b, c, d]) != "0"
        for a in range(4)
        for b in range(4)
        for c in range(4)
        for d in range(4)
    )
    assert has_nonzero, "Schwarzschild Weyl tensor should have nonzero components"
