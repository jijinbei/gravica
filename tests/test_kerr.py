"""Tests for Kerr metric."""

from gravica.simplify import simplify
from symbolica import Expression, S
from conftest import check_inverse_identity, assert_zero


def test_kerr_metric_inverse(kerr_metric):
    check_inverse_identity(kerr_metric)


def test_kerr_not_diagonal(kerr_metric):
    """Kerr has off-diagonal g_{t,phi} component."""
    assert not kerr_metric.is_diagonal()


def test_kerr_reduces_to_schwarzschild(kerr_metric):
    """Setting a=0 in Kerr should give Schwarzschild metric components."""
    a_sym = S("a")
    zero = Expression.num(0)
    from gravica.metrics.schwarzschild import schwarzschild

    schw = schwarzschild()

    for i in range(4):
        for j in range(4):
            kerr_comp = kerr_metric[i, j].replace(a_sym, zero)
            kerr_comp = simplify(kerr_comp)
            schw_comp = schw[i, j]
            # Schwarzschild uses r_s, Kerr uses 2*M; substitute r_s = 2*M
            r_s = S("r_s")
            M = S("M")
            schw_sub = schw_comp.replace(r_s, Expression.num(2) * M)
            schw_sub = simplify(schw_sub)
            diff = simplify(kerr_comp - schw_sub)
            assert_zero(
                kerr_comp - schw_sub,
                f"g[{i},{j}]: Kerr(a=0)={kerr_comp} != Schw={schw_sub}, diff={diff}",
            )
