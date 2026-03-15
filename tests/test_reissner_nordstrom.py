"""Tests for the Reissner-Nordström metric."""

from symbolica import Expression, S

from gravica.ricci import ricci_scalar
from gravica.simplify import simplify, str_is_zero
from conftest import assert_zero, check_inverse_identity


def test_rn_dimensions(rn_metric):
    assert rn_metric.dim == 4


def test_rn_diagonal(rn_metric):
    assert rn_metric.is_diagonal()


def test_rn_inverse(rn_metric):
    check_inverse_identity(rn_metric)


def test_rn_reduces_to_schwarzschild(rn_metric):
    """When r_Q → 0, should reduce to Schwarzschild."""
    r = S('r')
    r_s = S('r_s')
    r_Q = S('r_Q')

    g_tt = rn_metric[0, 0]
    # Replace r_Q with 0
    g_tt_schw = g_tt.replace(r_Q, Expression.num(0))
    expected = (r - r_s) / r
    diff = simplify(g_tt_schw - expected)
    assert str_is_zero(diff), f"g_tt with r_Q=0 should be (r-r_s)/r, got {g_tt_schw}"


def test_rn_einstein_nonzero(rn_einstein):
    """Reissner-Nordström is NOT vacuum — Einstein tensor is nonzero (electromagnetic stress-energy)."""
    n = rn_einstein.dim
    has_nonzero = False
    for a in range(n):
        for b in range(n):
            if not str_is_zero(rn_einstein[a, b]):
                has_nonzero = True
                break
    assert has_nonzero, "RN Einstein tensor should be nonzero"
