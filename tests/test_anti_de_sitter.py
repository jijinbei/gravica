"""Tests for the anti-de Sitter metric."""

from symbolica import Expression, S

from gravica.ricci import ricci_scalar
from gravica.simplify import simplify, str_is_zero
from conftest import assert_zero, check_inverse_identity


def test_ads_dimensions(ads_metric):
    assert ads_metric.dim == 4


def test_ads_diagonal(ads_metric):
    assert ads_metric.is_diagonal()


def test_ads_inverse(ads_metric):
    check_inverse_identity(ads_metric)


def test_ads_ricci_scalar(ads_ricci):
    """Anti-de Sitter with (+,-,-,-) signature: R = 12/l²."""
    R = ricci_scalar(ads_ricci)
    l = S('l')
    expected = Expression.num(12) / l**2
    diff = simplify(R - expected)
    assert str_is_zero(diff), f"Expected R = -12/l², got {R}"


def test_ads_reduces_to_minkowski(ads_metric):
    """When l → ∞ (r²/l² → 0), the metric function f → 1."""
    l = S('l')
    g_tt = ads_metric[0, 0]
    # f = 1 + r²/l²; as l → ∞, f → 1 (Minkowski)
    # Check structure: f = (l² + r²) / l²
    r = S('r')
    expected_f = (l**2 + r**2) / l**2
    diff = simplify(g_tt - expected_f)
    assert str_is_zero(diff), f"g_tt should be (l²+r²)/l², got {g_tt}"
