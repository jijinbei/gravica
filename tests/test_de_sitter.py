"""Tests for the de Sitter metric."""

from symbolica import Expression, S

from gravica.ricci import ricci_scalar
from gravica.simplify import simplify, str_is_zero
from tests.conftest import assert_zero, check_inverse_identity


def test_ds_dimensions(ds_metric):
    assert ds_metric.dim == 4


def test_ds_diagonal(ds_metric):
    assert ds_metric.is_diagonal()


def test_ds_inverse(ds_metric):
    check_inverse_identity(ds_metric)


def test_ds_ricci_scalar(ds_ricci):
    """de Sitter with (+,-,-,-) signature: R = -4Λ."""
    R = ricci_scalar(ds_ricci)
    Lambda = S('Lambda')
    expected = Expression.num(-4) * Lambda
    diff = simplify(R - expected)
    assert str_is_zero(diff), f"Expected R = 4Λ, got {R}"


def test_ds_reduces_to_minkowski(ds_metric):
    """When Λ → 0, should reduce to Minkowski."""
    Lambda = S('Lambda')
    g_tt = ds_metric[0, 0]
    g_tt_flat = simplify(g_tt.replace(Lambda, Expression.num(0)))
    assert str(g_tt_flat) == "1", f"g_tt with Λ=0 should be 1, got {g_tt_flat}"
