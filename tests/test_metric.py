"""Tests for metric tensor operations."""

from atlas.metric import MetricTensor, symbolic_det, symbolic_inverse
from atlas.simplify import simplify, is_zero
from symbolica import Expression


def _check_inverse_identity(metric: MetricTensor):
    """Verify g_{ac} g^{cb} = δ^b_a."""
    n = metric.dim
    g = metric.components
    g_inv = metric.inverse
    for a in range(n):
        for b in range(n):
            val = Expression.num(0)
            for c in range(n):
                val = val + g[a][c] * g_inv[c][b]
            val = simplify(val)
            expected = "1" if a == b else "0"
            assert str(val) == expected, f"g_{{a{a}c}} g^{{c{b}}} = {val}, expected {expected}"


def test_minkowski_inverse(mink):
    _check_inverse_identity(mink)


def test_schwarzschild_inverse(schw):
    _check_inverse_identity(schw)


def test_minkowski_det(mink):
    assert str(mink.det) == "-1"


def test_schwarzschild_det(schw):
    # det(g) = -r^4 sin²θ for Schwarzschild
    det = schw.det
    det_str = str(simplify(det))
    assert "sin" in det_str and "r" in det_str


def test_minkowski_diagonal(mink):
    assert mink.is_diagonal()


def test_schwarzschild_diagonal(schw):
    assert schw.is_diagonal()
