"""Tests for metric tensor operations."""

from gravica.metric import MetricTensor, symbolic_det, symbolic_inverse
from gravica.simplify import simplify, is_zero
from conftest import check_inverse_identity


def test_minkowski_inverse(mink):
    check_inverse_identity(mink)


def test_schwarzschild_inverse(schw):
    check_inverse_identity(schw)


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
