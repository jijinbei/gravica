"""Tests for FLRW metric."""

from conftest import check_inverse_identity


def test_flrw_metric_inverse(flrw_metric):
    check_inverse_identity(flrw_metric)


def test_flrw_diagonal(flrw_metric):
    assert flrw_metric.is_diagonal()


def test_flrw_christoffel_nonzero(flrw_christoffel):
    """Gamma^r_{tr} = a'(t)/a(t), i.e. der(1,a(t))/a(t) in Symbolica."""
    # Gamma^1_{01} = Gamma^r_{tr}
    gamma_r_tr = flrw_christoffel[1, 0, 1]
    s = str(gamma_r_tr)
    assert s != "0", f"Gamma^r_tr should be nonzero, got {s}"
    # Should contain derivative of a(t)
    assert "der" in s or "a" in s, f"Gamma^r_tr should involve a(t), got {s}"
