"""Shared fixtures for GR tests."""

import pytest
from symbolica import Expression
from gravica.metric import MetricTensor
from gravica.christoffel import ChristoffelSymbols
from gravica.riemann import RiemannTensor
from gravica.ricci import RicciTensor
from gravica.einstein import EinsteinTensor
from gravica.weyl import WeylTensor
from gravica.simplify import simplify, str_is_zero
from gravica.metrics.minkowski import minkowski
from gravica.metrics.schwarzschild import schwarzschild
from gravica.metrics.flrw import flrw
from gravica.metrics.kerr import kerr
from gravica.metrics.godel import godel
from gravica.metrics.reissner_nordstrom import reissner_nordstrom
from gravica.metrics.de_sitter import de_sitter
from gravica.metrics.anti_de_sitter import anti_de_sitter


def assert_zero(expr, msg=""):
    """Assert an expression is symbolically zero after simplification."""
    val = simplify(expr)
    assert str_is_zero(val), msg or f"Expected 0, got {val}"


def check_inverse_identity(metric: MetricTensor):
    """Verify g_{ac} g^{cb} = delta^b_a."""
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
            assert str(val) == expected, f"g g^-1 [{a}][{b}] = {val}, expected {expected}"


@pytest.fixture(scope="session")
def mink():
    return minkowski()


@pytest.fixture(scope="session")
def schw():
    return schwarzschild()


@pytest.fixture(scope="session")
def mink_christoffel(mink):
    return ChristoffelSymbols(mink)


@pytest.fixture(scope="session")
def schw_christoffel(schw):
    return ChristoffelSymbols(schw)


@pytest.fixture(scope="session")
def mink_riemann(mink_christoffel):
    return RiemannTensor(mink_christoffel)


@pytest.fixture(scope="session")
def schw_riemann(schw_christoffel):
    return RiemannTensor(schw_christoffel)


@pytest.fixture(scope="session")
def mink_ricci(mink_riemann):
    return RicciTensor(mink_riemann)


@pytest.fixture(scope="session")
def schw_ricci(schw_riemann):
    return RicciTensor(schw_riemann)


@pytest.fixture(scope="session")
def mink_einstein(mink_ricci):
    return EinsteinTensor(mink_ricci)


@pytest.fixture(scope="session")
def schw_einstein(schw_ricci):
    return EinsteinTensor(schw_ricci)


@pytest.fixture(scope="session")
def mink_weyl(mink_riemann, mink_ricci):
    return WeylTensor(mink_riemann, mink_ricci)


@pytest.fixture(scope="session")
def schw_weyl(schw_riemann, schw_ricci):
    return WeylTensor(schw_riemann, schw_ricci)


# --- FLRW ---

@pytest.fixture(scope="session")
def flrw_metric():
    return flrw()


@pytest.fixture(scope="session")
def flrw_christoffel(flrw_metric):
    return ChristoffelSymbols(flrw_metric)


@pytest.fixture(scope="session")
def flrw_riemann(flrw_christoffel):
    return RiemannTensor(flrw_christoffel)


@pytest.fixture(scope="session")
def flrw_ricci(flrw_riemann):
    return RicciTensor(flrw_riemann)


@pytest.fixture(scope="session")
def flrw_einstein(flrw_ricci):
    return EinsteinTensor(flrw_ricci)


# --- Kerr ---

@pytest.fixture(scope="session")
def kerr_metric():
    return kerr()


# --- Gödel ---

@pytest.fixture(scope="session")
def godel_metric():
    return godel()


@pytest.fixture(scope="session")
def godel_christoffel(godel_metric):
    return ChristoffelSymbols(godel_metric)


@pytest.fixture(scope="session")
def godel_riemann(godel_christoffel):
    return RiemannTensor(godel_christoffel)


@pytest.fixture(scope="session")
def godel_ricci(godel_riemann):
    return RicciTensor(godel_riemann)


# --- Reissner-Nordström ---

@pytest.fixture(scope="session")
def rn_metric():
    return reissner_nordstrom()


@pytest.fixture(scope="session")
def rn_christoffel(rn_metric):
    return ChristoffelSymbols(rn_metric)


@pytest.fixture(scope="session")
def rn_riemann(rn_christoffel):
    return RiemannTensor(rn_christoffel)


@pytest.fixture(scope="session")
def rn_ricci(rn_riemann):
    return RicciTensor(rn_riemann)


@pytest.fixture(scope="session")
def rn_einstein(rn_ricci):
    return EinsteinTensor(rn_ricci)


# --- de Sitter ---

@pytest.fixture(scope="session")
def ds_metric():
    return de_sitter()


@pytest.fixture(scope="session")
def ds_christoffel(ds_metric):
    return ChristoffelSymbols(ds_metric)


@pytest.fixture(scope="session")
def ds_riemann(ds_christoffel):
    return RiemannTensor(ds_christoffel)


@pytest.fixture(scope="session")
def ds_ricci(ds_riemann):
    return RicciTensor(ds_riemann)


# --- anti-de Sitter ---

@pytest.fixture(scope="session")
def ads_metric():
    return anti_de_sitter()


@pytest.fixture(scope="session")
def ads_christoffel(ads_metric):
    return ChristoffelSymbols(ads_metric)


@pytest.fixture(scope="session")
def ads_riemann(ads_christoffel):
    return RiemannTensor(ads_christoffel)


@pytest.fixture(scope="session")
def ads_ricci(ads_riemann):
    return RicciTensor(ads_riemann)
