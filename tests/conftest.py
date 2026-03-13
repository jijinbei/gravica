"""Shared fixtures for GR tests."""

import pytest
from symbolica import Expression
from atlas.metric import MetricTensor
from atlas.christoffel import ChristoffelSymbols
from atlas.riemann import RiemannTensor
from atlas.ricci import RicciTensor
from atlas.einstein import EinsteinTensor
from atlas.weyl import WeylTensor
from atlas.simplify import simplify, str_is_zero
from atlas.metrics.minkowski import minkowski
from atlas.metrics.schwarzschild import schwarzschild
from atlas.metrics.flrw import flrw
from atlas.metrics.kerr import kerr
from atlas.metrics.godel import godel


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
