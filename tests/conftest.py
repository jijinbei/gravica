"""Shared fixtures for GR tests."""

import pytest
from atlas.metric import MetricTensor
from atlas.christoffel import ChristoffelSymbols
from atlas.riemann import RiemannTensor
from atlas.ricci import RicciTensor
from atlas.einstein import EinsteinTensor
from atlas.weyl import WeylTensor
from atlas.metrics.minkowski import minkowski
from atlas.metrics.schwarzschild import schwarzschild
from atlas.metrics.flrw import flrw
from atlas.metrics.kerr import kerr
from atlas.metrics.godel import godel


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
