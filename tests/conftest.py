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
