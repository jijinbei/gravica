"""Tests for the Stress-Energy-Momentum tensor."""

from symbolica import S

from gravica.stress_energy import StressEnergyTensor
from gravica.simplify import simplify, str_is_zero


def test_stress_energy_minkowski(mink_einstein):
    """Minkowski: T_ab = 0 (vacuum, no cosmological constant)."""
    T = StressEnergyTensor(mink_einstein)
    n = T.dim
    for a in range(n):
        for b in range(n):
            assert str_is_zero(T[a, b]), f"T[{a},{b}] = {T[a, b]}, expected 0"


def test_stress_energy_schwarzschild(schw_einstein):
    """Schwarzschild: T_ab = 0 (vacuum solution)."""
    T = StressEnergyTensor(schw_einstein)
    n = T.dim
    for a in range(n):
        for b in range(n):
            assert str_is_zero(T[a, b]), f"T[{a},{b}] = {T[a, b]}, expected 0"


def test_stress_energy_with_cosmological_constant(mink_einstein, mink):
    """With Λ, T_ab = Λ g_ab for flat spacetime (G_ab = 0)."""
    Lambda = S("Lambda")
    T = StressEnergyTensor(mink_einstein, cosmological_constant=Lambda)
    g = mink.components
    n = T.dim
    for a in range(n):
        for b in range(n):
            expected = Lambda * g[a][b]
            diff = simplify(T[a, b] - expected)
            assert str_is_zero(diff), (
                f"T[{a},{b}] = {T[a, b]}, expected Λ·g[{a},{b}] = {expected}"
            )


def test_stress_energy_symmetric(rn_einstein):
    """Stress-energy tensor should be symmetric."""
    T = StressEnergyTensor(rn_einstein)
    n = T.dim
    for a in range(n):
        for b in range(a + 1, n):
            diff = simplify(T[a, b] - T[b, a])
            assert str_is_zero(diff), f"T[{a},{b}] != T[{b},{a}]"
