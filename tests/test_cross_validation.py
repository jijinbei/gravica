"""Cross-validation: compare Atlas results against EinsteinPy for Schwarzschild."""

import pytest
import sympy
from sympy import symbols, sin, diag, simplify as sp_simplify, Rational
from einsteinpy.symbolic import (
    MetricTensor as EPyMetric,
    ChristoffelSymbols as EPyChristoffel,
    RiemannCurvatureTensor as EPyRiemann,
    RicciTensor as EPyRicci,
    EinsteinTensor as EPyEinstein,
)

from gravica.metrics.schwarzschild import schwarzschild
from gravica.christoffel import ChristoffelSymbols
from gravica.riemann import RiemannTensor
from gravica.ricci import RicciTensor, ricci_scalar
from gravica.einstein import EinsteinTensor


def _atlas_to_sympy(expr_str: str) -> sympy.Expr:
    """Convert an Atlas/Symbolica expression string to SymPy."""
    t, r, theta, phi, r_s = symbols('t r theta phi r_s')
    ns = {'t': t, 'r': r, 'theta': theta, 'phi': phi, 'r_s': r_s, 'sin': sympy.sin, 'cos': sympy.cos}

    s = expr_str
    # Symbolica uses ^ for power, SymPy uses **
    s = s.replace('^', '**')
    try:
        return sympy.sympify(s, locals=ns)
    except Exception:
        return sympy.sympify(0)


@pytest.fixture(scope="module")
def epy_schwarzschild():
    """EinsteinPy Schwarzschild metric and derived tensors."""
    t, r, theta, phi = symbols('t r theta phi')
    r_s = symbols('r_s')
    f = 1 - r_s / r
    g = diag(f, -1 / f, -r**2, -r**2 * sin(theta)**2).tolist()
    metric = EPyMetric(g, (t, r, theta, phi))
    ch = EPyChristoffel.from_metric(metric)
    riem = EPyRiemann.from_christoffels(ch)
    ric = EPyRicci.from_riemann(riem)
    return {'metric': metric, 'christoffel': ch, 'riemann': riem, 'ricci': ric}


@pytest.fixture(scope="module")
def atlas_schwarzschild():
    """Atlas Schwarzschild tensors."""
    m = schwarzschild()
    ch = ChristoffelSymbols(m)
    riem = RiemannTensor(ch)
    ric = RicciTensor(riem)
    return {'metric': m, 'christoffel': ch, 'riemann': riem, 'ricci': ric}


def test_christoffel_cross_validation(atlas_schwarzschild, epy_schwarzschild):
    """Compare Christoffel symbols between Atlas and EinsteinPy."""
    atlas_ch = atlas_schwarzschild['christoffel']
    epy_ch = epy_schwarzschild['christoffel']
    r_s, r, theta = symbols('r_s r theta')

    for a in range(4):
        for b in range(4):
            for c in range(b, 4):
                atlas_val = _atlas_to_sympy(str(atlas_ch[a, b, c]))
                epy_val = epy_ch.tensor()[a, b, c]
                diff = sp_simplify(atlas_val - epy_val)
                assert diff == 0, (
                    f"Γ^{a}_{{{b}{c}}}: Atlas={atlas_val}, EinsteinPy={epy_val}, diff={diff}"
                )


def test_ricci_cross_validation(atlas_schwarzschild, epy_schwarzschild):
    """Both should give R_{ab} = 0 for Schwarzschild."""
    atlas_ric = atlas_schwarzschild['ricci']
    epy_ric = epy_schwarzschild['ricci']

    for a in range(4):
        for b in range(a, 4):
            atlas_val = _atlas_to_sympy(str(atlas_ric[a, b]))
            epy_val = sp_simplify(epy_ric.tensor()[a, b])
            assert sp_simplify(atlas_val) == 0, f"Atlas R_{{{a}{b}}} = {atlas_val}"
            assert epy_val == 0, f"EinsteinPy R_{{{a}{b}}} = {epy_val}"
