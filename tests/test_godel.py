"""Tests for Gödel metric."""

from atlas.simplify import simplify
from symbolica import Expression, S


def _check_inverse_identity(metric):
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


def test_godel_metric_inverse(godel_metric):
    _check_inverse_identity(godel_metric)


def test_godel_not_diagonal(godel_metric):
    """Gödel has off-diagonal g_{ty} component."""
    assert not godel_metric.is_diagonal()


def test_godel_christoffel_nonzero(godel_christoffel):
    """Gödel should have nontrivial Christoffel symbols."""
    has_nonzero = False
    for a in range(4):
        for b in range(4):
            for c in range(b, 4):
                if str(godel_christoffel[a, b, c]) != "0":
                    has_nonzero = True
                    break
            if has_nonzero:
                break
        if has_nonzero:
            break
    assert has_nonzero, "Gödel should have at least one nonzero Christoffel symbol"


def test_godel_ricci_scalar_constant(godel_ricci, godel_metric):
    """Ricci scalar for Gödel should be R = 1/a².

    Symbolica treats exp(2x) and exp(x)² as independent, so we substitute
    exp(2x) -> exp(x)² to allow cancellation.
    """
    from atlas.ricci import ricci_scalar

    R = ricci_scalar(godel_ricci)
    R_str = str(R)
    assert R_str != "0", f"Ricci scalar should be nonzero for Gödel, got {R_str}"

    x_sym = S('x')
    a_sym = S('a')
    exp_sym = S('exp')

    # Replace exp(2*x) -> exp(x)^2 (mathematically identical)
    e2x = exp_sym(Expression.num(2) * x_sym)
    ex2 = exp_sym(x_sym) ** 2
    R_simplified = simplify(R.replace(e2x, ex2))

    expected = Expression.num(1) / a_sym ** 2
    diff = simplify(R_simplified - expected)
    assert str(diff) == "0", (
        f"Ricci scalar: got {R_simplified}, expected 1/a², diff={diff}"
    )
