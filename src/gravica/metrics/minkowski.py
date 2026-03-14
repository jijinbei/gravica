"""Minkowski (flat spacetime) metric."""

from symbolica import Expression, S

from gravica.metric import MetricTensor, ZERO

ONE = Expression.num(1)
NEG = Expression.num(-1)


def minkowski() -> MetricTensor:
    """η_{ab} = diag(1, -1, -1, -1) in Cartesian coordinates (t, x, y, z)."""
    t, x, y, z = S('t'), S('x'), S('y'), S('z')
    g = [
        [ONE, ZERO, ZERO, ZERO],
        [ZERO, NEG, ZERO, ZERO],
        [ZERO, ZERO, NEG, ZERO],
        [ZERO, ZERO, ZERO, NEG],
    ]
    return MetricTensor(g, (t, x, y, z))
