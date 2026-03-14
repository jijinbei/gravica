"""Gödel metric."""

from symbolica import Expression, S

from gravica.metric import MetricTensor, ZERO

ONE = Expression.num(1)
NEG = Expression.num(-1)
TWO = Expression.num(2)


def godel() -> MetricTensor:
    """Gödel metric ds² = a²[dt² - dx² + (e^{2x}/2)dy² - dz² + 2e^x dt dy].

    Coordinates: (t, x, y, z)
    Parameter: a (related to angular velocity ω by a² = 1/(2ω²))

    Uses exp as a symbolic function: exp(x) = e^x.
    """
    t, x, y, z = S('t'), S('x'), S('y'), S('z')
    a = S('a')
    exp = S('exp')

    ex = exp(x)
    e2x = exp(TWO * x)

    a2 = a**2

    g = [
        [a2, ZERO, a2 * ex, ZERO],
        [ZERO, NEG * a2, ZERO, ZERO],
        [a2 * ex, ZERO, a2 * e2x / TWO, ZERO],
        [ZERO, ZERO, ZERO, NEG * a2],
    ]
    return MetricTensor(g, (t, x, y, z))
