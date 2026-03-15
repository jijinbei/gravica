"""de Sitter metric — spacetime with positive cosmological constant."""

from symbolica import Expression, S

from gravica.metric import MetricTensor, ZERO, ONE

NEG = Expression.num(-1)


def de_sitter() -> MetricTensor:
    r"""de Sitter metric in static coordinates :math:`(t, r, \theta, \varphi)`.

    .. math::

        ds^2 = \Bigl(1 - \frac{\Lambda r^2}{3}\Bigr)\,dt^2
             - \Bigl(1 - \frac{\Lambda r^2}{3}\Bigr)^{-1}\,dr^2
             - r^2\,d\Omega^2

    Parameter: :math:`\Lambda > 0` (cosmological constant).
    """
    t, r, theta, phi = S("t"), S("r"), S("theta"), S("phi")
    Lambda = S("Lambda")
    sin = S("sin")

    f = ONE - Lambda * r**2 / Expression.num(3)

    g = [
        [f, ZERO, ZERO, ZERO],
        [ZERO, NEG / f, ZERO, ZERO],
        [ZERO, ZERO, NEG * r**2, ZERO],
        [ZERO, ZERO, ZERO, NEG * r**2 * sin(theta) ** 2],
    ]
    return MetricTensor(g, (t, r, theta, phi))
