"""Anti-de Sitter metric — spacetime with negative cosmological constant."""

from symbolica import Expression, S

from gravica.metric import MetricTensor, ZERO, ONE

NEG = Expression.num(-1)


def anti_de_sitter() -> MetricTensor:
    r"""Anti-de Sitter metric in static coordinates :math:`(t, r, \theta, \varphi)`.

    .. math::

        ds^2 = \Bigl(1 + \frac{r^2}{l^2}\Bigr)\,dt^2
             - \Bigl(1 + \frac{r^2}{l^2}\Bigr)^{-1}\,dr^2
             - r^2\,d\Omega^2

    Parameter: :math:`l` (AdS radius, related to cosmological constant
    :math:`\Lambda = -3/l^2`).
    """
    t, r, theta, phi = S("t"), S("r"), S("theta"), S("phi")
    ads_l = S("l")
    sin = S("sin")

    f = ONE + r**2 / ads_l**2

    g = [
        [f, ZERO, ZERO, ZERO],
        [ZERO, NEG / f, ZERO, ZERO],
        [ZERO, ZERO, NEG * r**2, ZERO],
        [ZERO, ZERO, ZERO, NEG * r**2 * sin(theta) ** 2],
    ]
    return MetricTensor(g, (t, r, theta, phi))
