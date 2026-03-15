"""Schwarzschild metric in Schwarzschild coordinates."""

from symbolica import Expression, S

from gravica.metric import MetricTensor, ZERO

ONE = Expression.num(1)
NEG = Expression.num(-1)


def schwarzschild() -> MetricTensor:
    r"""Schwarzschild metric in coordinates :math:`(t, r, \theta, \varphi)`.

    .. math::

        ds^2 = \Bigl(1 - \frac{r_s}{r}\Bigr) dt^2
             - \Bigl(1 - \frac{r_s}{r}\Bigr)^{-1} dr^2
             - r^2\,d\Omega^2

    Uses the rational form :math:`(r - r_s)/r` for better symbolic cancellation.
    Parameter: :math:`r_s = 2GM/c^2`.
    """
    t, r, theta, phi = S('t'), S('r'), S('theta'), S('phi')
    r_s = S('r_s')
    sin = S('sin')

    # Use (r - r_s)/r form for clean polynomial cancellation
    f = (r - r_s) / r  # = 1 - r_s/r

    g = [
        [f, ZERO, ZERO, ZERO],
        [ZERO, NEG * r / (r - r_s), ZERO, ZERO],
        [ZERO, ZERO, NEG * r**2, ZERO],
        [ZERO, ZERO, ZERO, NEG * r**2 * sin(theta)**2],
    ]
    return MetricTensor(g, (t, r, theta, phi))
