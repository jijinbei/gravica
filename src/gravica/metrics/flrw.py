"""Friedmann-Lemaître-Robertson-Walker metric."""

from symbolica import Expression, S

from gravica.metric import MetricTensor, ZERO

ONE = Expression.num(1)
NEG = Expression.num(-1)


def flrw() -> MetricTensor:
    r"""FLRW metric in coordinates :math:`(t, r, \theta, \varphi)`.

    .. math::

        ds^2 = dt^2 - a(t)^2 \left[
            \frac{dr^2}{1 - k\,r^2} + r^2\,d\Omega^2
        \right]

    Parameters: :math:`k` (spatial curvature), :math:`a(t)` (scale factor).
    """
    t, r, theta, phi = S("t"), S("r"), S("theta"), S("phi")
    k = S("k")
    a_sym = S("a")
    sin = S("sin")

    a = a_sym(t)  # a(t) — scale factor as a function of t
    a2 = a**2

    g = [
        [ONE, ZERO, ZERO, ZERO],
        [ZERO, NEG * a2 / (ONE - k * r**2), ZERO, ZERO],
        [ZERO, ZERO, NEG * a2 * r**2, ZERO],
        [ZERO, ZERO, ZERO, NEG * a2 * r**2 * sin(theta) ** 2],
    ]
    return MetricTensor(g, (t, r, theta, phi))
