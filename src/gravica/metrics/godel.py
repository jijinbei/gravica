"""Gödel metric."""

from symbolica import Expression, S

from gravica.metric import MetricTensor, ZERO

ONE = Expression.num(1)
NEG = Expression.num(-1)
TWO = Expression.num(2)


def godel() -> MetricTensor:
    r"""Gödel metric in coordinates :math:`(t, x, y, z)`.

    .. math::

        ds^2 = a^2 \left[
            dt^2 - dx^2 + \tfrac{1}{2}\,e^{2x}\,dy^2 - dz^2
            + 2\,e^{x}\,dt\,dy
        \right]

    Parameter: :math:`a` (related to angular velocity :math:`\omega` by
    :math:`a^2 = 1/(2\omega^2)`).

    Uses ``exp`` as a symbolic function: :math:`\exp(x) = e^x`.
    """
    t, x, y, z = S("t"), S("x"), S("y"), S("z")
    a = S("a")
    exp = S("exp")

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
