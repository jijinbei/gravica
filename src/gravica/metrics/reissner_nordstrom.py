"""Reissner-Nordström metric for a charged, non-rotating black hole."""

from symbolica import Expression, S

from gravica.metric import MetricTensor, ZERO
from gravica.metrics._symbols import _greek

NEG = Expression.num(-1)


def reissner_nordstrom() -> MetricTensor:
    r"""Reissner-Nordström metric in coordinates :math:`(t, r, \theta, \varphi)`.

    .. math::

        ds^2 = f(r)\,dt^2 - f(r)^{-1}\,dr^2 - r^2\,d\Omega^2

    where :math:`f(r) = 1 - r_s/r + r_Q^2/r^2 = (r^2 - r_s\,r + r_Q^2)/r^2`.

    Parameters: :math:`r_s = 2GM/c^2`, :math:`r_Q^2 = GQ^2/(4\pi\varepsilon_0 c^4)`.
    """
    t, r = S("t"), S("r")
    theta = _greek("theta")
    phi = _greek("phi")
    r_s = S("r_s")
    r_Q = S("r_Q")
    sin = S("sin")

    # f(r) = (r² - r_s r + r_Q²) / r²
    f = (r**2 - r_s * r + r_Q**2) / r**2

    g = [
        [f, ZERO, ZERO, ZERO],
        [ZERO, NEG / f, ZERO, ZERO],
        [ZERO, ZERO, NEG * r**2, ZERO],
        [ZERO, ZERO, ZERO, NEG * r**2 * sin(theta) ** 2],
    ]
    return MetricTensor(g, (t, r, theta, phi))
