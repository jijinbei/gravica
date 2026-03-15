"""Kerr metric in Boyer-Lindquist coordinates."""

from symbolica import Expression, S

from gravica.metric import MetricTensor, ZERO
from gravica.metrics._symbols import _greek

ONE = Expression.num(1)
NEG = Expression.num(-1)
TWO = Expression.num(2)


def kerr() -> MetricTensor:
    r"""Kerr metric in Boyer--Lindquist coordinates :math:`(t, r, \theta, \varphi)`.

    Parameters: :math:`M` (mass), :math:`a = J/M` (spin).

    .. math::

        \Sigma = r^2 + a^2 \cos^2\theta, \qquad
        \Delta = r^2 - 2Mr + a^2
    """
    t, r = S("t"), S("r")
    theta = _greek("theta")
    phi = _greek("phi")
    M, a = S("M"), S("a")
    sin = S("sin")
    cos = S("cos")

    sin_th = sin(theta)
    cos_th = cos(theta)
    sin2 = sin_th**2
    cos2 = cos_th**2

    Sigma = r**2 + a**2 * cos2
    Delta = r**2 - TWO * M * r + a**2

    g_tt = ONE - TWO * M * r / Sigma
    g_tphi = TWO * M * a * r * sin2 / Sigma
    g_rr = NEG * Sigma / Delta
    g_thth = NEG * Sigma
    g_phiphi = NEG * (r**2 + a**2 + TWO * M * a**2 * r * sin2 / Sigma) * sin2

    g = [
        [g_tt, ZERO, ZERO, NEG * g_tphi],
        [ZERO, g_rr, ZERO, ZERO],
        [ZERO, ZERO, g_thth, ZERO],
        [NEG * g_tphi, ZERO, ZERO, g_phiphi],
    ]
    return MetricTensor(g, (t, r, theta, phi))
