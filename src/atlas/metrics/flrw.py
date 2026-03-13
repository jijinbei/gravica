"""Friedmann-Lemaître-Robertson-Walker metric."""

from symbolica import Expression, S

from atlas.metric import MetricTensor, ZERO

ONE = Expression.num(1)
NEG = Expression.num(-1)


def flrw() -> MetricTensor:
    """FLRW metric ds² = dt² - a(t)²[dr²/(1-kr²) + r²dΩ²].

    Coordinates: (t, r, θ, φ)
    Parameters: k (curvature), a(t) (scale factor as function)
    """
    t, r, theta, phi = S('t'), S('r'), S('theta'), S('phi')
    k = S('k')
    a_sym = S('a')
    sin = S('sin')

    a = a_sym(t)  # a(t) — scale factor as a function of t
    a2 = a**2

    g = [
        [ONE, ZERO, ZERO, ZERO],
        [ZERO, NEG * a2 / (ONE - k * r**2), ZERO, ZERO],
        [ZERO, ZERO, NEG * a2 * r**2, ZERO],
        [ZERO, ZERO, ZERO, NEG * a2 * r**2 * sin(theta)**2],
    ]
    return MetricTensor(g, (t, r, theta, phi))
