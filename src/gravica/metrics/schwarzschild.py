"""Schwarzschild metric in Schwarzschild coordinates."""

from symbolica import Expression, S

from gravica.metric import MetricTensor, ZERO

ONE = Expression.num(1)
NEG = Expression.num(-1)


def schwarzschild() -> MetricTensor:
    """Schwarzschild metric ds² = (1-r_s/r)dt² - (1-r_s/r)^{-1}dr² - r²dΩ².

    Uses rational form (r-r_s)/r for better symbolic cancellation.
    Coordinates: (t, r, θ, φ)
    Parameter: r_s = 2GM/c²
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
