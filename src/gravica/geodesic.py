"""Geodesic equations — symbolic generation of equations of motion."""

from __future__ import annotations

from collections import namedtuple

from symbolica import Expression, S

from gravica.christoffel import ChristoffelSymbols
from gravica.simplify import str_is_zero


GeodesicEquations = namedtuple(
    "GeodesicEquations", ["equations", "velocities", "accelerations"]
)


def geodesic_equations(christoffel: ChristoffelSymbols) -> GeodesicEquations:
    r"""Generate the geodesic equations for the given spacetime.

    .. math::

        \frac{d^2 x^a}{d\tau^2}
        + \Gamma^a_{\ bc}\,\frac{dx^b}{d\tau}\,\frac{dx^c}{d\tau} = 0

    Returns the **left-hand side** for each coordinate as a
    :class:`GeodesicEquations` namedtuple ``(equations, velocities, accelerations)``.

    Velocity symbols are named ``d<coord>_`` (e.g. ``dt_``, ``dr_``) and
    acceleration symbols ``dd<coord>_`` (e.g. ``ddt_``, ``ddr_``).
    """
    coords = christoffel.coords
    n = christoffel.dim

    # Create velocity and acceleration symbols
    coord_names = [str(c) for c in coords]
    velocities = tuple(S(f"d{name}_") for name in coord_names)
    accelerations = tuple(S(f"dd{name}_") for name in coord_names)

    equations = []
    for a in range(n):
        # d²x^a/dτ² term
        eq = accelerations[a]

        # + Γ^a_{bc} (dx^b/dτ)(dx^c/dτ)
        for b in range(n):
            for c in range(b, n):  # Symmetric in b,c
                gamma = christoffel[a, b, c]
                if str_is_zero(gamma):
                    continue
                if b == c:
                    eq = eq + gamma * velocities[b] * velocities[c]
                else:
                    # Factor of 2 from symmetry
                    eq = eq + Expression.num(2) * gamma * velocities[b] * velocities[c]

        equations.append(eq)

    return GeodesicEquations(
        equations=tuple(equations),
        velocities=velocities,
        accelerations=accelerations,
    )
