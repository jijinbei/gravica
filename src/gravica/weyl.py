"""Weyl conformal tensor."""

from __future__ import annotations

from symbolica import Expression

from gravica.riemann import RiemannTensor
from gravica.ricci import RicciTensor, ricci_scalar
from gravica.metric import ZERO
from gravica.simplify import simplify


class WeylTensor:
    r"""Weyl tensor :math:`C_{abcd}` (fully covariant).

    .. math::

        C_{abcd} = R_{abcd}
            - \frac{2}{n-2}\bigl(g_{a[c}\,R_{d]b} - g_{b[c}\,R_{d]a}\bigr)
            + \frac{2}{(n-1)(n-2)}\,R\;g_{a[c}\,g_{d]b}
    """

    def __init__(self, riemann: RiemannTensor, ricci: RicciTensor):
        self.riemann = riemann
        self.ricci = ricci
        self.metric = riemann.metric
        self.dim = riemann.dim
        self._components: list[list[list[list[Expression]]]] | None = None

    @property
    def components(self) -> list[list[list[list[Expression]]]]:
        r""":math:`C_{abcd}` indexed as ``[a][b][c][d]``."""
        if self._components is None:
            self._compute()
        return self._components  # type: ignore

    def _compute(self) -> None:
        g = self.metric.components
        R = self.riemann
        Ric = self.ricci
        R_scalar = ricci_scalar(Ric)
        n = self.dim

        if n < 3:
            raise ValueError("Weyl tensor is only defined for dimension >= 3")

        coeff1 = Expression.num(1) / Expression.num(n - 2)
        coeff2 = Expression.num(1) / (Expression.num(n - 1) * Expression.num(n - 2))

        self._components = [
            [[[ZERO for _ in range(n)] for _ in range(n)] for _ in range(n)]
            for _ in range(n)
        ]

        for a in range(n):
            for b in range(n):
                for c in range(n):
                    for d in range(n):
                        R_abcd = R.fully_covariant(a, b, c, d)

                        # -1/(n-2) * (g_{ac} R_{db} - g_{ad} R_{cb}
                        #           - g_{bc} R_{da} + g_{bd} R_{ca})
                        ricci_part = (
                            g[a][c] * Ric[d, b]
                            - g[a][d] * Ric[c, b]
                            - g[b][c] * Ric[d, a]
                            + g[b][d] * Ric[c, a]
                        )

                        # +1/((n-1)(n-2)) * R * (g_{ac} g_{db} - g_{ad} g_{cb})
                        scalar_part = R_scalar * (g[a][c] * g[d][b] - g[a][d] * g[c][b])

                        val = R_abcd - coeff1 * ricci_part + coeff2 * scalar_part
                        self._components[a][b][c][d] = simplify(val)

    def __getitem__(self, idx: tuple[int, int, int, int]) -> Expression:
        r""":math:`C_{abcd}` = ``weyl[a, b, c, d]``."""
        return self.components[idx[0]][idx[1]][idx[2]][idx[3]]
