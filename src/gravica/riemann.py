"""Riemann curvature tensor."""

from __future__ import annotations

from symbolica import Expression

from gravica.christoffel import ChristoffelSymbols
from gravica.metric import ZERO
from gravica.simplify import simplify, str_is_zero


class RiemannTensor:
    """Riemann curvature tensor R^a_{bcd}.

    R^a_{bcd} = ∂_c Γ^a_{db} - ∂_d Γ^a_{cb} + Γ^a_{ce} Γ^e_{db} - Γ^a_{de} Γ^e_{cb}
    """

    def __init__(self, christoffel: ChristoffelSymbols):
        self.christoffel = christoffel
        self.metric = christoffel.metric
        self.dim = christoffel.dim
        self.coords = christoffel.coords
        self._components: list[list[list[list[Expression]]]] | None = None

    @property
    def components(self) -> list[list[list[list[Expression]]]]:
        """R^a_{bcd} indexed as [a][b][c][d]."""
        if self._components is None:
            self._compute()
        return self._components  # type: ignore

    def _compute(self) -> None:
        gamma = self.christoffel.second_kind
        coords = self.coords
        n = self.dim

        self._components = [
            [[[ZERO for _ in range(n)] for _ in range(n)] for _ in range(n)]
            for _ in range(n)
        ]

        for a in range(n):
            for b in range(n):
                for c in range(n):
                    for d in range(c + 1, n):  # Antisymmetric in c,d
                        # ∂_c Γ^a_{db} - ∂_d Γ^a_{cb}
                        val = (
                            gamma[a][d][b].derivative(coords[c])
                            - gamma[a][c][b].derivative(coords[d])
                        )

                        # + Γ^a_{ce} Γ^e_{db} - Γ^a_{de} Γ^e_{cb}
                        for e in range(n):
                            g_ace = gamma[a][c][e]
                            g_ade = gamma[a][d][e]
                            g_edb = gamma[e][d][b]
                            g_ecb = gamma[e][c][b]

                            if not str_is_zero(g_ace) and not str_is_zero(g_edb):
                                val = val + g_ace * g_edb
                            if not str_is_zero(g_ade) and not str_is_zero(g_ecb):
                                val = val - g_ade * g_ecb

                        val = simplify(val)
                        self._components[a][b][c][d] = val
                        # Antisymmetric: R^a_{bdc} = -R^a_{bcd}
                        self._components[a][b][d][c] = simplify(Expression.num(-1) * val)

    def __getitem__(self, idx: tuple[int, int, int, int]) -> Expression:
        """R^a_{bcd} = riemann[a, b, c, d]."""
        return self.components[idx[0]][idx[1]][idx[2]][idx[3]]

    def fully_covariant(self, a: int, b: int, c: int, d: int) -> Expression:
        """R_{abcd} = g_{ae} R^e_{bcd}."""
        g = self.metric.components
        n = self.dim
        val = ZERO
        for e in range(n):
            if str_is_zero(g[a][e]):
                continue
            r_ebcd = self[e, b, c, d]
            if str_is_zero(r_ebcd):
                continue
            val = val + g[a][e] * r_ebcd
        return simplify(val)
