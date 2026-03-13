"""Christoffel symbols of the first and second kind."""

from __future__ import annotations

from symbolica import Expression

from atlas.metric import MetricTensor, Grid, ZERO
from atlas.simplify import simplify


class ChristoffelSymbols:
    """Christoffel symbols computed from a metric tensor.

    Γ_{abc} = ½(∂_b g_{ac} + ∂_c g_{ab} - ∂_a g_{bc})   (first kind)
    Γ^a_{bc} = g^{ad} Γ_{dbc}                              (second kind)
    """

    def __init__(self, metric: MetricTensor):
        self.metric = metric
        self.dim = metric.dim
        self.coords = metric.coords
        self._first: list[list[list[Expression]]] | None = None
        self._second: list[list[list[Expression]]] | None = None

    @property
    def first_kind(self) -> list[list[list[Expression]]]:
        """Γ_{abc} indexed as [a][b][c]."""
        if self._first is None:
            self._compute_first_kind()
        return self._first  # type: ignore

    @property
    def second_kind(self) -> list[list[list[Expression]]]:
        """Γ^a_{bc} indexed as [a][b][c]."""
        if self._second is None:
            self._compute_second_kind()
        return self._second  # type: ignore

    def _compute_first_kind(self) -> None:
        g = self.metric.components
        coords = self.coords
        n = self.dim
        half = Expression.num(1) / Expression.num(2)

        self._first = [[[ZERO for _ in range(n)] for _ in range(n)] for _ in range(n)]

        for a in range(n):
            for b in range(n):
                for c in range(b, n):  # Symmetric in b,c
                    # Γ_{abc} = ½(∂_b g_{ac} + ∂_c g_{ab} - ∂_a g_{bc})
                    term = (
                        g[a][c].derivative(coords[b])
                        + g[a][b].derivative(coords[c])
                        - g[b][c].derivative(coords[a])
                    )
                    val = simplify(half * term)
                    self._first[a][b][c] = val
                    self._first[a][c][b] = val  # Symmetric in b,c

    def _compute_second_kind(self) -> None:
        first = self.first_kind
        g_inv = self.metric.inverse
        n = self.dim

        self._second = [[[ZERO for _ in range(n)] for _ in range(n)] for _ in range(n)]

        for a in range(n):
            for b in range(n):
                for c in range(b, n):  # Symmetric in b,c
                    # Γ^a_{bc} = g^{ad} Γ_{dbc}
                    val = ZERO
                    for d in range(n):
                        g_ad = g_inv[a][d]
                        if str(g_ad) == "0":
                            continue
                        gamma_dbc = first[d][b][c]
                        if str(gamma_dbc) == "0":
                            continue
                        val = val + g_ad * gamma_dbc
                    val = simplify(val)
                    self._second[a][b][c] = val
                    self._second[a][c][b] = val  # Symmetric in b,c

    def __getitem__(self, idx: tuple[int, int, int]) -> Expression:
        """Γ^a_{bc} = christoffel[a, b, c]."""
        return self.second_kind[idx[0]][idx[1]][idx[2]]
