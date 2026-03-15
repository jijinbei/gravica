"""Ricci tensor and Ricci scalar."""

from __future__ import annotations

from symbolica import Expression

from gravica.christoffel import ChristoffelSymbols
from gravica.riemann import RiemannTensor
from gravica.metric import MetricTensor, ZERO
from gravica.simplify import simplify, str_is_zero


class RicciTensor:
    r"""Ricci tensor :math:`R_{ab} = R^c_{\ acb}`."""

    def __init__(self, riemann: RiemannTensor):
        self.riemann = riemann
        self.metric = riemann.metric
        self.dim = riemann.dim
        self._components: list[list[Expression]] | None = None

    @classmethod
    def from_metric(cls, metric: MetricTensor) -> RicciTensor:
        """Build from a :class:`~gravica.metric.MetricTensor`."""
        return cls(RiemannTensor(ChristoffelSymbols(metric)))

    @property
    def components(self) -> list[list[Expression]]:
        r""":math:`R_{ab}` indexed as ``[a][b]``."""
        if self._components is None:
            self._compute()
        return self._components  # type: ignore

    def _compute(self) -> None:
        R = self.riemann
        n = self.dim
        self._components = [[ZERO for _ in range(n)] for _ in range(n)]

        for a in range(n):
            for b in range(a, n):  # Symmetric
                val = ZERO
                for c in range(n):
                    r_cacb = R[c, a, c, b]
                    if not str_is_zero(r_cacb):
                        val = val + r_cacb
                val = simplify(val)
                self._components[a][b] = val
                self._components[b][a] = val

    def __getitem__(self, idx: tuple[int, int]) -> Expression:
        r""":math:`R_{ab}` = ``ricci[a, b]``."""
        return self.components[idx[0]][idx[1]]


def ricci_scalar(ricci: RicciTensor) -> Expression:
    r"""Ricci scalar :math:`R = g^{ab}\,R_{ab}`."""
    g_inv = ricci.metric.inverse
    n = ricci.dim
    val = ZERO
    for a in range(n):
        for b in range(n):
            g_ab = g_inv[a][b]
            r_ab = ricci[a, b]
            if not str_is_zero(g_ab) and not str_is_zero(r_ab):
                val = val + g_ab * r_ab
    return simplify(val)
