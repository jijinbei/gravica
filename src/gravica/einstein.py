"""Einstein tensor."""

from __future__ import annotations

from symbolica import Expression

from gravica.christoffel import ChristoffelSymbols
from gravica.riemann import RiemannTensor
from gravica.ricci import RicciTensor, ricci_scalar
from gravica.metric import MetricTensor, ZERO
from gravica.simplify import simplify


class EinsteinTensor:
    r"""Einstein tensor :math:`G_{ab} = R_{ab} - \tfrac{1}{2}\,g_{ab}\,R`."""

    def __init__(self, ricci: RicciTensor):
        self.ricci = ricci
        self.metric = ricci.metric
        self.dim = ricci.dim
        self._components: list[list[Expression]] | None = None
        self._scalar: Expression | None = None

    @classmethod
    def from_metric(cls, metric: MetricTensor) -> EinsteinTensor:
        """Build from a :class:`~gravica.metric.MetricTensor`."""
        return cls(RicciTensor(RiemannTensor(ChristoffelSymbols(metric))))

    @property
    def scalar(self) -> Expression:
        if self._scalar is None:
            self._scalar = ricci_scalar(self.ricci)
        return self._scalar

    @property
    def components(self) -> list[list[Expression]]:
        r""":math:`G_{ab}` indexed as ``[a][b]``."""
        if self._components is None:
            self._compute()
        return self._components  # type: ignore

    def _compute(self) -> None:
        g = self.metric.components
        R = self.ricci
        R_scalar = self.scalar
        n = self.dim
        half = Expression.num(1) / Expression.num(2)

        self._components = [[ZERO for _ in range(n)] for _ in range(n)]

        for a in range(n):
            for b in range(a, n):  # Symmetric
                val = simplify(R[a, b] - half * g[a][b] * R_scalar)
                self._components[a][b] = val
                self._components[b][a] = val

    def __getitem__(self, idx: tuple[int, int]) -> Expression:
        r""":math:`G_{ab}` = ``einstein[a, b]``."""
        return self.components[idx[0]][idx[1]]
