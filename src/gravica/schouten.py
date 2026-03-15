"""Schouten tensor — used in conformal geometry."""

from __future__ import annotations

from symbolica import Expression

from gravica.ricci import RicciTensor, ricci_scalar
from gravica.metric import ZERO
from gravica.simplify import simplify


class SchoutenTensor:
    r"""Schouten tensor :math:`S_{ab} = \frac{1}{n-2}\bigl(R_{ab} - \frac{R\,g_{ab}}{2(n-1)}\bigr)`.

    Used in conformal geometry; the Weyl tensor can be expressed in terms
    of the Schouten tensor and the metric.
    """

    def __init__(self, ricci: RicciTensor):
        self.ricci = ricci
        self.metric = ricci.metric
        self.dim = ricci.dim
        self._components: list[list[Expression]] | None = None
        self._scalar: Expression | None = None

    @property
    def scalar(self) -> Expression:
        if self._scalar is None:
            self._scalar = ricci_scalar(self.ricci)
        return self._scalar

    @property
    def components(self) -> list[list[Expression]]:
        r""":math:`S_{ab}` indexed as ``[a][b]``."""
        if self._components is None:
            self._compute()
        return self._components  # type: ignore

    def _compute(self) -> None:
        g = self.metric.components
        R = self.ricci
        R_scalar = self.scalar
        n = self.dim

        if n < 3:
            raise ValueError("Schouten tensor is only defined for dimension >= 3")

        coeff = Expression.num(1) / Expression.num(n - 2)
        scalar_coeff = Expression.num(1) / (Expression.num(2) * Expression.num(n - 1))

        self._components = [[ZERO for _ in range(n)] for _ in range(n)]

        for a in range(n):
            for b in range(a, n):  # Symmetric
                val = simplify(coeff * (R[a, b] - scalar_coeff * R_scalar * g[a][b]))
                self._components[a][b] = val
                self._components[b][a] = val

    def __getitem__(self, idx: tuple[int, int]) -> Expression:
        r""":math:`S_{ab}` = ``schouten[a, b]``."""
        return self.components[idx[0]][idx[1]]
