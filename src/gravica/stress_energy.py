"""Stress-energy-momentum tensor from the Einstein field equations."""

from __future__ import annotations

from symbolica import Expression

from gravica.einstein import EinsteinTensor
from gravica.metric import ZERO
from gravica.simplify import simplify


class StressEnergyTensor:
    r"""Stress-energy-momentum tensor from the Einstein field equations.

    .. math::

        T_{ab} = \frac{1}{8\pi G}\bigl(G_{ab} + \Lambda\,g_{ab}\bigr)

    When :math:`\Lambda = 0` (default), simplifies to
    :math:`T_{ab} = G_{ab} / (8\pi G)`.

    The result is expressed with an overall factor of :math:`1/(8\pi G)`,
    represented by the symbol ``_8piG_inv``.  For vacuum solutions the
    components are exactly zero regardless of this prefactor.
    """

    def __init__(
        self,
        einstein: EinsteinTensor,
        cosmological_constant: Expression | None = None,
    ):
        self.einstein = einstein
        self.metric = einstein.metric
        self.dim = einstein.dim
        self.cosmological_constant = cosmological_constant
        self._components: list[list[Expression]] | None = None

    @property
    def components(self) -> list[list[Expression]]:
        r""":math:`8\pi G\,T_{ab}` indexed as ``[a][b]``."""
        if self._components is None:
            self._compute()
        return self._components  # type: ignore

    def _compute(self) -> None:
        G = self.einstein
        g = self.metric.components
        n = self.dim
        Lambda = self.cosmological_constant

        self._components = [[ZERO for _ in range(n)] for _ in range(n)]

        for a in range(n):
            for b in range(a, n):  # Symmetric
                val = G[a, b]
                if Lambda is not None:
                    val = val + Lambda * g[a][b]
                val = simplify(val)
                self._components[a][b] = val
                self._components[b][a] = val

    def __getitem__(self, idx: tuple[int, int]) -> Expression:
        r""":math:`8\pi G\,T_{ab}` = ``stress_energy[a, b]``."""
        return self.components[idx[0]][idx[1]]
