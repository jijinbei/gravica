"""Metric tensor and symbolic matrix utilities for GR."""

from __future__ import annotations

from symbolica import Expression

from atlas.simplify import simplify

# Type alias: a 4x4 (or NxN) grid of Expressions
Grid = list[list[Expression]]

ZERO = Expression.num(0)
ONE = Expression.num(1)


def _minor(matrix: Grid, row: int, col: int) -> Grid:
    """Return the (row, col) minor of a square matrix."""
    return [
        [matrix[i][j] for j in range(len(matrix)) if j != col]
        for i in range(len(matrix))
        if i != row
    ]


def symbolic_det(matrix: Grid) -> Expression:
    """Compute the determinant of a symbolic square matrix via Laplace expansion."""
    n = len(matrix)
    if n == 1:
        return matrix[0][0]
    if n == 2:
        return simplify(matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0])

    det = ZERO
    for j in range(n):
        if str(matrix[0][j]) == "0":
            continue
        cofactor = symbolic_det(_minor(matrix, 0, j))
        if j % 2 == 0:
            det = det + matrix[0][j] * cofactor
        else:
            det = det - matrix[0][j] * cofactor
    return simplify(det)


def symbolic_inverse(matrix: Grid) -> Grid:
    """Compute the inverse of a symbolic square matrix using cofactors."""
    n = len(matrix)
    det = symbolic_det(matrix)
    det_str = str(det)
    if det_str == "0":
        raise ValueError("Matrix is singular (determinant is zero)")

    inv_det = ONE / det

    cofactors: Grid = []
    for i in range(n):
        row = []
        for j in range(n):
            m = symbolic_det(_minor(matrix, i, j))
            sign = ONE if (i + j) % 2 == 0 else Expression.num(-1)
            row.append(simplify(sign * m * inv_det))
        cofactors.append(row)

    # Transpose to get adjugate
    result: Grid = [[cofactors[j][i] for j in range(n)] for i in range(n)]
    return result


class MetricTensor:
    """A metric tensor g_{ab} for a given coordinate system.

    Parameters
    ----------
    components : Grid
        NxN symmetric matrix of Expression entries.
    coords : tuple of Expression
        Coordinate symbols, e.g. (t, r, theta, phi).
    """

    def __init__(self, components: Grid, coords: tuple[Expression, ...]):
        n = len(coords)
        if len(components) != n or any(len(row) != n for row in components):
            raise ValueError(f"Components must be {n}x{n}, got {len(components)} rows")
        self.components = components
        self.coords = coords
        self.dim = n
        self._inverse: Grid | None = None
        self._det: Expression | None = None

    def __getitem__(self, idx: tuple[int, int]) -> Expression:
        """g_{ab} = metric[a, b]."""
        return self.components[idx[0]][idx[1]]

    @property
    def inverse(self) -> Grid:
        """Compute g^{ab} (cached)."""
        if self._inverse is None:
            self._inverse = symbolic_inverse(self.components)
        return self._inverse

    def inv(self, a: int, b: int) -> Expression:
        """g^{ab}."""
        return self.inverse[a][b]

    @property
    def det(self) -> Expression:
        """det(g_{ab}) (cached)."""
        if self._det is None:
            self._det = symbolic_det(self.components)
        return self._det

    def is_diagonal(self) -> bool:
        for i in range(self.dim):
            for j in range(self.dim):
                if i != j and str(self.components[i][j]) != "0":
                    return False
        return True
