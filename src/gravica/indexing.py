"""Index raising and lowering utilities for tensors."""

from __future__ import annotations


from gravica.metric import MetricTensor, Grid, ZERO
from gravica.simplify import simplify, str_is_zero


def raise_index_2d(metric: MetricTensor, T_lower: Grid, which: int) -> Grid:
    r"""Raise one index of a rank-2 covariant tensor.

    Parameters
    ----------
    metric : MetricTensor
    T_lower : Grid
        :math:`T_{ab}` as a 2D list.
    which : int
        Which index to raise: ``0`` gives :math:`T^a{}_b = g^{ac}\,T_{cb}`,
        ``1`` gives :math:`T_a{}^b = g^{bc}\,T_{ac}`.

    Returns
    -------
    Grid
        The mixed tensor as a 2D list.
    """
    n = metric.dim
    g_inv = metric.inverse
    result: Grid = [[ZERO for _ in range(n)] for _ in range(n)]

    for a in range(n):
        for b in range(n):
            val = ZERO
            if which == 0:
                # T^a_b = g^{ac} T_{cb}
                for c in range(n):
                    if str_is_zero(g_inv[a][c]):
                        continue
                    if str_is_zero(T_lower[c][b]):
                        continue
                    val = val + g_inv[a][c] * T_lower[c][b]
            elif which == 1:
                # T_a^b = g^{bc} T_{ac}
                for c in range(n):
                    if str_is_zero(g_inv[b][c]):
                        continue
                    if str_is_zero(T_lower[a][c]):
                        continue
                    val = val + g_inv[b][c] * T_lower[a][c]
            else:
                raise ValueError(f"which must be 0 or 1, got {which}")
            result[a][b] = simplify(val)

    return result


def lower_index_2d(metric: MetricTensor, T_upper: Grid, which: int) -> Grid:
    r"""Lower one index of a rank-2 contravariant tensor.

    Parameters
    ----------
    metric : MetricTensor
    T_upper : Grid
        :math:`T^{ab}` as a 2D list.
    which : int
        Which index to lower: ``0`` gives :math:`T_a{}^b = g_{ac}\,T^{cb}`,
        ``1`` gives :math:`T^a{}_b = g_{bc}\,T^{ac}`.

    Returns
    -------
    Grid
        The mixed tensor as a 2D list.
    """
    n = metric.dim
    g = metric.components
    result: Grid = [[ZERO for _ in range(n)] for _ in range(n)]

    for a in range(n):
        for b in range(n):
            val = ZERO
            if which == 0:
                # T_a^b = g_{ac} T^{cb}
                for c in range(n):
                    if str_is_zero(g[a][c]):
                        continue
                    if str_is_zero(T_upper[c][b]):
                        continue
                    val = val + g[a][c] * T_upper[c][b]
            elif which == 1:
                # T^a_b = g_{bc} T^{ac}
                for c in range(n):
                    if str_is_zero(g[b][c]):
                        continue
                    if str_is_zero(T_upper[a][c]):
                        continue
                    val = val + g[b][c] * T_upper[a][c]
            else:
                raise ValueError(f"which must be 0 or 1, got {which}")
            result[a][b] = simplify(val)

    return result
