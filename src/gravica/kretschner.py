"""Kretschner scalar — a coordinate-invariant curvature invariant."""

from __future__ import annotations

from symbolica import Expression

from gravica.riemann import RiemannTensor
from gravica.metric import ZERO
from gravica.simplify import simplify, str_is_zero


def kretschner_scalar(riemann: RiemannTensor) -> Expression:
    r"""Kretschner scalar :math:`K = R_{abcd}\,R^{abcd}`.

    Detects true curvature singularities (coordinate-independent).
    For Schwarzschild: :math:`K = 12\,r_s^2 / r^6`.
    """
    n = riemann.dim
    g_inv = riemann.metric.inverse

    # Precompute R_{abcd} (fully covariant)
    R_cov = [
        [[[riemann.fully_covariant(a, b, c, d) for d in range(n)]
          for c in range(n)]
         for b in range(n)]
        for a in range(n)
    ]

    # K = R_{abcd} R^{abcd}
    #   = R_{abcd} g^{ae} g^{bf} g^{cg} g^{dh} R_{efgh}
    # Precompute R^{abcd} by raising all 4 indices
    val = ZERO

    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    r_abcd = R_cov[a][b][c][d]
                    if str_is_zero(r_abcd):
                        continue

                    # Compute R^{abcd}
                    r_upper = ZERO
                    for e in range(n):
                        g_ae = g_inv[a][e]
                        if str_is_zero(g_ae):
                            continue
                        for f in range(n):
                            g_bf = g_inv[b][f]
                            if str_is_zero(g_bf):
                                continue
                            for g in range(n):
                                g_cg = g_inv[c][g]
                                if str_is_zero(g_cg):
                                    continue
                                for h in range(n):
                                    g_dh = g_inv[d][h]
                                    if str_is_zero(g_dh):
                                        continue
                                    r_efgh = R_cov[e][f][g][h]
                                    if str_is_zero(r_efgh):
                                        continue
                                    r_upper = r_upper + (
                                        g_ae * g_bf * g_cg * g_dh * r_efgh
                                    )

                    if not str_is_zero(r_upper):
                        val = val + r_abcd * r_upper

    return simplify(val)
