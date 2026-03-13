"""Tests for Christoffel symbols."""

from atlas.simplify import simplify


def test_minkowski_all_zero(mink_christoffel):
    """All Christoffel symbols vanish for flat spacetime."""
    ch = mink_christoffel
    for a in range(4):
        for b in range(4):
            for c in range(4):
                assert str(ch[a, b, c]) == "0", f"Γ^{a}_{{{b}{c}}} != 0"


def test_schwarzschild_symmetry(schw_christoffel):
    """Γ^a_{bc} = Γ^a_{cb} (symmetric in lower indices)."""
    ch = schw_christoffel
    for a in range(4):
        for b in range(4):
            for c in range(b + 1, 4):
                diff = simplify(ch[a, b, c] - ch[a, c, b])
                assert str(diff) == "0", f"Γ^{a}_{{{b}{c}}} != Γ^{a}_{{{c}{b}}}"


def test_schwarzschild_known_values(schw_christoffel):
    """Check known Schwarzschild Christoffel symbols."""
    ch = schw_christoffel
    from symbolica import S, Expression
    r, r_s = S('r'), S('r_s')

    # Γ^r_{tt} = r_s(r - r_s) / (2r³)
    gamma_r_tt = ch[1, 0, 0]
    expected = r_s * (r - r_s) / (Expression.num(2) * r**3)
    diff = simplify(gamma_r_tt - expected)
    assert str(diff) == "0", f"Γ^r_tt = {gamma_r_tt}, expected {expected}"

    # Γ^θ_{rθ} = 1/r
    gamma_th_rth = ch[2, 1, 2]
    diff = simplify(gamma_th_rth - Expression.num(1) / r)
    assert str(diff) == "0", f"Γ^θ_rθ = {gamma_th_rth}"
