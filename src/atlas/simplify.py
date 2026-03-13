"""Simplification utilities for symbolic GR expressions."""

from symbolica import Expression

_ZERO = Expression.num(0)


def str_is_zero(expr: Expression) -> bool:
    """Fast zero check via string representation (for optimization skips)."""
    return str(expr) == "0"


def simplify(expr: Expression) -> Expression:
    """Best-effort simplification: together then cancel.

    Using together() first puts everything over a common denominator,
    which allows cancel() to find and eliminate common polynomial factors.
    """
    try:
        result = expr.together().cancel()
    except Exception:
        try:
            result = expr.expand().cancel()
        except Exception:
            result = expr
    return result


def is_zero(expr: Expression, coords: tuple[Expression, ...] | None = None) -> bool:
    """Check if a symbolic expression is zero.

    First tries symbolic simplification. Falls back to random numerical
    substitution if that is inconclusive.
    """
    s = simplify(expr)
    text = str(s)
    if str_is_zero(s):
        return True

    if coords is not None:
        import random
        random.seed(42)
        for _ in range(3):
            vals = {c: Expression.num(random.randint(2, 97)) for c in coords}
            try:
                numerical = s
                for var, val in vals.items():
                    numerical = numerical.replace(var, val)
                numerical = simplify(numerical)
                if not str_is_zero(numerical):
                    return False
            except Exception:
                pass
        return True

    return False
