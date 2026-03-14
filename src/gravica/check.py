"""Validation helpers for tensor computations."""

from gravica.simplify import str_is_zero


def zero(tensor, simplify_fn):
    """Check if all components of a rank-2 tensor are zero after simplification."""
    n = tensor.dim
    for a in range(n):
        for b in range(n):
            val = simplify_fn(tensor[a, b])
            if not str_is_zero(val):
                return False
    return True


def antisymmetry(tensor, simplify_fn):
    """Check T^a_{bcd} + T^a_{bdc} = 0 for all index combinations of a rank-4 tensor."""
    n = tensor.dim
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    s = simplify_fn(tensor[a, b, c, d] + tensor[a, b, d, c])
                    if not str_is_zero(s):
                        return False
    return True
