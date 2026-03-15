"""Shared symbol constructors with LaTeX custom_print for metric definitions."""

from symbolica import S, PrintMode

_cache = {}


def _greek(name):
    """Create a symbol that renders as a LaTeX Greek letter.

    Cached so that repeated calls with the same name return the same
    symbol object, avoiding Symbolica's restriction on duplicate
    custom_print registrations.
    """
    if name in _cache:
        return _cache[name]

    tex = f"\\{name}"

    def _print(expr, mode, **kw):
        if mode == PrintMode.Latex:
            return tex

    sym = S(name, custom_print=_print)
    _cache[name] = sym
    return sym
