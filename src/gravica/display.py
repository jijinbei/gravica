"""Display helpers for Jupyter notebooks — turn tensor objects into readable tables."""

from gravica.simplify import str_is_zero


def _fix_der_latex(s):
    r"""Replace Symbolica's ``der\!\left(N, FUNC\right)`` with dot notation.

    - ``der(1, f)`` → ``\dot{f}``
    - ``der(2, f)`` → ``\ddot{f}``
    """
    result = s
    while True:
        idx = result.find(r"der\!\left(")
        if idx == -1:
            break
        start = idx
        pos = idx + len(r"der\!\left(")

        comma = result.index(",", pos)
        order = result[pos:comma].strip()

        pos = comma + 1
        depth = 1
        while depth > 0 and pos < len(result):
            if result[pos:].startswith(r"\left("):
                depth += 1
                pos += len(r"\left(")
            elif result[pos:].startswith(r"\right)"):
                depth -= 1
                if depth == 0:
                    break
                pos += len(r"\right)")
            else:
                pos += 1

        end = pos + len(r"\right)")
        func_latex = result[comma + 1 : pos].strip()

        if order == "1":
            replacement = r"\dot{" + func_latex + "}"
        elif order == "2":
            replacement = r"\ddot{" + func_latex + "}"
        else:
            replacement = (
                r"\frac{d^{" + order + r"}}{dt^{" + order + r"}}" + func_latex
            )

        result = result[:start] + replacement + result[end:]

    return result


def _nonzero_components_2d(tensor, coord_names, symmetric=False):
    """Return list of (label, value) for non-zero rank-2 tensor components.

    If symmetric=True, only yields (a, b) with a <= b.
    """
    dim = len(coord_names)
    items = []
    for a in range(dim):
        b_start = a if symmetric else 0
        for b in range(b_start, dim):
            val = tensor[a, b]
            if not str_is_zero(val):
                label = f"{coord_names[a]} {coord_names[b]}"
                items.append((label, val))
    return items


def _nonzero_components_3d(christoffel, coord_names):
    """Return list of (label, value) for non-zero Christoffel symbol components.

    Uses b <= c symmetry.
    """
    dim = len(coord_names)
    items = []
    for a in range(dim):
        for b in range(dim):
            for c in range(b, dim):
                val = christoffel[a, b, c]
                if not str_is_zero(val):
                    label = (
                        f"^{{{coord_names[a]}}}_{{{coord_names[b]} {coord_names[c]}}}"
                    )
                    items.append((label, val))
    return items


def _nonzero_components_4d(riemann, coord_names):
    """Return list of (label, value) for non-zero Riemann tensor components.

    Uses c < d antisymmetry.
    """
    dim = len(coord_names)
    items = []
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                for d in range(c + 1, dim):
                    val = riemann[a, b, c, d]
                    if not str_is_zero(val):
                        label = f"^{{{coord_names[a]}}}{{}}_{{{coord_names[b]} {coord_names[c]} {coord_names[d]}}}"
                        items.append((label, val))
    return items


# --- Type-to-config mapping for components_table ---

_TENSOR_TABLE_CONFIG = {}


def _register_table_config(cls_name, symbol, style):
    _TENSOR_TABLE_CONFIG[cls_name] = (symbol, style)


_register_table_config("MetricTensor", "g", "sub")
_register_table_config("ChristoffelSymbols", "\\Gamma", "christoffel")
_register_table_config("RiemannTensor", "R", "riemann")
_register_table_config("RicciTensor", "R", "sub")
_register_table_config("EinsteinTensor", "G", "sub")
_register_table_config("WeylTensor", "C", "riemann")
_register_table_config("SchoutenTensor", "S", "sub")
_register_table_config("StressEnergyTensor", "T", "sub")


# --- Unified public API ---


def nonzero_components(tensor, coord_names):
    r"""Return list of (label, value) for non-zero components, auto-dispatched by tensor type.

    Dispatches based on tensor class:

    - :class:`~gravica.metric.MetricTensor` -- rank-2, symmetric
    - :class:`~gravica.ricci.RicciTensor`, :class:`~gravica.einstein.EinsteinTensor` -- rank-2, non-symmetric
    - :class:`~gravica.christoffel.ChristoffelSymbols` -- rank-3
    - :class:`~gravica.riemann.RiemannTensor`, :class:`~gravica.weyl.WeylTensor` -- rank-4
    """
    cls_name = type(tensor).__name__
    if cls_name == "MetricTensor":
        return _nonzero_components_2d(tensor, coord_names, symmetric=True)
    elif cls_name in (
        "RicciTensor",
        "EinsteinTensor",
        "SchoutenTensor",
        "StressEnergyTensor",
    ):
        return _nonzero_components_2d(tensor, coord_names, symmetric=False)
    elif cls_name == "ChristoffelSymbols":
        return _nonzero_components_3d(tensor, coord_names)
    elif cls_name in ("RiemannTensor", "WeylTensor"):
        return _nonzero_components_4d(tensor, coord_names)
    else:
        raise TypeError(f"Unknown tensor type: {cls_name}")


def components_table(items, tensor=None, *, tensor_symbol=None, index_style=None):
    r"""Format a list of (label, value) pairs as an IPython Markdown table.

    Parameters
    ----------
    items : list of (str, expression)
    tensor : tensor object, optional
        If provided, *tensor_symbol* and *index_style* are inferred from its type.
    tensor_symbol : str, optional
        The symbol for the tensor, e.g. ``"g"``, ``"\\Gamma"``, ``"R"``.
        Used as fallback when *tensor* is None.
    index_style : str, optional
        ``"sub"`` for subscript labels like :math:`g_{tt}`,
        ``"christoffel"`` for mixed like :math:`\Gamma^a_{\ bc}`,
        ``"riemann"`` for mixed like :math:`R^a{}_{bcd}`.
        Used as fallback when *tensor* is None.
    """
    from IPython.display import Markdown

    if tensor is not None:
        cls_name = type(tensor).__name__
        config = _TENSOR_TABLE_CONFIG.get(cls_name)
        if config is not None:
            tensor_symbol = tensor_symbol or config[0]
            index_style = index_style or config[1]

    # Defaults for fully manual usage
    if tensor_symbol is None:
        tensor_symbol = "T"
    if index_style is None:
        index_style = "sub"

    if not items:
        return Markdown(f"All components of ${tensor_symbol}$ are zero.")

    rows = []
    for label, val in items:
        if index_style == "sub":
            tex = f"${tensor_symbol}_{{{label}}}$"
        elif index_style == "super":
            tex = f"${tensor_symbol}^{{{label}}}$"
        else:
            # label already contains the index structure
            tex = f"${tensor_symbol}{label}$"
        if hasattr(val, "to_latex"):
            from gravica.simplify import simplify

            val_latex = _fix_der_latex(
                simplify(val).to_latex().removeprefix("$$").removesuffix("$$")
            )
        else:
            val_latex = str(val)
        rows.append(f"| {tex} | ${val_latex}$ |")

    table = "| Component | Value |\n|-----------|-------|\n" + "\n".join(rows)
    return Markdown(table)
