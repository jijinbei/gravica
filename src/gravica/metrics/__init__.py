r"""Pre-defined spacetime metrics.

This subpackage provides factory functions that return
:class:`~gravica.metric.MetricTensor` instances for well-known spacetimes:

* :func:`minkowski` -- flat Minkowski metric
* :func:`schwarzschild` -- Schwarzschild black-hole metric
* :func:`kerr` -- Kerr rotating black-hole metric
* :func:`flrw` -- Friedmann--Lemaître--Robertson--Walker cosmological metric
* :func:`godel` -- Gödel rotating-universe metric
* :func:`reissner_nordstrom` -- Reissner-Nordström charged black-hole metric
* :func:`de_sitter` -- de Sitter metric (positive cosmological constant)
* :func:`anti_de_sitter` -- anti-de Sitter metric (negative cosmological constant)
"""

from gravica.metrics.minkowski import minkowski
from gravica.metrics.schwarzschild import schwarzschild
from gravica.metrics.kerr import kerr
from gravica.metrics.flrw import flrw
from gravica.metrics.godel import godel
from gravica.metrics.reissner_nordstrom import reissner_nordstrom
from gravica.metrics.de_sitter import de_sitter
from gravica.metrics.anti_de_sitter import anti_de_sitter

__all__ = [
    "minkowski",
    "schwarzschild",
    "kerr",
    "flrw",
    "godel",
    "reissner_nordstrom",
    "de_sitter",
    "anti_de_sitter",
]
