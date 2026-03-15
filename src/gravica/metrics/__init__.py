r"""Pre-defined spacetime metrics.

This subpackage provides factory functions that return
:class:`~gravica.metric.MetricTensor` instances for well-known spacetimes:

* :func:`minkowski` -- flat Minkowski metric
* :func:`schwarzschild` -- Schwarzschild black-hole metric
* :func:`kerr` -- Kerr rotating black-hole metric
* :func:`flrw` -- Friedmann--Lemaître--Robertson--Walker cosmological metric
* :func:`godel` -- Gödel rotating-universe metric
"""

from gravica.metrics.minkowski import minkowski
from gravica.metrics.schwarzschild import schwarzschild
from gravica.metrics.kerr import kerr
from gravica.metrics.flrw import flrw
from gravica.metrics.godel import godel

__all__ = ["minkowski", "schwarzschild", "kerr", "flrw", "godel"]
