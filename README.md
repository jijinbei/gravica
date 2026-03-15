# Gravica

[![PyPI](https://img.shields.io/pypi/v/gravica)](https://pypi.org/project/gravica/)

General Relativity computation library built on [Symbolica](https://symbolica.io/) (Rust-powered CAS).

Gravica computes the full GR tensor chain — from metric tensor to Einstein tensor — using Symbolica's high-performance symbolic algebra engine, achieving **23x–4300x speedup** over EinsteinPy/SymPy.

## Features

- Full GR computation chain: Metric → Christoffel → Riemann → Ricci → Einstein → Weyl
- Curvature invariants: Kretschner scalar, Ricci scalar
- Schouten tensor, Stress-Energy-Momentum tensor
- Geodesic equation generator
- Index raising/lowering utilities
- Built-in metrics: Minkowski, Schwarzschild, Kerr, FLRW, Reissner-Nordström, de Sitter, anti-de Sitter, Gödel
- Lazy evaluation with caching
- Cross-validated against EinsteinPy

## Documentation

API reference: [https://site.jijinbei.jp/gravica/](https://site.jijinbei.jp/gravica/)

## Quick Start

```bash
uv add gravica
```

```python
from gravica.metrics.schwarzschild import schwarzschild
from gravica.christoffel import ChristoffelSymbols
from gravica.riemann import RiemannTensor
from gravica.ricci import RicciTensor, ricci_scalar
from gravica.einstein import EinsteinTensor

metric = schwarzschild()
christoffel = ChristoffelSymbols(metric)
riemann = RiemannTensor(christoffel)
ricci = RicciTensor(riemann)
einstein = EinsteinTensor(ricci)

# Verify vacuum solution: G_ab = 0
for a in range(4):
    for b in range(4):
        assert str(einstein[a, b]) == "0"
```

## Benchmarks: Gravica vs EinsteinPy

All benchmarks measured on the same machine. Median of 3 runs with GC disabled.

### Speedup Heatmap

![Speedup Heatmap](docs/bench_heatmap.png)

### Absolute Time Comparison

![Time Comparison](docs/bench_time_comparison.png)

### Summary

| Computation | Minkowski | Schwarzschild | FLRW |
|---|---|---|---|
| Christoffel | **33x** | **23x** | **23x** |
| Riemann | **91x** | **149x** | **147x** |
| Ricci | **72x** | **40x** | **296x** |
| Ricci Scalar | **472x** | **2021x** | **544x** |
| Einstein | **1391x** | **4342x** | **1029x** |

> Metric inverse is ~0.3–0.5x (Python cofactor overhead), but this is amortized by the massive speedups in downstream computations.

### Reproduce

```bash
uv run benchmarks/run_benchmarks.py    # Run benchmarks
uv run benchmarks/plot_benchmarks.py   # Generate charts
```

## Architecture

```
MetricTensor → ChristoffelSymbols → RiemannTensor → RicciTensor → EinsteinTensor
                      ↓                   ↓              ↓       → WeylTensor
               GeodesicEquations   KretschnerScalar  SchoutenTensor
                                                              ↓
                                                     StressEnergyTensor
```

| Module | Computes |
|---|---|
| `metric.py` | $g_{ab}$, $g^{ab}$, $\det(g)$ |
| `christoffel.py` | $\Gamma^a_{\ bc} = g^{ad}\,\tfrac{1}{2}(\partial_b\,g_{ac} + \partial_c\,g_{ab} - \partial_a\,g_{bc})$ |
| `riemann.py` | $R^a_{\ bcd}$, $R_{abcd}$, $R^{abcd}$ |
| `ricci.py` | $R_{ab} = R^c_{\ acb}$, $R = g^{ab}\,R_{ab}$ |
| `einstein.py` | $G_{ab} = R_{ab} - \tfrac{1}{2}\,g_{ab}\,R$ |
| `weyl.py` | $C_{abcd}$ (Weyl conformal tensor) |
| `kretschner.py` | $K = R_{abcd}\,R^{abcd}$ (Kretschner scalar) |
| `geodesic.py` | $\ddot{x}^a + \Gamma^a_{\ bc}\,\dot{x}^b\,\dot{x}^c = 0$ |
| `schouten.py` | $S_{ab} = \tfrac{1}{n-2}\bigl(R_{ab} - \tfrac{R\,g_{ab}}{2(n-1)}\bigr)$ |
| `stress_energy.py` | $8\pi G\,T_{ab} = G_{ab} + \Lambda\,g_{ab}$ |
| `indexing.py` | Index raising / lowering for rank-2 tensors |

## Tests

```bash
uv run pytest
```

Verified properties:
- **Minkowski**: All tensors $= 0$
- **Schwarzschild**: $R_{ab} = 0$, $G_{ab} = 0$ (vacuum), $K = 12\,r_s^2/r^6$
- **Riemann symmetries**: $R^a_{\ bcd} = -R^a_{\ bdc}$
- **Christoffel known values**: $\Gamma^r_{\ tt} = r_s(r-r_s)/(2r^3)$
- **de Sitter / anti-de Sitter**: Ricci scalar matches analytic values
- **Geodesic equations**: Free particle in Minkowski
- **Index roundtrip**: Raise then lower recovers original tensor
- **EinsteinPy cross-validation**: Christoffel and Ricci match

## License

MIT
