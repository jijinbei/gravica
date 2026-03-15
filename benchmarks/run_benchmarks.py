"""Performance benchmarks: Gravica (Symbolica) vs EinsteinPy (SymPy).

Usage: uv run python benchmarks/run_benchmarks.py
"""

import gc
import json
import time
from pathlib import Path

import sympy
from sympy import symbols, sin, diag

# ── Metric definitions ─────────────────────────────────────────────


def _sympy_schwarzschild():
    t, r, theta, phi, r_s = symbols("t r theta phi r_s")
    f = 1 - r_s / r
    g = diag(f, -1 / f, -(r**2), -(r**2) * sin(theta) ** 2).tolist()
    return g, (t, r, theta, phi)


def _sympy_minkowski():
    t, x, y, z = symbols("t x y z")
    g = diag(1, -1, -1, -1).tolist()
    return g, (t, x, y, z)


def _sympy_flrw():
    t, r, theta, phi, k = symbols("t r theta phi k")
    a = sympy.Function("a")(t)
    g = diag(
        1, -(a**2) / (1 - k * r**2), -(a**2) * r**2, -(a**2) * r**2 * sin(theta) ** 2
    ).tolist()
    return g, (t, r, theta, phi)


# ── Benchmark runners ──────────────────────────────────────────────


def bench_atlas(metric_name: str) -> dict:
    """Run Gravica benchmarks for a given metric."""
    from gravica.christoffel import ChristoffelSymbols
    from gravica.riemann import RiemannTensor
    from gravica.ricci import RicciTensor, ricci_scalar
    from gravica.einstein import EinsteinTensor
    from gravica.weyl import WeylTensor

    if metric_name == "schwarzschild":
        from gravica.metrics.schwarzschild import schwarzschild

        m = schwarzschild()
    elif metric_name == "minkowski":
        from gravica.metrics.minkowski import minkowski

        m = minkowski()
    elif metric_name == "flrw":
        from gravica.metrics.flrw import flrw

        m = flrw()
    else:
        raise ValueError(f"Unknown metric: {metric_name}")

    results = {}

    # Metric inverse
    gc.disable()
    t0 = time.perf_counter()
    _ = m.inverse
    results["metric_inverse"] = time.perf_counter() - t0
    gc.enable()

    # Christoffel
    gc.disable()
    t0 = time.perf_counter()
    ch = ChristoffelSymbols(m)
    _ = ch.second_kind
    results["christoffel"] = time.perf_counter() - t0
    gc.enable()

    # Riemann
    gc.disable()
    t0 = time.perf_counter()
    riem = RiemannTensor(ch)
    _ = riem.components
    results["riemann"] = time.perf_counter() - t0
    gc.enable()

    # Ricci
    gc.disable()
    t0 = time.perf_counter()
    ric = RicciTensor(riem)
    _ = ric.components
    results["ricci"] = time.perf_counter() - t0
    gc.enable()

    # Ricci scalar
    gc.disable()
    t0 = time.perf_counter()
    ricci_scalar(ric)
    results["ricci_scalar"] = time.perf_counter() - t0
    gc.enable()

    # Einstein
    gc.disable()
    t0 = time.perf_counter()
    ein = EinsteinTensor(ric)
    _ = ein.components
    results["einstein"] = time.perf_counter() - t0
    gc.enable()

    # Weyl
    gc.disable()
    t0 = time.perf_counter()
    weyl = WeylTensor(riem, ric)
    _ = weyl.components
    results["weyl"] = time.perf_counter() - t0
    gc.enable()

    return results


def bench_einsteinpy(metric_name: str) -> dict:
    """Run EinsteinPy benchmarks for a given metric."""
    from einsteinpy.symbolic import (
        MetricTensor,
        ChristoffelSymbols,
        RiemannCurvatureTensor,
        RicciTensor,
        RicciScalar,
        EinsteinTensor,
    )

    if metric_name == "schwarzschild":
        g, coords = _sympy_schwarzschild()
    elif metric_name == "minkowski":
        g, coords = _sympy_minkowski()
    elif metric_name == "flrw":
        g, coords = _sympy_flrw()
    else:
        raise ValueError(f"Unknown metric: {metric_name}")

    results = {}

    # Metric inverse
    gc.disable()
    t0 = time.perf_counter()
    metric = MetricTensor(g, coords)
    _ = metric.tensor()  # Force evaluation
    results["metric_inverse"] = time.perf_counter() - t0
    gc.enable()

    # Christoffel
    gc.disable()
    t0 = time.perf_counter()
    ch = ChristoffelSymbols.from_metric(metric)
    _ = ch.tensor()
    results["christoffel"] = time.perf_counter() - t0
    gc.enable()

    # Riemann
    gc.disable()
    t0 = time.perf_counter()
    riem = RiemannCurvatureTensor.from_christoffels(ch)
    _ = riem.tensor()
    results["riemann"] = time.perf_counter() - t0
    gc.enable()

    # Ricci
    gc.disable()
    t0 = time.perf_counter()
    ric = RicciTensor.from_riemann(riem)
    _ = ric.tensor()
    results["ricci"] = time.perf_counter() - t0
    gc.enable()

    # Ricci scalar
    gc.disable()
    t0 = time.perf_counter()
    R = RicciScalar.from_riccitensor(ric)
    _ = R.tensor()
    results["ricci_scalar"] = time.perf_counter() - t0
    gc.enable()

    # Einstein (uses from_metric in current EinsteinPy)
    gc.disable()
    t0 = time.perf_counter()
    ein = EinsteinTensor.from_metric(metric)
    _ = ein.tensor()
    results["einstein"] = time.perf_counter() - t0
    gc.enable()

    # Weyl: EinsteinPy doesn't have Weyl, skip
    results["weyl"] = None

    return results


def run_single(metric_name: str) -> dict:
    """Benchmark a single metric with warmup + 3 runs (median)."""
    import statistics

    print(f"\n{'=' * 60}")
    print(f"  Benchmarking: {metric_name}")
    print(f"{'=' * 60}")

    # Gravica
    atlas_runs = []
    for i in range(3):
        if i == 0:
            print("  Gravica warmup...", end="", flush=True)
        r = bench_atlas(metric_name)
        if i == 0:
            print(" done")
        atlas_runs.append(r)

    # EinsteinPy
    epy_runs = []
    for i in range(3):
        if i == 0:
            print("  EinsteinPy warmup...", end="", flush=True)
        r = bench_einsteinpy(metric_name)
        if i == 0:
            print(" done")
        epy_runs.append(r)

    steps = [
        "metric_inverse",
        "christoffel",
        "riemann",
        "ricci",
        "ricci_scalar",
        "einstein",
        "weyl",
    ]
    result = {}
    for step in steps:
        atlas_times = [run[step] for run in atlas_runs]
        atlas_median = statistics.median(atlas_times)

        epy_times = [run[step] for run in epy_runs if run[step] is not None]
        epy_median = statistics.median(epy_times) if epy_times else None

        result[step] = {
            "gravica": atlas_median,
            "einsteinpy": epy_median,
            "speedup": epy_median / atlas_median
            if epy_median and atlas_median > 0
            else None,
        }

    return result


def format_time(seconds: float | None) -> str:
    if seconds is None:
        return "N/A"
    if seconds < 0.001:
        return f"{seconds * 1e6:.0f}μs"
    if seconds < 1:
        return f"{seconds * 1e3:.1f}ms"
    return f"{seconds:.2f}s"


def main():
    metrics = ["minkowski", "schwarzschild", "flrw"]
    all_results = {}

    for metric_name in metrics:
        try:
            all_results[metric_name] = run_single(metric_name)
        except Exception as e:
            print(f"  ERROR: {e}")
            all_results[metric_name] = {"error": str(e)}

    # Print Markdown table
    steps = [
        "metric_inverse",
        "christoffel",
        "riemann",
        "ricci",
        "ricci_scalar",
        "einstein",
        "weyl",
    ]

    print(f"\n\n{'=' * 80}")
    print("RESULTS")
    print(f"{'=' * 80}\n")

    print("| Computation | Metric | Gravica | EinsteinPy | Speedup |")
    print("|---|---|---|---|---|")

    for step in steps:
        for metric_name in metrics:
            r = all_results.get(metric_name, {}).get(step, {})
            if isinstance(r, dict) and "gravica" in r:
                atlas_t = format_time(r["gravica"])
                epy_t = format_time(r["einsteinpy"])
                speedup = f"{r['speedup']:.1f}x" if r.get("speedup") else "N/A"
                print(f"| {step} | {metric_name} | {atlas_t} | {epy_t} | {speedup} |")

    # Save JSON
    out_path = Path("benchmarks/results.json")
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nResults saved to {out_path}")


if __name__ == "__main__":
    main()
