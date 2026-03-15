# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
uv sync                                        # install deps
uv run pytest                                  # all tests (tests/ -v by default)
uv run pytest tests/test_riemann.py            # single file
uv run pytest tests/test_riemann.py::test_foo  # single test
uv run benchmarks/run_benchmarks.py            # run benchmarks
uv run benchmarks/plot_benchmarks.py           # generate charts
uv run ruff check .                            # lint (ruff)
uv run ruff check . --fix                      # auto-fix lint warnings
uv run ruff format .                           # format (ruff)
uv build                                       # build dist/
uv run twine upload --repository testpypi dist/* # publish to TestPyPI
uv run twine upload dist/*                     # publish to PyPI
```

## Architecture

Pipeline: `MetricTensor → ChristoffelSymbols → RiemannTensor → RicciTensor → EinsteinTensor / WeylTensor`

Each tensor class takes its upstream tensor in `__init__`, lazy-computes via `@property`, and supports `__getitem__` indexing. Symmetries are exploited to halve computation.

All symbolic math uses `symbolica.Expression` (not SymPy). Zero-skipping via `str_is_zero()` is critical for performance. Simplification is `together().cancel()` with `expand().cancel()` fallback.

Test fixtures in `conftest.py` are `scope="session"`. Use `assert_zero()` for zero-checks in tests.
