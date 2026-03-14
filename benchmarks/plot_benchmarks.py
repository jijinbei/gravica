"""Generate benchmark charts from results.json.

Usage: uv run python benchmarks/plot_benchmarks.py
"""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

RESULTS_PATH = Path("benchmarks/results.json")
OUTPUT_DIR = Path("docs")

STEP_LABELS = {
    "metric_inverse": "Metric\nInverse",
    "christoffel": "Christoffel",
    "riemann": "Riemann",
    "ricci": "Ricci",
    "ricci_scalar": "Ricci\nScalar",
    "einstein": "Einstein",
}

METRIC_LABELS = {
    "minkowski": "Minkowski",
    "schwarzschild": "Schwarzschild",
    "flrw": "FLRW",
}

COLORS = {
    "gravica": "#2563EB",
    "einsteinpy": "#DC2626",
}


def load_results() -> dict:
    with open(RESULTS_PATH) as f:
        return json.load(f)


def to_ms(seconds: float | None) -> float | None:
    if seconds is None:
        return None
    return seconds * 1000


def plot_time_comparison(results: dict) -> None:
    """Bar chart: absolute time comparison per computation step, grouped by metric."""
    steps = list(STEP_LABELS.keys())
    metrics = list(METRIC_LABELS.keys())

    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)
    fig.suptitle("Computation Time: Gravica (Symbolica) vs EinsteinPy (SymPy)", fontsize=14, fontweight="bold", y=1.02)

    for ax, metric in zip(axes, metrics):
        atlas_times = []
        epy_times = []
        for step in steps:
            data = results[metric][step]
            atlas_times.append(to_ms(data["gravica"]))
            epy_times.append(to_ms(data["einsteinpy"]))

        x = np.arange(len(steps))
        width = 0.35

        bars_atlas = ax.bar(x - width / 2, atlas_times, width, label="Gravica", color=COLORS["gravica"], zorder=3)
        bars_epy = ax.bar(x + width / 2, epy_times, width, label="EinsteinPy", color=COLORS["einsteinpy"], zorder=3)

        ax.set_title(METRIC_LABELS[metric], fontsize=12, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels([STEP_LABELS[s] for s in steps], fontsize=8)
        ax.set_yscale("log")
        ax.set_ylabel("Time (ms)" if ax == axes[0] else "")
        ax.grid(axis="y", alpha=0.3, zorder=0)
        ax.set_axisbelow(True)

    axes[0].legend(fontsize=10, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "bench_time_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUTPUT_DIR / 'bench_time_comparison.png'}")


def plot_speedup(results: dict) -> None:
    """Horizontal bar chart: speedup factor per computation step for each metric."""
    steps = list(STEP_LABELS.keys())
    metrics = list(METRIC_LABELS.keys())

    fig, ax = plt.subplots(figsize=(12, 6))

    y_positions = []
    y_labels = []
    colors = []
    speedups = []
    metric_colors = {"minkowski": "#3B82F6", "schwarzschild": "#F59E0B", "flrw": "#10B981"}

    group_gap = 0.6
    bar_height = 0.25
    y = 0

    for step in reversed(steps):
        for metric in metrics:
            data = results[metric][step]
            sp = data.get("speedup")
            if sp is not None and sp > 0:
                y_positions.append(y)
                y_labels.append(f"{METRIC_LABELS[metric]}")
                colors.append(metric_colors[metric])
                speedups.append(sp)
                y += bar_height + 0.05
        y += group_gap

    bars = ax.barh(y_positions, speedups, height=bar_height, color=colors, zorder=3, edgecolor="white", linewidth=0.5)

    # Add value labels
    for bar, sp in zip(bars, speedups):
        if sp >= 1:
            label = f"{sp:.0f}x" if sp >= 10 else f"{sp:.1f}x"
            ax.text(bar.get_width() * 1.02, bar.get_y() + bar.get_height() / 2,
                    label, va="center", fontsize=8, fontweight="bold")

    # Add step group labels on the left
    y = 0
    step_label_positions = []
    for step in reversed(steps):
        count = sum(1 for m in metrics if results[m][step].get("speedup") is not None and results[m][step]["speedup"] > 0)
        if count > 0:
            mid = y + (count - 1) * (bar_height + 0.05) / 2
            step_label_positions.append((mid, STEP_LABELS[step].replace("\n", " ")))
            y += count * (bar_height + 0.05) + group_gap

    ax.set_yticks([p for p, _ in step_label_positions])
    ax.set_yticklabels([l for _, l in step_label_positions], fontsize=10, fontweight="bold")

    ax.set_xscale("log")
    ax.set_xlabel("Speedup (Gravica / EinsteinPy)", fontsize=11)
    ax.set_title("Gravica Speedup over EinsteinPy", fontsize=14, fontweight="bold")
    ax.axvline(x=1, color="gray", linestyle="--", linewidth=0.8, zorder=2)
    ax.grid(axis="x", alpha=0.3, zorder=0)
    ax.set_axisbelow(True)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=metric_colors[m], label=METRIC_LABELS[m]) for m in metrics]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=10)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "bench_speedup.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUTPUT_DIR / 'bench_speedup.png'}")


def plot_summary_table(results: dict) -> None:
    """Heatmap-style summary of speedup values."""
    steps = [s for s in STEP_LABELS.keys()]
    metrics = list(METRIC_LABELS.keys())

    data = []
    for step in steps:
        row = []
        for metric in metrics:
            sp = results[metric][step].get("speedup")
            row.append(sp if sp and sp > 0 else float("nan"))
        data.append(row)

    data_arr = np.array(data)

    fig, ax = plt.subplots(figsize=(8, 5))

    # Use log scale for color mapping
    from matplotlib.colors import LogNorm
    valid = data_arr[~np.isnan(data_arr)]
    vmin = max(valid.min(), 0.1)
    vmax = valid.max()

    cmap = plt.cm.RdYlGn
    im = ax.imshow(data_arr, cmap=cmap, norm=LogNorm(vmin=vmin, vmax=vmax), aspect="auto")

    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels([METRIC_LABELS[m] for m in metrics], fontsize=11)
    ax.set_yticks(range(len(steps)))
    ax.set_yticklabels([STEP_LABELS[s].replace("\n", " ") for s in steps], fontsize=10)

    # Annotate cells
    for i in range(len(steps)):
        for j in range(len(metrics)):
            val = data_arr[i, j]
            if np.isnan(val):
                text = "N/A"
                color = "gray"
            elif val < 1:
                text = f"{val:.1f}x"
                color = "white"
            elif val < 10:
                text = f"{val:.1f}x"
                color = "black"
            else:
                text = f"{val:.0f}x"
                color = "black" if val < 500 else "white"
            ax.text(j, i, text, ha="center", va="center", fontsize=11, fontweight="bold", color=color)

    ax.set_title("Speedup: Gravica over EinsteinPy\n(higher = Gravica is faster)", fontsize=13, fontweight="bold")
    cbar = fig.colorbar(im, ax=ax, label="Speedup factor", shrink=0.8)
    cbar.set_label("Speedup (log scale)", fontsize=10)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "bench_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUTPUT_DIR / 'bench_heatmap.png'}")


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    results = load_results()
    plot_time_comparison(results)
    plot_speedup(results)
    plot_summary_table(results)
    print("\nAll charts generated.")


if __name__ == "__main__":
    main()
