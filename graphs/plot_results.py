"""
plot_results.py
---------------
Generates all graphs for the DNA Sequence Similarity project.

Usage (run from the DNASequenceSimilarity/ directory):
    python graphs/plot_results.py

Requires:
    pip install matplotlib pandas
"""

import os
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")          # non-interactive backend (saves to file)
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ── output directory ────────────────────────────────────────────────────────
GRAPH_DIR = os.path.dirname(os.path.abspath(__file__))

def savefig(name):
    path = os.path.join(GRAPH_DIR, name)
    plt.savefig(path, dpi=150, bbox_inches="tight")
    print(f"  Saved: {path}")
    plt.close()

# ── colour scheme ────────────────────────────────────────────────────────────
PALETTE = ["#4C72B0", "#DD8452", "#55A868", "#C44E52",
           "#8172B3", "#937860", "#DA8BC3", "#8C8C8C"]
plt.rcParams.update({
    "figure.facecolor": "#0D1117",
    "axes.facecolor":   "#161B22",
    "axes.edgecolor":   "#30363D",
    "axes.labelcolor":  "#C9D1D9",
    "xtick.color":      "#C9D1D9",
    "ytick.color":      "#C9D1D9",
    "text.color":       "#C9D1D9",
    "grid.color":       "#30363D",
    "grid.linestyle":   "--",
    "grid.alpha":       0.5,
    "legend.facecolor": "#161B22",
    "legend.edgecolor": "#30363D",
    "font.family":      "sans-serif",
    "axes.titlesize":   13,
    "axes.labelsize":   11,
})

# ────────────────────────────────────────────────────────────────────────────
# 1. Calibration: time_small vs time_large per k
# ────────────────────────────────────────────────────────────────────────────
def plot_calibration():
    path = "calibration_log.csv"
    if not os.path.exists(path):
        print(f"  [skip] {path} not found")
        return

    df = pd.read_csv(path)
    # Filter for 50000 bp primary calibration
    if "calib_len" in df.columns:
        sub50 = df[df["calib_len"] == 50000]
    else:
        sub50 = df

    small = sub50[(sub50["strategy"] == "small") & (sub50["avg_time_ms"] > 0)].sort_values("k")
    large = sub50[sub50["strategy"] == "large"].sort_values("k")

    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.plot(small["k"], small["avg_time_ms"], "o-", color=PALETTE[0],
            label="Strategy A: Direct 4^k Array", linewidth=2.0, markersize=5)
    ax.plot(large["k"], large["avg_time_ms"], "s-", color=PALETTE[1],
            label="Strategy B: Integer Sort + Merge-Walk", linewidth=2.0, markersize=5)

    # Detect crossover
    merged = pd.merge(small[["k","avg_time_ms"]], large[["k","avg_time_ms"]],
                      on="k", suffixes=("_s","_l"))
    cross = merged[merged["avg_time_ms_l"] < merged["avg_time_ms_s"] * 0.95]
    if not cross.empty:
        k_thresh = int(cross.iloc[0]["k"])
        ax.axvline(k_thresh, color="#E34A33", linestyle="--", linewidth=1.8,
                   label=f"Experimental threshold (k={k_thresh})")
        ax.annotate(f"Crossover @ k={k_thresh}",
                    xy=(k_thresh, cross.iloc[0]["avg_time_ms_l"]),
                    xytext=(k_thresh + 2, cross.iloc[0]["avg_time_ms_l"] + 15),
                    arrowprops=dict(facecolor="#E34A33", shrink=0.08, width=1.5, headwidth=7),
                    fontweight="bold", color="#E34A33")

    ax.axvspan(13.5, 32.5, color="#f0f0f0", alpha=0.5, label="Memory Wall (>268 MB array)")
    ax.set_xlabel("k (k-mer length)")
    ax.set_ylabel("Average Time (ms)")
    ax.set_title("Calibration Crossover: Strategy A (Direct Array) vs Strategy B (Integer Sort)")
    ax.legend(loc="upper left")
    ax.grid(True)
    ax.set_xticks(range(1, 33))
    savefig("calibration_crossover.png")


# ────────────────────────────────────────────────────────────────────────────
# 2. Experiment A: Execution time vs Sequence length for selected k
# ────────────────────────────────────────────────────────────────────────────
def plot_exp_a():
    path = "results_benchmark.csv"
    if not os.path.exists(path):
        print(f"  [skip] {path} not found")
        return

    df = pd.read_csv(path)
    avg = df.groupby(["seq_len", "k"])["total_ms"].mean().reset_index()

    selected_k = [3, 8, 16, 24, 32]
    fig, ax = plt.subplots(figsize=(10, 6))

    for i, k in enumerate(selected_k):
        sub = avg[avg["k"] == k].sort_values("seq_len")
        ax.plot(sub["seq_len"], sub["total_ms"], "o-",
                color=PALETTE[i], label=f"k={k}", linewidth=1.8, markersize=5)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Sequence Length (log scale)")
    ax.set_ylabel("Total Execution Time  ms  (log scale)")
    ax.set_title("Experiment A — Execution Time vs Sequence Length")
    ax.legend(title="k value")
    ax.grid(True, which="both")
    savefig("exp_a_time_vs_length.png")

# ────────────────────────────────────────────────────────────────────────────
# 3. Experiment B: Execution time vs k for each sequence length
# ────────────────────────────────────────────────────────────────────────────
def plot_exp_b():
    path = "results_benchmark.csv"
    if not os.path.exists(path):
        print(f"  [skip] {path} not found")
        return

    df = pd.read_csv(path)
    avg = df.groupby(["seq_len", "k"])["total_ms"].mean().reset_index()
    lengths = sorted(avg["seq_len"].unique())

    fig, ax = plt.subplots(figsize=(12, 6))
    for i, L in enumerate(lengths):
        sub = avg[avg["seq_len"] == L].sort_values("k")
        ax.plot(sub["k"], sub["total_ms"], "o-",
                color=PALETTE[i % len(PALETTE)],
                label=f"len={L:,}", linewidth=1.8, markersize=4)

    ax.set_xlabel("k  (k-mer length)")
    ax.set_ylabel("Average Execution Time  (ms)")
    ax.set_title("Experiment B — Execution Time vs k for Different Sequence Lengths")
    ax.legend(title="Sequence length")
    ax.set_xticks(range(1, 33))
    ax.grid(True)
    savefig("exp_b_time_vs_k.png")

# ────────────────────────────────────────────────────────────────────────────
# 4. Phase breakdown: counting time vs Jaccard time
# ────────────────────────────────────────────────────────────────────────────
def plot_phase_breakdown():
    path = "results_benchmark.csv"
    if not os.path.exists(path):
        return

    df = pd.read_csv(path)
    max_len = df["seq_len"].max()
    sub = df[df["seq_len"] == max_len].groupby("k")[
        ["time_count_ms", "time_jaccard_ms"]].mean().reset_index()

    x = sub["k"]
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.stackplot(x, sub["time_count_ms"], sub["time_jaccard_ms"],
                 labels=["k-mer counting", "Jaccard computation"],
                 colors=[PALETTE[0], PALETTE[2]], alpha=0.85)
    ax.set_xlabel("k  (k-mer length)")
    ax.set_ylabel("Time (ms)")
    ax.set_title(f"Phase Breakdown — Sequence Length = {max_len:,}")
    ax.legend(loc="upper left")
    ax.set_xticks(range(1, 33))
    ax.grid(True)
    savefig("phase_breakdown.png")

# ────────────────────────────────────────────────────────────────────────────
# 5. Experiment C: Similarity vs mutation percentage for each k
# ────────────────────────────────────────────────────────────────────────────
def plot_exp_c():
    path = "results_mutation.csv"
    if not os.path.exists(path):
        print(f"  [skip] {path} not found")
        return

    df = pd.read_csv(path)
    avg = df.groupby(["k", "mutation_pct"])["similarity"].mean().reset_index()
    k_vals = sorted(avg["k"].unique())

    fig, ax = plt.subplots(figsize=(11, 6))
    highlight_k = {1, 2, 4, 8, 12, 16, 20, 24, 28, 32}
    cmap = matplotlib.colormaps.get_cmap("viridis")

    for k in k_vals:
        sub = avg[avg["k"] == k].sort_values("mutation_pct")
        norm_k = (k - 1) / 31.0
        c = cmap(norm_k)
        if k in highlight_k or len(k_vals) <= 10:
            ax.plot(sub["mutation_pct"], sub["similarity"], "o-",
                    color=c, label=f"k={k}", linewidth=2.0, markersize=5)
        else:
            ax.plot(sub["mutation_pct"], sub["similarity"], "-",
                    color=c, alpha=0.35, linewidth=1.0)

    ax.set_xlabel("Mutation Percentage  (%)")
    ax.set_ylabel("Generalized Jaccard Similarity")
    ax.set_title("Experiment C — Sensitivity Decay across all k in [1, 32]")
    ax.set_ylim(-0.05, 1.05)
    ax.legend(title="Selected k values", ncol=2, loc="upper right")
    ax.grid(True)
    savefig("exp_c_mutation_vs_similarity.png")


# ────────────────────────────────────────────────────────────────────────────
# 6. Heatmap: avg total_ms  [seq_len x k]
# ────────────────────────────────────────────────────────────────────────────
def plot_heatmap():
    path = "results_benchmark.csv"
    if not os.path.exists(path):
        return

    df = pd.read_csv(path)
    pivot = df.groupby(["seq_len","k"])["total_ms"].mean().unstack(level="k")

    fig, ax = plt.subplots(figsize=(16, 5))
    im = ax.imshow(pivot.values, aspect="auto", cmap="inferno",
                   norm=matplotlib.colors.LogNorm(
                       vmin=max(pivot.values.min(), 1e-3),
                       vmax=pivot.values.max()))

    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, fontsize=8)
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels([f"{v:,}" for v in pivot.index])
    ax.set_xlabel("k  (k-mer length)")
    ax.set_ylabel("Sequence Length")
    ax.set_title("Heatmap — Average Total Time (ms) per [seq_len x k]  (log colour scale)")

    cb = plt.colorbar(im, ax=ax, pad=0.01)
    cb.set_label("Time (ms)")
    savefig("heatmap_time.png")

# ────────────────────────────────────────────────────────────────────────────
# 7. Strategy comparison at the threshold boundary
# ────────────────────────────────────────────────────────────────────────────
def plot_strategy_comparison():
    path = "results_benchmark.csv"
    if not os.path.exists(path):
        return

    df = pd.read_csv(path)
    avg = df.groupby(["seq_len","k","strategy"])["total_ms"].mean().reset_index()

    lengths = sorted(avg["seq_len"].unique())
    pick_len = lengths[len(lengths)//2]
    sub = avg[avg["seq_len"] == pick_len]

    small_df = sub[sub["strategy"] == "small"].sort_values("k")
    large_df = sub[sub["strategy"] == "large"].sort_values("k")

    fig, ax = plt.subplots(figsize=(10, 5))
    if not small_df.empty:
        ax.bar(small_df["k"] - 0.2, small_df["total_ms"], width=0.4,
               color=PALETTE[0], label="Small-k strategy", alpha=0.85)
    if not large_df.empty:
        ax.bar(large_df["k"] + 0.2, large_df["total_ms"], width=0.4,
               color=PALETTE[1], label="Large-k strategy", alpha=0.85)

    ax.set_xlabel("k  (k-mer length)")
    ax.set_ylabel("Avg Total Time  (ms)")
    ax.set_title(f"Strategy Comparison — Sequence Length = {pick_len:,}")
    ax.legend()
    ax.set_xticks(range(1, 33))
    ax.grid(True, axis="y")
    savefig("strategy_comparison.png")


# ────────────────────────────────────────────────────────────────────────────
# MAIN
# ────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    # Move to DNASequenceSimilarity/ directory (parent of graphs/)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    os.chdir(project_dir)

    print("Generating graphs...")
    plot_calibration()
    plot_exp_a()
    plot_exp_b()
    plot_phase_breakdown()
    plot_exp_c()
    plot_heatmap()
    plot_strategy_comparison()
    print("\nDone. All graphs saved to graphs/")
