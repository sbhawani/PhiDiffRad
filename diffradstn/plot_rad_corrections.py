#!/usr/bin/env python3
"""
plot_rad_corrections.py
-----------------------
Reads diffrad_rc_results.csv produced by RunMDiffradNew and makes
publication-quality plots of the radiative correction δ_RC.

Usage:
    python3 plot_rad_corrections.py                          # default CSV name
    python3 plot_rad_corrections.py diffrad_rc_results.csv  # explicit path
    python3 plot_rad_corrections.py --period rgasp19_inb    # tag in filenames

Outputs (all saved to the same directory as the CSV):
    rad_corr_vs_tprime.pdf   δ_RC vs t'  (one curve per Q² bin)
    rad_corr_vs_tabs.pdf     δ_RC vs |t| (one curve per Q² bin)
    rad_corr_vs_Q2.pdf       δ_RC vs Q²  (one curve per t' bin)
    rad_corr_heatmap.pdf     2D heatmap δ_RC(Q², t')
    rad_corr_summary.pdf     all four panels on one page
"""

import sys
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.cm import get_cmap
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

# ── Argument parsing ─────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Plot MDIFFRAD radiative corrections")
parser.add_argument("csv", nargs="?", default="diffrad_rc_results.csv",
                    help="Path to diffrad_rc_results.csv")
parser.add_argument("--period", default="",
                    help="Run-period tag appended to output filenames (e.g. rgasp19_inb)")
parser.add_argument("--outdir", default="",
                    help="Output directory (default: same as CSV)")
args = parser.parse_args()

csv_path = args.csv
if not os.path.isfile(csv_path):
    sys.exit(f"[ERROR] CSV not found: {csv_path}")

out_dir = args.outdir if args.outdir else os.path.dirname(os.path.abspath(csv_path))
os.makedirs(out_dir, exist_ok=True)

tag = f"_{args.period}" if args.period else ""

def out_path(name):
    return os.path.join(out_dir, name.replace(".pdf", f"{tag}.pdf"))

# ── Load CSV ─────────────────────────────────────────────────────────────────
# Skip comment lines starting with #
df = pd.read_csv(csv_path, comment="#")
df.columns = df.columns.str.strip()

print(f"[plot_rad_corrections] Loaded {len(df)} rows from {csv_path}")
print(f"  Columns: {list(df.columns)}")

# Drop rows where rad_corr is NaN (failed MC jobs)
n_before = len(df)
df = df.dropna(subset=["rad_corr"])
if len(df) < n_before:
    print(f"  WARNING: dropped {n_before - len(df)} rows with NaN rad_corr")

# ── Metadata from CSV header comments ────────────────────────────────────────
beam_mom = None
W2_range = None
with open(csv_path) as f:
    for line in f:
        if not line.startswith("#"):
            break
        if "beam_momentum_GeV" in line:
            beam_mom = float(line.split("=")[1].strip())
        if "W2_range_GeV2" in line:
            parts = line.split("=")[1].strip().split()
            W2_range = (float(parts[0]), float(parts[1]))

beam_label = f"$E_{{beam}}$ = {beam_mom:.2f} GeV" if beam_mom else ""
w2_label   = (f"$W^2 \\in [{W2_range[0]:.1f},\\ {W2_range[1]:.1f}]$ GeV$^2$"
              if W2_range else "")
kin_label  = ",  ".join(filter(None, [beam_label, w2_label]))

# ── Bin index sets ────────────────────────────────────────────────────────────
q2_bins  = sorted(df["iQ2"].unique())
tp_bins  = sorted(df["itp"].unique())
nQ2      = len(q2_bins)
nTp      = len(tp_bins)

# ── Colour palettes ───────────────────────────────────────────────────────────
# Q² bins: blues/reds from a diverging map
q2_colors  = [get_cmap("RdYlBu")(i / max(nQ2 - 1, 1)) for i in range(nQ2)]
# t' bins: warm→cool
tp_colors  = [get_cmap("plasma")(i / max(nTp - 1, 1)) for i in range(nTp)]

markers_q2 = ["o", "s", "^", "D", "v", "P", "*", "X"]
markers_tp = ["o", "s", "^", "D", "v", "P", "*", "X", "h"]

# ── Helper: Q² and t' bin labels ─────────────────────────────────────────────
def q2_label(iQ2):
    row = df[df["iQ2"] == iQ2].iloc[0]
    return f"${row['Q2_lo']:.2f} \\leq Q^2 < {row['Q2_hi']:.2f}$ GeV$^2$"

def tp_label(itp):
    row = df[df["itp"] == itp].iloc[0]
    return f"${row['tprime_lo']:.3f} \\leq t' < {row['tprime_hi']:.3f}$ GeV$^2$"

# ── Shared axis limits ────────────────────────────────────────────────────────
rc_all  = df["rad_corr"].values
sig_all = df["rad_corr_sigma"].fillna(0).values
ylo = max(0.5, (rc_all - sig_all).min() - 0.05)
yhi = min(1.5, (rc_all + sig_all).max() + 0.05)

# ── Style ─────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":       "serif",
    "font.size":         11,
    "axes.labelsize":    12,
    "axes.titlesize":    12,
    "legend.fontsize":   8,
    "xtick.direction":   "in",
    "ytick.direction":   "in",
    "xtick.top":         True,
    "ytick.right":       True,
    "axes.linewidth":    1.2,
    "lines.linewidth":   1.5,
    "errorbar.capsize":  3,
})

# ── Plot 1: δ_RC vs t' ────────────────────────────────────────────────────────
def plot_vs_tprime(ax):
    for i, iQ2 in enumerate(q2_bins):
        sub = df[df["iQ2"] == iQ2].sort_values("tprime_c")
        ax.errorbar(sub["tprime_c"], sub["rad_corr"],
                    xerr=0.5*(sub["tprime_hi"] - sub["tprime_lo"]),
                    yerr=sub["rad_corr_sigma"].fillna(0),
                    fmt=markers_q2[i % len(markers_q2)],
                    color=q2_colors[i], label=q2_label(iQ2),
                    markersize=5, linewidth=1.2, capsize=3)
        # shaded ±σ band
        ax.fill_between(sub["tprime_c"],
                        sub["rad_corr"] - sub["rad_corr_sigma"].fillna(0),
                        sub["rad_corr"] + sub["rad_corr_sigma"].fillna(0),
                        color=q2_colors[i], alpha=0.12)
    # bin-edge guide lines
    for e in sorted(df["tprime_lo"].unique()):
        ax.axvline(e, color="gray", lw=0.6, ls=":")
    ax.axhline(1.0, color="black", lw=0.8, ls="--", alpha=0.5)
    ax.set_xlabel("$t' = |t| - |t_{\\min}|$ (GeV$^2$)")
    ax.set_ylabel("Radiative Correction $\\delta_{RC}$")
    ax.set_title("$\\delta_{RC}$ vs $t'$")
    ax.set_ylim(ylo, yhi)
    ax.legend(loc="best", framealpha=0.7, ncol=1)

# ── Plot 2: δ_RC vs |t| ───────────────────────────────────────────────────────
def plot_vs_tabs(ax):
    for i, iQ2 in enumerate(q2_bins):
        sub = df[df["iQ2"] == iQ2].sort_values("t_abs_c")
        ax.errorbar(sub["t_abs_c"], sub["rad_corr"],
                    xerr=0.5*(sub["t_abs_hi"] - sub["t_abs_lo"]),
                    yerr=sub["rad_corr_sigma"].fillna(0),
                    fmt=markers_q2[i % len(markers_q2)],
                    color=q2_colors[i], label=q2_label(iQ2),
                    markersize=5, linewidth=1.2, capsize=3)
        ax.fill_between(sub["t_abs_c"],
                        sub["rad_corr"] - sub["rad_corr_sigma"].fillna(0),
                        sub["rad_corr"] + sub["rad_corr_sigma"].fillna(0),
                        color=q2_colors[i], alpha=0.12)
    for e in sorted(df["t_abs_lo"].unique()):
        ax.axvline(e, color="gray", lw=0.6, ls=":")
    ax.axhline(1.0, color="black", lw=0.8, ls="--", alpha=0.5)
    ax.set_xlabel("$|t|$ (GeV$^2$)")
    ax.set_ylabel("Radiative Correction $\\delta_{RC}$")
    ax.set_title("$\\delta_{RC}$ vs $|t|$")
    ax.set_ylim(ylo, yhi)
    ax.legend(loc="best", framealpha=0.7, ncol=1)

# ── Plot 3: δ_RC vs Q² ───────────────────────────────────────────────────────
def plot_vs_Q2(ax):
    for j, itp in enumerate(tp_bins):
        sub = df[df["itp"] == itp].sort_values("Q2_c")
        ax.errorbar(sub["Q2_c"], sub["rad_corr"],
                    xerr=0.5*(sub["Q2_hi"] - sub["Q2_lo"]),
                    yerr=sub["rad_corr_sigma"].fillna(0),
                    fmt=markers_tp[j % len(markers_tp)],
                    color=tp_colors[j], label=tp_label(itp),
                    markersize=5, linewidth=1.2, capsize=3)
        ax.fill_between(sub["Q2_c"],
                        sub["rad_corr"] - sub["rad_corr_sigma"].fillna(0),
                        sub["rad_corr"] + sub["rad_corr_sigma"].fillna(0),
                        color=tp_colors[j], alpha=0.10)
    for e in sorted(df["Q2_lo"].unique()):
        ax.axvline(e, color="gray", lw=0.6, ls=":")
    ax.axhline(1.0, color="black", lw=0.8, ls="--", alpha=0.5)
    ax.set_xlabel("$Q^2$ (GeV$^2$)")
    ax.set_ylabel("Radiative Correction $\\delta_{RC}$")
    ax.set_title("$\\delta_{RC}$ vs $Q^2$")
    ax.set_ylim(ylo, yhi)
    # Two-column legend for many t' bins
    ax.legend(loc="best", framealpha=0.7, ncol=2, fontsize=7)

# ── Plot 4: 2D heatmap ───────────────────────────────────────────────────────
def plot_heatmap(ax):
    # Build 2D array: rows = Q² bins, cols = t' bins
    grid = np.full((nQ2, nTp), np.nan)
    for _, row in df.iterrows():
        i = q2_bins.index(int(row["iQ2"]))
        j = tp_bins.index(int(row["itp"]))
        grid[i, j] = row["rad_corr"]

    im = ax.imshow(grid, aspect="auto", origin="lower",
                   cmap="RdYlGn", vmin=ylo, vmax=yhi,
                   extent=[-0.5, nTp - 0.5, -0.5, nQ2 - 0.5])
    plt.colorbar(im, ax=ax, label="$\\delta_{RC}$", pad=0.02)

    # Annotate cells
    for i in range(nQ2):
        for j in range(nTp):
            if not np.isnan(grid[i, j]):
                ax.text(j, i, f"{grid[i, j]:.3f}",
                        ha="center", va="center", fontsize=7,
                        color="black" if 0.85 < grid[i, j] < 1.15 else "white")

    # Axis labels: t' bin centre for x, Q² bin centre for y
    tp_centres = [df[df["itp"] == itp]["tprime_c"].iloc[0] for itp in tp_bins]
    q2_centres = [df[df["iQ2"] == iQ2]["Q2_c"].iloc[0]     for iQ2 in q2_bins]
    ax.set_xticks(range(nTp))
    ax.set_xticklabels([f"{v:.2f}" for v in tp_centres], rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(nQ2))
    ax.set_yticklabels([f"{v:.2f}" for v in q2_centres], fontsize=8)
    ax.set_xlabel("$t'$ bin centre (GeV$^2$)")
    ax.set_ylabel("$Q^2$ bin centre (GeV$^2$)")
    ax.set_title("$\\delta_{RC}(Q^2,\\ t')$ heatmap")

# ── Individual PDFs ───────────────────────────────────────────────────────────
for fname, plot_fn, figsize in [
    ("rad_corr_vs_tprime.pdf", plot_vs_tprime, (7, 5)),
    ("rad_corr_vs_tabs.pdf",   plot_vs_tabs,   (7, 5)),
    ("rad_corr_vs_Q2.pdf",     plot_vs_Q2,     (7, 5)),
    ("rad_corr_heatmap.pdf",   plot_heatmap,   (8, 5)),
]:
    fig, ax = plt.subplots(figsize=figsize)
    plot_fn(ax)
    if kin_label:
        fig.text(0.5, 0.01, kin_label, ha="center", fontsize=9, style="italic")
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    p = out_path(fname)
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {p}")

# ── Summary: all four panels on one page ─────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
plot_vs_tprime(axes[0, 0])
plot_vs_tabs  (axes[0, 1])
plot_vs_Q2    (axes[1, 0])
plot_heatmap  (axes[1, 1])
fig.suptitle(f"MDIFFRAD Radiative Corrections  —  {kin_label}", fontsize=12)
fig.tight_layout()
p = out_path("rad_corr_summary.pdf")
fig.savefig(p, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {p}")

print("\n[plot_rad_corrections] Done.")
