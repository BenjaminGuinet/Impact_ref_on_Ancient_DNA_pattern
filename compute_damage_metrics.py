#!/usr/bin/env python3
"""
Compute damage stochasticity and damage intensity from DamageProfiler output.
Optionally plot superimposed damage profiles (5p left / 3p right).

Damage stochasticity:
    Sum of absolute differences in C>T (5p) / G>A (3p) frequencies at each
    position between an observed profile and a 1M-reads reference profile
    at the same divergence level.

Damage intensity:
    5p: mean C>T at first 3 positions  minus  mean of all other mismatches across all positions
    3p: mean G>A at last  3 positions  minus  mean of all other mismatches across all positions

The 3p file path is derived automatically from the 5p path by replacing
'5p_freq_misincorporations' with '3p_freq_misincorporations'.

Usage:
    python3 compute_damage_metrics.py \\
        --sample_5p <5p_file> \\
        --ref_5p    <1M_ref_5p_file> \\
        --stats     <stats_tsv> \\
        --out       metrics.tsv \\
        --plot      damage_plot.pdf
"""

import argparse
import sys
import os
import textwrap
import numpy as np
import pandas as pd


# ── Column definitions ────────────────────────────────────────────────────────
CT_COL       = "C>T"
GA_COL       = "G>A"
ALL_MISMATCH = ["C>T", "G>A", "A>C", "A>G", "A>T",
                "C>A", "C>G", "G>C", "G>T",
                "T>A", "T>C", "T>G"]
OTHER_MM_5P  = [c for c in ALL_MISMATCH if c != "C>T"]
OTHER_MM_3P  = [c for c in ALL_MISMATCH if c != "G>A"]


# ── File helpers ──────────────────────────────────────────────────────────────

def derive_3p_path(path_5p: str) -> str:
    return path_5p.replace("5p_freq_misincorporations", "3p_freq_misincorporations")


def read_damage_profile(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#")
    df.columns = df.columns.str.strip()
    df = df.sort_values("Pos").reset_index(drop=True)
    return df


def read_mapping_stats(path: str) -> dict:
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return {}
    return {str(k): v for k, v in df.iloc[0].to_dict().items()}


# ── Metrics ───────────────────────────────────────────────────────────────────

def damage_intensity_5p(df: pd.DataFrame, n_edge: int = 3) -> dict:
    edge    = df[df["Pos"] < n_edge]
    ct_edge = float(edge[CT_COL].mean())
    bg      = float(df[OTHER_MM_5P].values.mean())
    return {
        "5p_ct_mean_edge"     : round(ct_edge,      6),
        "5p_bg_mean_other_mm" : round(bg,           6),
        "5p_damage_intensity" : round(ct_edge - bg, 6),
    }


def damage_intensity_3p(df: pd.DataFrame, n_edge: int = 3) -> dict:
    edge    = df[df["Pos"] < n_edge]
    ga_edge = float(edge[GA_COL].mean())
    bg      = float(df[OTHER_MM_3P].values.mean())
    return {
        "3p_ga_mean_edge"     : round(ga_edge,      6),
        "3p_bg_mean_other_mm" : round(bg,           6),
        "3p_damage_intensity" : round(ga_edge - bg, 6),
    }


def damage_stochasticity(df_sample: pd.DataFrame,
                         df_ref:    pd.DataFrame,
                         signal_col: str,
                         label: str) -> dict:
    merged = pd.merge(
        df_sample[["Pos", signal_col]].rename(columns={signal_col: "sample"}),
        df_ref   [["Pos", signal_col]].rename(columns={signal_col: "ref"}),
        on="Pos", how="inner"
    )
    if merged.empty:
        return {f"{label}_n_pos_compared": 0,
                f"{label}_stochasticity":  float("nan")}
    abs_diff = (merged["sample"] - merged["ref"]).abs()
    return {
        f"{label}_n_pos_compared" : int(len(merged)),
        f"{label}_stochasticity"  : round(float(abs_diff.sum()), 6),
    }


def nan_stoch(label: str) -> dict:
    return {f"{label}_n_pos_compared": float("nan"),
            f"{label}_stochasticity":  float("nan")}


# ── PyDamage ──────────────────────────────────────────────────────────────────

PYDAMAGE_COLS = [
    "predicted_accuracy",
    "null_model_p0", "null_model_p0_stdev",
    "damage_model_p", "damage_model_p_stdev",
    "damage_model_pmin", "damage_model_pmin_stdev",
    "damage_model_pmax", "damage_model_pmax_stdev",
    "pvalue", "qvalue", "RMSE",
    "nb_reads_aligned", "coverage", "reflen",
]

def read_pydamage(path: str) -> dict:
    """
    Read the first data row of a PyDamage CSV and return a dict of selected
    columns prefixed with 'pydamage_'.
    """
    try:
        df = pd.read_csv(path)
        if df.empty:
            return {f"pydamage_{c}": float("nan") for c in PYDAMAGE_COLS}
        row = df.iloc[0]
        return {
            f"pydamage_{c}": row[c] if c in row.index else float("nan")
            for c in PYDAMAGE_COLS
        }
    except Exception as e:
        print(f"[WARN] Could not read PyDamage file {path}: {e}", file=sys.stderr)
        return {f"pydamage_{c}": float("nan") for c in PYDAMAGE_COLS}


# ── Per-position expanders ────────────────────────────────────────────────────

def _mm_tag(col: str) -> str:
    """Friendly column tag: 'C>T' -> 'CtoT', etc."""
    return col.replace(">", "to")


def per_position_5p(df5: pd.DataFrame) -> dict:
    """
    Return per-position values for the 5p end.
    CtoT-0, CtoT-1, ... for the primary signal.
    other_5p-0, other_5p-1, ... for the mean of all other mismatches.
    Position 0 = first base of read.
    """
    df = df5.sort_values("Pos")
    out = {}
    for _, row in df.iterrows():
        pos = int(row["Pos"])
        out[f"CtoT-{pos}"] = round(float(row[CT_COL]), 6)
        other_mean = float(df.loc[df["Pos"] == row["Pos"], OTHER_MM_5P].values.mean())
        out[f"other_5p-{pos}"] = round(other_mean, 6)
    return out


def per_position_3p(df3: pd.DataFrame) -> dict:
    """
    Return per-position values for the 3p end.
    GtoA-0, GtoA-1, ... for the primary signal.
    other_3p-0, other_3p-1, ... for the mean of all other mismatches.
    Position 0 = last base of read.
    """
    df = df3.sort_values("Pos")
    out = {}
    for _, row in df.iterrows():
        pos = int(row["Pos"])
        out[f"GtoA-{pos}"] = round(float(row[GA_COL]), 6)
        other_mean = float(df.loc[df["Pos"] == row["Pos"], OTHER_MM_3P].values.mean())
        out[f"other_3p-{pos}"] = round(other_mean, 6)
    return out


# backward-compat aliases
def per_position_ct(df5: pd.DataFrame) -> dict:
    return per_position_5p(df5)


def per_position_ga(df3: pd.DataFrame) -> dict:
    return per_position_3p(df3)


# ── Plot ──────────────────────────────────────────────────────────────────────

def plot_damage(df5_sample, df3_sample,
                df5_ref,    df3_ref,
                label_sample, label_ref,
                metrics: dict,
                out_path: str,
                n_reads_sample: int = None,
                n_reads_ref: int = None):
    """
    Two-panel damage plot:
      Left  — 5p end: C>T in colour, all other mismatches in grey
      Right — 3p end: G>A in colour, all other mismatches in grey
    Colour scheme:
      reference (1M reads) -> red
      sample               -> green  (or blue if no ref)
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines
    import matplotlib.ticker as ticker

    GREY       = "#BBBBBB"
    RED        = "#E03030"
    GREEN      = "#2EAA2E"
    GREY_ALPHA = 0.55
    LW_MAIN    = 2.2
    LW_OTHER   = 0.9

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), sharey=False)
    fig.subplots_adjust(wspace=0.35)

    def _plot_panel(ax, df_sample, df_ref, signal_col, other_cols,
                    xlabel, title, neg_x=False):
        def _pos(arr):
            return -(arr + 1) if neg_x else arr

        pos_s = _pos(df_sample["Pos"].values)

        # Grey: all other mismatches (sample)
        for col in other_cols:
            ax.plot(pos_s, df_sample[col].values,
                    color=GREY, lw=LW_OTHER, alpha=GREY_ALPHA, zorder=1)

        # Grey: all other mismatches (ref, if available)
        if df_ref is not None:
            pos_r = _pos(df_ref["Pos"].values)
            for col in other_cols:
                ax.plot(pos_r, df_ref[col].values,
                        color=GREY, lw=LW_OTHER, alpha=GREY_ALPHA, zorder=1)

        # Red: reference signal (1M reads)
        if df_ref is not None:
            pos_r = _pos(df_ref["Pos"].values)
            lbl_ref = "\n".join(textwrap.wrap(label_ref or "", width=35))
            if n_reads_ref is not None:
                lbl_ref += f"\nn={n_reads_ref:,} reads"
            ax.plot(pos_r, df_ref[signal_col].values,
                    color=RED, lw=LW_MAIN, zorder=3, label=lbl_ref)

        # Green: sample signal
        color_sample = GREEN if df_ref is not None else RED
        lbl_sample = "\n".join(textwrap.wrap(label_sample or "", width=35))
        if n_reads_sample is not None:
            lbl_sample += f"\nn={n_reads_sample:,} reads"
        ax.plot(pos_s, df_sample[signal_col].values,
                color=color_sample, lw=LW_MAIN, zorder=4, label=lbl_sample)

        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_ylabel("Substitution frequency", fontsize=10)
        ax.set_title(title, fontsize=11, fontweight="bold")

        x_min, x_max = pos_s.min(), pos_s.max()
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(bottom=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=9)
        ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        ax.legend(fontsize=7.5, frameon=False, labelspacing=0.6)

    # 5p panel
    _plot_panel(
        axes[0], df5_sample, df5_ref,
        signal_col = CT_COL,
        other_cols = OTHER_MM_5P,
        xlabel     = "Position from 5' end",
        title      = "5' end  --  C->T",
        neg_x      = False,
    )

    # 3p panel
    _plot_panel(
        axes[1], df3_sample, df3_ref,
        signal_col = GA_COL,
        other_cols = OTHER_MM_3P,
        xlabel     = "Position from 3' end",
        title      = "3' end  --  G->A",
        neg_x      = True,
    )

    # Metrics annotation box
    lines = []
    if not np.isnan(metrics.get("5p_stochasticity", float("nan"))):
        lines.append(f"5p stochasticity : {metrics['5p_stochasticity']:.4f}")
    if not np.isnan(metrics.get("3p_stochasticity", float("nan"))):
        lines.append(f"3p stochasticity : {metrics['3p_stochasticity']:.4f}")
    lines.append(f"5p damage intensity : {metrics.get('5p_damage_intensity', float('nan')):.4f}")
    lines.append(f"3p damage intensity : {metrics.get('3p_damage_intensity', float('nan')):.4f}")
    if lines:
        fig.text(0.5, -0.04, "   |   ".join(lines),
                 ha="center", va="top", fontsize=8.5,
                 color="#333333",
                 bbox=dict(boxstyle="round,pad=0.4", fc="#F5F5F5",
                           ec="#CCCCCC", lw=0.8))

    # Legend for grey lines
    grey_handle = mlines.Line2D([], [], color=GREY, lw=LW_OTHER,
                                label="Other mismatches")
    fig.legend(handles=[grey_handle], loc="lower right",
               fontsize=8, frameon=False,
               bbox_to_anchor=(0.98, -0.02))

    sample_str = label_sample[:60] + ("..." if len(label_sample) > 60 else "")
    fig.suptitle(sample_str, fontsize=10, y=1.02, color="#444444")

    plt.tight_layout()
    plt.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"[INFO] Plot saved to: {out_path}", file=sys.stderr)


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample_5p", required=True,
                   help="5p misincorporation file for the sample.")
    p.add_argument("--sample_3p", required=False, default=None,
                   help="3p misincorporation file for the sample. "
                        "If omitted, derived from --sample_5p filename.")
    p.add_argument("--ref_5p",    required=False, default=None,
                   help="5p misincorporation file of the 1M-reads reference.")
    p.add_argument("--ref_3p",    required=False, default=None,
                   help="3p misincorporation file of the 1M-reads reference. "
                        "If omitted, derived from --ref_5p filename.")
    p.add_argument("--pydamage",  required=False, default=None,
                   help="PyDamage CSV for this sample (Pydamage_*.csv).")
    p.add_argument("--stats",     required=False, default=None,
                   help="mapping_stats.tsv produced by the pipeline.")
    p.add_argument("--out",       required=False, default=None,
                   help="Output TSV path. If omitted, prints to stdout.")
    p.add_argument("--plot",      required=False, default=None,
                   help="Output plot path (.pdf / .png / .svg). "
                        "If omitted, no plot is produced.")
    p.add_argument("--n_edge",    type=int, default=3,
                   help="Number of edge positions for intensity (default: 3).")
    p.add_argument("--label",     default=None,
                   help="Override sample label (default: derived from filename).")
    p.add_argument("--label_ref", default=None,
                   help="Override reference label (default: derived from ref filename).")
    return p.parse_args()


def main():
    args = parse_args()

    # Labels
    label_sample = args.label or os.path.basename(args.sample_5p) \
        .replace("_5p_freq_misincorporations.txt", "") \
        .replace("_5p_freq_misincorporations.tsv", "")

    label_ref = None
    if args.ref_5p:
        label_ref = args.label_ref or os.path.basename(args.ref_5p) \
            .replace("_5p_freq_misincorporations.txt", "") \
            .replace("_5p_freq_misincorporations.tsv", "")

    row = {"sample": label_sample}

    # Mapping stats
    if args.stats:
        row.update(read_mapping_stats(args.stats))

    # Load profiles
    df5_sample = read_damage_profile(args.sample_5p)

    path_3p_sample = args.sample_3p or derive_3p_path(args.sample_5p)
    if os.path.exists(path_3p_sample):
        df3_sample    = read_damage_profile(path_3p_sample)
        has_3p_sample = True
    else:
        print(f"[WARN] 3p sample file not found: {path_3p_sample}", file=sys.stderr)
        df3_sample    = None
        has_3p_sample = False

    df5_ref = None
    df3_ref = None
    if args.ref_5p:
        df5_ref = read_damage_profile(args.ref_5p)
        path_3p_ref = args.ref_3p or derive_3p_path(args.ref_5p)
        if os.path.exists(path_3p_ref):
            df3_ref = read_damage_profile(path_3p_ref)
        else:
            print(f"[WARN] 3p ref file not found: {path_3p_ref}", file=sys.stderr)

    # Damage intensity
    row.update(damage_intensity_5p(df5_sample, n_edge=args.n_edge))

    if has_3p_sample:
        row.update(damage_intensity_3p(df3_sample, n_edge=args.n_edge))
    else:
        row.update({"3p_ga_mean_edge"     : float("nan"),
                    "3p_bg_mean_other_mm" : float("nan"),
                    "3p_damage_intensity" : float("nan")})

    # Damage stochasticity
    if df5_ref is not None:
        row.update(damage_stochasticity(df5_sample, df5_ref, CT_COL, "5p"))
        if has_3p_sample and df3_ref is not None:
            row.update(damage_stochasticity(df3_sample, df3_ref, GA_COL, "3p"))
        else:
            row.update(nan_stoch("3p"))
    else:
        print("[INFO] --ref_5p not provided -> stochasticity set to NaN.",
              file=sys.stderr)
        row.update(nan_stoch("5p"))
        row.update(nan_stoch("3p"))

    # PyDamage metrics
    if args.pydamage:
        row.update(read_pydamage(args.pydamage))
    else:
        row.update({f"pydamage_{c}": float("nan") for c in PYDAMAGE_COLS})

    # Per-position mismatches (5p: all, 3p: all)
    row.update(per_position_5p(df5_sample))
    if has_3p_sample:
        row.update(per_position_3p(df3_sample))

    # Print pretty summary
    print("\n" + "=" * 70)
    print(f"  Damage metrics  --  {label_sample}")
    print("=" * 70)
    for col, val in row.items():
        print(f"  {col:<42} {val}")
    print("=" * 70 + "\n")

    # TSV output
    df_out = pd.DataFrame([row])
    if args.out:
        df_out.to_csv(args.out, sep="\t", index=False)
        print(f"[INFO] Table written to: {args.out}", file=sys.stderr)
    else:
        sys.stdout.write(df_out.to_csv(sep="\t", index=False))

    # Plot
    if args.plot:
        if not has_3p_sample:
            print("[WARN] Cannot plot: 3p file missing for sample.", file=sys.stderr)
        else:
            n_reads_sample = int(row["n_subsampled"]) if "n_subsampled" in row and str(row["n_subsampled"]).isdigit() else None
            try:
                n_reads_sample = int(float(row.get("n_subsampled", 0))) or None
            except (TypeError, ValueError):
                n_reads_sample = None
            n_reads_ref = 1_000_000 if df5_ref is not None else None
            plot_damage(
                df5_sample     = df5_sample,
                df3_sample     = df3_sample,
                df5_ref        = df5_ref,
                df3_ref        = df3_ref,
                label_sample   = label_sample,
                label_ref      = label_ref or "reference (1M reads)",
                metrics        = row,
                out_path       = args.plot,
                n_reads_sample = n_reads_sample,
                n_reads_ref    = n_reads_ref,
            )


if __name__ == "__main__":
    main()