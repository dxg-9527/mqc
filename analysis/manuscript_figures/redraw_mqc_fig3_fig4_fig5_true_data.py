#!/usr/bin/env python3
"""Redraw MQC Figure 3/4/5 panels from curated source data.

The script uses the curated workbook produced for the manuscript and writes
publication-style PNG/PDF figures plus CSV summaries for cross-checking.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


PAPER_DIR = Path("/home/dengxg/project/mqc/others/定量分析/result/论文")
WORKBOOK = PAPER_DIR / "MQC图2-5原始数据整理.xlsx"


mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
    }
)


DB_ORDER = ["CarveMe", "BiGG", "ModelSEED", "VMH"]
DB_COLORS = {
    "CarveMe": "#E6862E",
    "BiGG": "#2F73B7",
    "ModelSEED": "#54A65B",
    "VMH": "#D85863",
}
STAGE_COLORS = {"Initial": "#2369B3", "Middle": "#F28E2B", "Final": "#2E8B43"}
INK = "#20262E"
MUTED = "#68727D"
GRID = "#D7DCE2"


def save_all(fig: mpl.figure.Figure, stem: str) -> None:
    # PNG + PDF are the manuscript-facing deliverables here.  SVG export becomes
    # unnecessarily fragile for these dense point layers on this machine.
    for suffix in [".png", ".pdf"]:
        kwargs = {"bbox_inches": "tight", "pad_inches": 0.04}
        if suffix == ".png":
            kwargs["dpi"] = 600
        fig.savefig(PAPER_DIR / f"{stem}{suffix}", **kwargs)
    plt.close(fig)


def pct_in_range(series: pd.Series, lo: float, hi: float) -> float:
    s = pd.to_numeric(series, errors="coerce").dropna()
    if len(s) == 0:
        return np.nan
    return float(((s >= lo) & (s <= hi)).mean() * 100)


def fmt_compact(value: float) -> str:
    """Compact numeric labels that stay legible in small manuscript panels."""
    if not np.isfinite(value):
        return "NA"
    if abs(value) >= 100:
        return f"{value:.0f}"
    if abs(value) >= 10:
        return f"{value:.1f}"
    if abs(value) >= 1:
        return f"{value:.2f}"
    return f"{value:.2g}"


def export_fig3_data() -> pd.DataFrame:
    energy = pd.read_excel(WORKBOOK, sheet_name="图3b_能量原始值")
    redox = pd.read_excel(WORKBOOK, sheet_name="图3c_还原力原始值")
    energy["panel"] = "b_energy"
    redox["panel"] = "c_redox"
    long_df = pd.concat([energy, redox], ignore_index=True)
    long_df.to_csv(PAPER_DIR / "MQC图3bc_true_data_long.csv", index=False)

    summary = (
        long_df.groupby(["panel", "metric", "database", "stage"], observed=True)["value"]
        .agg(n="count", median="median", mean="mean", q25=lambda s: s.quantile(0.25), q75=lambda s: s.quantile(0.75), max="max")
        .reset_index()
    )
    wide = summary.pivot_table(index=["panel", "metric", "database"], columns="stage", values=["n", "median", "mean", "q25", "q75", "max"])
    wide.columns = [f"{a}_{b}".lower() for a, b in wide.columns]
    wide = wide.reset_index()
    wide["median_reduction_pct"] = np.where(
        wide["median_initial"] > 0,
        (1 - wide["median_final"] / wide["median_initial"]) * 100,
        np.nan,
    )
    wide.to_csv(PAPER_DIR / "MQC图3bc_true_data_summary.csv", index=False)
    return wide


def draw_paired_median_grid(
    fig: mpl.figure.Figure,
    outer_grid: mpl.gridspec.GridSpec,
    summary: pd.DataFrame,
    panel: str,
    metrics: list[str],
    panel_label: str,
    title: str,
    cmap_name: str,
) -> None:
    sub = outer_grid.subgridspec(len(metrics), len(DB_ORDER), wspace=0.18, hspace=0.30)
    panel_df = summary[summary["panel"] == panel].copy()
    norm = Normalize(vmin=0, vmax=np.nanpercentile(panel_df["median_initial"], 95))
    cmap = mpl.colormaps[cmap_name]
    axes = []

    for r, metric in enumerate(metrics):
        for c, db in enumerate(DB_ORDER):
            ax = fig.add_subplot(sub[r, c])
            axes.append(ax)
            row = panel_df[(panel_df["metric"] == metric) & (panel_df["database"] == db)]
            if row.empty:
                ax.axis("off")
                continue
            rec = row.iloc[0]
            initial = float(rec["median_initial"])
            final = float(rec["median_final"])
            bg = cmap(norm(initial))
            ax.set_facecolor((*bg[:3], 0.18))
            ax.plot([0, 1], [initial, final], color="#7A7F85", lw=1.1, zorder=1)
            ax.scatter([0], [initial], s=28, color="#1E63B5", edgecolor="white", linewidth=0.5, zorder=3)
            ax.scatter([1], [final], s=28, color="#2B8C3E", edgecolor="white", linewidth=0.5, zorder=3)
            ax.set_yscale("log")
            ax.set_ylim(0.45, 1300)
            ax.set_xlim(-0.35, 1.35)
            ax.set_xticks([0, 1])
            if r == len(metrics) - 1:
                ax.set_xticklabels(["Initial", "Final"], fontsize=6.4)
            else:
                ax.set_xticklabels([])
                ax.tick_params(axis="x", length=0)
            ax.grid(axis="y", color=GRID, lw=0.6, ls="--", alpha=0.75)
            ax.text(
                0.05,
                0.94,
                f"{fmt_compact(initial)} -> {fmt_compact(final)}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=5.8,
                color=INK,
            )
            ax.text(
                0.98,
                0.06,
                f"{rec['median_reduction_pct']:.1f}%\nreduced",
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=5.8,
                color=MUTED,
            )
            if r == 0:
                ax.set_title(db, fontsize=7.8, fontweight="bold", color=DB_COLORS[db], pad=3)
            if c == 0:
                ax.set_ylabel(metric, fontsize=7.8, fontweight="bold", rotation=0, labelpad=22, va="center")
                ax.tick_params(axis="y", labelsize=6.2)
            else:
                ax.set_yticklabels([])
                ax.tick_params(axis="y", length=0)

    left = axes[0].get_position().x0
    top = max(ax.get_position().y1 for ax in axes)
    fig.text(left - 0.052, top + 0.055, panel_label, fontsize=18, fontweight="bold", va="top")
    fig.text(left - 0.010, top + 0.055, title, fontsize=10.8, fontweight="bold", va="top")


def draw_fig3() -> None:
    summary = export_fig3_data()
    fig = plt.figure(figsize=(14.6, 7.2))
    outer = fig.add_gridspec(1, 2, left=0.055, right=0.985, top=0.79, bottom=0.13, wspace=0.11)
    draw_paired_median_grid(
        fig,
        outer[0],
        summary,
        "b_energy",
        ["ATP", "GTP", "CTP"],
        "b",
        "Suppression of unbounded energy production",
        "Blues",
    )
    draw_paired_median_grid(
        fig,
        outer[1],
        summary,
        "c_redox",
        ["NADH", "NADPH", "Q8H2"],
        "c",
        "Suppression of unbounded redox production",
        "Greens",
    )
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#1E63B5", markeredgecolor="white", markersize=7, label="Initial median"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#2B8C3E", markeredgecolor="white", markersize=7, label="Final median"),
        Line2D([0], [0], color="#7A7F85", lw=1.2, label="Median shift"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8.2, bbox_to_anchor=(0.52, 0.035))
    fig.text(0.055, 0.93, "MQC reduces spurious energy and redox production across GEM resources", fontsize=13.5, fontweight="bold")
    fig.text(0.055, 0.897, "Cells show model-level medians; y-axes use a log scale to retain both extreme initial values and near-zero final values.", fontsize=7.8, color=MUTED)
    save_all(fig, "MQC图3bc_true_data_restyle")


def export_fig4_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    growth = pd.read_excel(WORKBOOK, sheet_name="图4a_生长分布")
    cme = pd.read_excel(WORKBOOK, sheet_name="图4b_CarveMe折线")
    vmh = pd.read_excel(WORKBOOK, sheet_name="图4b_VMH折线")
    line_df = pd.concat([cme, vmh], ignore_index=True)
    growth.to_csv(PAPER_DIR / "MQC图4a_true_data_long.csv", index=False)
    line_df.to_csv(PAPER_DIR / "MQC图4b_true_data_long.csv", index=False)

    summary = (
        growth.groupby(["stage", "database"], observed=True)["growth_rate_h1"]
        .agg(n="count", median="median", mean="mean", q25=lambda s: s.quantile(0.25), q75=lambda s: s.quantile(0.75), max="max")
        .reset_index()
    )
    summary["pct_0_1"] = [
        pct_in_range(growth[(growth["stage"] == r.stage) & (growth["database"] == r.database)]["growth_rate_h1"], 0, 1)
        for r in summary.itertuples()
    ]
    line_summary = []
    for db, df in line_df.groupby("database"):
        line_summary.append(
            {
                "database": db,
                "n": len(df),
                "initial_median": df["initial_growth_h1"].median(),
                "final_median": df["final_growth_h1"].median(),
                "initial_mean": df["initial_growth_h1"].mean(),
                "final_mean": df["final_growth_h1"].mean(),
                "initial_pct_0_1": pct_in_range(df["initial_growth_h1"], 0, 1),
                "final_pct_0_1": pct_in_range(df["final_growth_h1"], 0, 1),
            }
        )
    line_summary = pd.DataFrame(line_summary)
    summary.to_csv(PAPER_DIR / "MQC图4a_growth_distribution_summary.csv", index=False)
    line_summary.to_csv(PAPER_DIR / "MQC图4b_growth_ecdf_summary.csv", index=False)
    return growth, line_df


def violin_with_points(ax, values_by_db: dict[str, np.ndarray], stage: str, y_max: float | None = None) -> None:
    positions = np.arange(1, len(DB_ORDER) + 1)
    data = []
    for db in DB_ORDER:
        vals = np.asarray(values_by_db.get(db, []), dtype=float)
        vals = vals[np.isfinite(vals)]
        if y_max is not None:
            vals = np.clip(vals, 0, y_max)
        data.append(vals)
    viol = ax.violinplot(data, positions=positions, widths=0.78, showmeans=False, showmedians=False, showextrema=False)
    for body, db in zip(viol["bodies"], DB_ORDER):
        body.set_facecolor(DB_COLORS[db])
        body.set_edgecolor(DB_COLORS[db])
        body.set_alpha(0.25)
        body.set_linewidth(0.8)
    rng = np.random.default_rng(7)
    for pos, db, vals in zip(positions, DB_ORDER, data):
        if len(vals) == 0:
            continue
        sample = vals if len(vals) <= 700 else rng.choice(vals, size=700, replace=False)
        jitter = rng.normal(0, 0.045, len(sample))
        ax.scatter(np.full(len(sample), pos) + jitter, sample, s=4.2, color=DB_COLORS[db], alpha=0.42, linewidths=0, rasterized=True)
        q1, med, q3 = np.quantile(vals, [0.25, 0.5, 0.75])
        ax.plot([pos - 0.16, pos + 0.16], [med, med], color=INK, lw=1.1)
        ax.add_patch(plt.Rectangle((pos - 0.14, q1), 0.28, q3 - q1, ec=INK, fc="white", lw=0.8, alpha=0.70))
    ax.set_xticks(positions)
    ax.set_xticklabels(DB_ORDER, rotation=0)
    ax.set_title(stage, fontsize=9.8, fontweight="bold")
    ax.grid(axis="y", color=GRID, lw=0.6, ls="--", alpha=0.8)


def plot_ecdf(ax, values: np.ndarray, color: str, label: str, ls: str = "-") -> None:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    # A small display floor keeps true zero-growth models visible on log axes.
    vals = np.clip(vals, 1e-3, None)
    vals = np.sort(vals)
    y = np.arange(1, len(vals) + 1) / len(vals)
    ax.plot(vals, y, color=color, lw=2.0, ls=ls, label=label)


def draw_fig4() -> None:
    growth, line_df = export_fig4_data()
    fig = plt.figure(figsize=(14.4, 8.8))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.08, 0.92], left=0.065, right=0.985, top=0.88, bottom=0.12, wspace=0.16, hspace=0.38)

    y_limits = {"Initial": (0, 230), "Middle": (0, 2.0), "Final": (0, 2.2)}
    for i, stage in enumerate(["Initial", "Middle", "Final"]):
        ax = fig.add_subplot(gs[0, i])
        vals = {db: growth[(growth["stage"] == stage) & (growth["database"] == db)]["growth_rate_h1"].dropna().to_numpy() for db in DB_ORDER}
        violin_with_points(ax, vals, stage, y_max=y_limits[stage][1])
        ax.set_ylim(*y_limits[stage])
        if i == 0:
            ax.set_ylabel("Output (1/h)")
        else:
            ax.set_ylabel("")
        if stage == "Middle":
            ax.text(0.98, 0.95, "Values >2 clipped", transform=ax.transAxes, fontsize=6.4, color=MUTED, va="top", ha="right")
        if stage == "Final":
            ax.axhspan(0, 1, color="#DDEFE2", alpha=0.65, zorder=0)
            ax.text(0.03, 0.47, "0-1 1/h\nbiologically plausible window", transform=ax.transAxes, fontsize=6.3, color="#2E8B43", va="center")

    for j, db in enumerate(["CarveMe", "VMH"]):
        ax = fig.add_subplot(gs[1, j + 0])
        df = line_df[line_df["database"] == db]
        color = DB_COLORS[db]
        ax.axvspan(1e-3, 1, color="#DDEFE2", alpha=0.55, zorder=0)
        plot_ecdf(ax, df["initial_growth_h1"].to_numpy(), color=color, label="Initial", ls="-")
        plot_ecdf(ax, df["final_growth_h1"].to_numpy(), color=color, label="Final", ls="--")
        ax.axvline(1, color="#2E8B43", lw=1.0, ls=":")
        ax.set_xscale("log")
        ax.set_xlim(1e-3, 250)
        ax.set_ylim(0, 1.02)
        ax.grid(color=GRID, lw=0.6, ls="--", alpha=0.8)
        ax.set_title(f"{db} growth-rate ECDF", fontsize=9.8, fontweight="bold", color=color)
        ax.set_xlabel("Growth rate (1/h, log scale)")
        if j == 0:
            ax.set_ylabel("Cumulative fraction")
        pct = pct_in_range(df["final_growth_h1"], 0, 1)
        med_i = df["initial_growth_h1"].median()
        med_f = df["final_growth_h1"].median()
        ax.text(
            0.03,
            0.18,
            f"Final in 0-1 1/h: {pct:.1f}%\nMedian: {med_i:.1f} -> {med_f:.2f}",
            transform=ax.transAxes,
            fontsize=8.0,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=GRID, alpha=0.92),
        )
        ax.legend(frameon=False, loc="lower right", fontsize=8)

    ax_summary = fig.add_subplot(gs[1, 2])
    final = growth[growth["stage"] == "Final"].copy()
    pct_rows = []
    for db in DB_ORDER:
        vals = final[final["database"] == db]["growth_rate_h1"]
        pct_rows.append((db, pct_in_range(vals, 0, 1), float(pd.to_numeric(vals, errors="coerce").median())))
    y_pos = np.arange(len(pct_rows))
    pcts = [row[1] for row in pct_rows]
    colors = [DB_COLORS[row[0]] for row in pct_rows]
    ax_summary.axvspan(80, 100, color="#DDEFE2", alpha=0.55, zorder=0)
    ax_summary.barh(y_pos, pcts, color=colors, alpha=0.82, height=0.55)
    for y, (db, pct, med) in zip(y_pos, pct_rows):
        ax_summary.text(min(pct + 1.6, 99.2), y, f"{pct:.1f}%", va="center", ha="left", fontsize=8.0, color=INK)
        ax_summary.text(2.0, y, f"median {med:.2f}", va="center", ha="left", fontsize=6.5, color="white", fontweight="bold")
    ax_summary.set_yticks(y_pos)
    ax_summary.set_yticklabels([row[0] for row in pct_rows], fontsize=8.0)
    ax_summary.invert_yaxis()
    ax_summary.set_xlim(0, 100)
    ax_summary.set_xlabel("Final models within 0-1 1/h (%)")
    ax_summary.set_title("Final-range summary", fontsize=9.8, fontweight="bold")
    ax_summary.grid(axis="x", color=GRID, lw=0.6, ls="--", alpha=0.8)

    fig.text(0.065, 0.94, "Biomass-yield distributions during iterative MQC curation", fontsize=13.5, fontweight="bold")
    fig.text(0.065, 0.90, "a", fontsize=18, fontweight="bold")
    fig.text(0.065, 0.475, "b", fontsize=18, fontweight="bold")
    save_all(fig, "MQC图4ab_true_data_restyle")


@dataclass
class ParityDataset:
    panel: str
    title: str
    sheet: str
    category_col: str
    x_col: str
    initial_col: str
    final_col: str
    y_label: str
    csv_name: str


FIG5_DATASETS = [
    ParityDataset(
        panel="a",
        title="B. subtilis substrate growth",
        sheet="图5a_Bacillus",
        category_col="condition",
        x_col="experimental_growth_h1",
        initial_col="initial_model_growth_h1",
        final_col="control_model_growth_h1",
        y_label="Simulated growth rate (h$^{-1}$)",
        csv_name="MQC图5a_Bacillus_true_data.csv",
    ),
    ParityDataset(
        panel="b",
        title="E. coli uptake/secretion fluxes",
        sheet="图5b_Ecoli通量",
        category_col="measurement",
        x_col="experimental_flux",
        initial_col="initial_model_flux",
        final_col="control_model_flux",
        y_label="Simulated flux",
        csv_name="MQC图5b_Ecoli_flux_true_data.csv",
    ),
    ParityDataset(
        panel="c",
        title="E. coli substrate growth",
        sheet="图5c_Ecoli生长",
        category_col="substrate",
        x_col="experimental_growth_h1",
        initial_col="initial_model_growth_h1",
        final_col="control_model_growth_h1",
        y_label="Simulated growth rate (h$^{-1}$)",
        csv_name="MQC图5c_Ecoli_growth_true_data.csv",
    ),
]


def color_map_for_categories(categories: list[str]) -> dict[str, str]:
    base = [
        "#2F73B7",
        "#D85863",
        "#E6862E",
        "#54A65B",
        "#8E67C6",
        "#16A6A8",
        "#C49A00",
        "#7F5539",
        "#6C8EBF",
        "#E377C2",
        "#6BAED6",
        "#B15928",
    ]
    return {cat: base[i % len(base)] for i, cat in enumerate(categories)}


def mae(a: pd.Series, b: pd.Series) -> float:
    return float(np.mean(np.abs(pd.to_numeric(a, errors="coerce") - pd.to_numeric(b, errors="coerce"))))


def corr_rank(a: pd.Series, b: pd.Series) -> float:
    x = pd.to_numeric(a, errors="coerce")
    y = pd.to_numeric(b, errors="coerce")
    ok = x.notna() & y.notna()
    if ok.sum() < 3:
        return np.nan
    return float(x[ok].rank().corr(y[ok].rank()))


def draw_parity_row(fig: mpl.figure.Figure, grid, ds: ParityDataset, row_idx: int, all_summary: list[dict]) -> None:
    df = pd.read_excel(WORKBOOK, sheet_name=ds.sheet)
    df.to_csv(PAPER_DIR / ds.csv_name, index=False)
    categories = [str(x) for x in df[ds.category_col].dropna().unique()]
    cmap = color_map_for_categories(categories)
    sub = grid.subgridspec(1, len(DB_ORDER), wspace=0.26)
    row_axes = []

    for c, db in enumerate(DB_ORDER):
        ax = fig.add_subplot(sub[0, c])
        row_axes.append(ax)
        ddb = df[df["database"] == db].copy()
        x = pd.to_numeric(ddb[ds.x_col], errors="coerce")
        yi = pd.to_numeric(ddb[ds.initial_col], errors="coerce")
        yf = pd.to_numeric(ddb[ds.final_col], errors="coerce")
        lim_max = np.nanmax([x.max(), yi.max(), yf.max()])
        lim_min = 0
        lim_pad = 0.08 * lim_max if lim_max > 0 else 0.1
        ax.plot([lim_min, lim_max + lim_pad], [lim_min, lim_max + lim_pad], color="#7A7A7A", lw=1.0, zorder=0)

        for _, rec in ddb.iterrows():
            cat = str(rec[ds.category_col])
            color = cmap.get(cat, "#555555")
            xv = float(rec[ds.x_col])
            iv = float(rec[ds.initial_col])
            fv = float(rec[ds.final_col])
            ax.plot([xv, xv], [iv, fv], color="#A2A8AE", lw=0.8, alpha=0.75, zorder=1)
            ax.scatter([xv], [iv], s=30, facecolor="white", edgecolor=color, linewidth=1.0, marker="o", zorder=3)
            ax.scatter([xv], [fv], s=34, facecolor=color, edgecolor="white", linewidth=0.45, marker="o", zorder=4)

        before = mae(x, yi)
        after = mae(x, yf)
        improvement = (before - after) / before * 100 if before > 0 else np.nan
        rho_after = corr_rank(x, yf)
        all_summary.append(
            {
                "panel": ds.panel,
                "dataset": ds.title,
                "database": db,
                "n": len(ddb),
                "mae_initial": before,
                "mae_after_mqc": after,
                "mae_reduction_pct": improvement,
                "spearman_after": rho_after,
            }
        )

        ax.set_xlim(lim_min, lim_max + lim_pad)
        ax.set_ylim(lim_min, lim_max + lim_pad)
        ax.grid(color=GRID, lw=0.55, ls="--", alpha=0.70)
        ax.set_title(db, fontsize=8.6, fontweight="bold", color=DB_COLORS[db], pad=2)
        if c == 0:
            ax.set_ylabel(ds.y_label, fontsize=8.0)
        else:
            ax.set_ylabel("")
        ax.set_xlabel("Experimental value", fontsize=7.5)
        ax.text(
            0.04,
            0.94,
            f"MAE {before:.2g}->{after:.2g}\n" + (f"{improvement:.0f}% lower" if np.isfinite(improvement) else ""),
            transform=ax.transAxes,
            va="top",
            fontsize=6.2,
            color=INK,
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec=GRID, alpha=0.92),
        )

    left = row_axes[0].get_position().x0
    top = max(ax.get_position().y1 for ax in row_axes)
    fig.text(left - 0.035, top + 0.044, ds.panel, fontsize=18, fontweight="bold", va="top")
    fig.text(left + 0.000, top + 0.044, ds.title, fontsize=11.0, fontweight="bold", va="top")
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor="#333", markersize=6, label="Initial"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#333", markeredgecolor="white", markersize=6, label="After MQC"),
    ]
    cat_handles = [Line2D([0], [0], marker="o", color="none", markerfacecolor=cmap[cat], markeredgecolor="white", markersize=5, label=cat) for cat in categories[:8]]
    if len(categories) > 8:
        cat_handles.append(Line2D([0], [0], color="none", label=f"+{len(categories)-8} categories"))
    row_axes[-1].legend(handles=handles + cat_handles, loc="center left", bbox_to_anchor=(1.015, 0.5), frameon=False, fontsize=6.0, handletextpad=0.45)


def draw_fig5() -> None:
    fig = plt.figure(figsize=(15.8, 11.8))
    gs = fig.add_gridspec(3, 1, left=0.060, right=0.885, top=0.885, bottom=0.07, hspace=0.48)
    summary_rows: list[dict] = []
    for i, ds in enumerate(FIG5_DATASETS):
        draw_parity_row(fig, gs[i, 0], ds, i, summary_rows)

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(PAPER_DIR / "MQC图5_true_data_mae_summary.csv", index=False)

    fig.text(0.060, 0.965, "MQC improves agreement with experimental phenotypes across three validation datasets", fontsize=13.5, fontweight="bold")
    fig.text(
        0.060,
        0.935,
        "Open points show initial models, filled points show MQC-corrected models; grey vertical segments connect each before/after pair at the same experimental value.",
        fontsize=7.8,
        color=MUTED,
    )
    save_all(fig, "MQC图5_true_data_restyle")


def main() -> None:
    draw_fig3()
    draw_fig4()
    draw_fig5()
    print(PAPER_DIR / "MQC图3bc_true_data_restyle.png")
    print(PAPER_DIR / "MQC图4ab_true_data_restyle.png")
    print(PAPER_DIR / "MQC图5_true_data_restyle.png")


if __name__ == "__main__":
    main()
