#!/usr/bin/env python3
"""Generate every Fig. 2 plot directly from source tables.

Overview, module prevalence, selected biomass issues and per-model QC burden.
ATP and NADH profiles are now part of Fig. 3.
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np
import pandas as pd

from panel_style import (
    DATA_DIR,
    DB_COLORS,
    DB_ORDER,
    GRID,
    INK,
    apply_style,
    centered_title,
    publication_axes,
    save_panel,
)


EXPECTED_COUNTS = {"CarveMe": 5585, "BiGG": 73, "ModelSEED": 79, "VMH": 806}
MODULES = [
    ("Reducing equivalents", "Reducing\nequivalents"),
    ("Energy", "Energy"),
    ("Net metabolite production", "Net metabolite\nproduction"),
    ("Yield", "Yield"),
    ("Biomass", "Biomass"),
]
BIOMASS_FAILURES = [
    ("Missing biomass reaction", "Missing biomass\nreaction"),
    ("Growth infeasible", "Growth\ninfeasibility"),
    ("Not normalized to 1 gDW", "1 gDW\nnormalization"),
]


def vector_heatmap(ax, values, cmap):
    norm = Normalize(vmin=0, vmax=100)
    for row in range(values.shape[0]):
        for col in range(values.shape[1]):
            value = float(values[row, col])
            ax.add_patch(
                plt.Rectangle(
                    (col - 0.5, row - 0.5),
                    1,
                    1,
                    facecolor=cmap(norm(value)),
                    edgecolor="white",
                    linewidth=1.4,
                )
            )
            ax.text(
                col,
                row,
                f"{value:.1f}%",
                ha="center",
                va="center",
                fontsize=8.5,
                color="white" if value >= 58 else INK,
                fontweight="bold",
            )
    ax.set_xlim(-0.5, values.shape[1] - 0.5)
    ax.set_ylim(values.shape[0] - 0.5, -0.5)
    scalar = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar.set_array(values)
    return scalar


def style_heatmap_axes(ax, row_labels, column_labels):
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels)
    ax.set_xticks(np.arange(len(column_labels)))
    ax.set_xticklabels(column_labels)
    ax.tick_params(
        axis="x",
        top=True,
        bottom=False,
        labeltop=True,
        labelbottom=False,
        length=0,
        pad=6,
    )
    ax.tick_params(axis="y", length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_overall() -> None:
    data = pd.read_csv(DATA_DIR / "fig2" / "overall_error_prevalence.csv")
    ordered = data.set_index("database").loc[DB_ORDER].reset_index()
    fig, ax = plt.subplots(figsize=(6.6, 3.7))
    y = np.arange(len(DB_ORDER))
    bars = ax.barh(
        y,
        ordered["error_rate_pct"],
        height=0.58,
        color=[DB_COLORS[x] for x in DB_ORDER],
        alpha=0.92,
    )
    for bar, record in zip(bars, ordered.itertuples()):
        ax.text(
            min(record.error_rate_pct - 2.0, 97.0),
            bar.get_y() + bar.get_height() / 2,
            f"{record.error_rate_pct:.1f}%",
            ha="right",
            va="center",
            color="white",
            fontweight="bold",
        )
    ax.set_yticks(y, DB_ORDER)
    ax.invert_yaxis()
    ax.set_xlim(0, 104)
    ax.set_xlabel("Models with at least one detected error (%)")
    centered_title(ax, "Models with detected errors")
    publication_axes(ax, grid_axis="x")
    fig.tight_layout()
    save_panel(fig, "Fig2", "Fig2a_overall_error_prevalence")


def draw_heatmap(
    source, stem, title, categories, cmap_colors, colorbar_width=0.035
) -> None:
    data = pd.read_csv(DATA_DIR / "fig2" / source)
    values = (
        data.pivot(
            index="database", columns=categories[0], values="model_prevalence_pct"
        )
        .reindex(index=DB_ORDER, columns=[key for key, _ in categories[1]])
        .to_numpy(dtype=float)
    )
    cmap = LinearSegmentedColormap.from_list(stem, cmap_colors)
    fig, ax = plt.subplots(figsize=(6.6, 3.7))
    image = vector_heatmap(ax, values, cmap)
    style_heatmap_axes(
        ax,
        [f"{db}  ($n$={EXPECTED_COUNTS[db]:,})" for db in DB_ORDER],
        [label for _, label in categories[1]],
    )
    centered_title(ax, title, pad=12)
    cbar = fig.colorbar(
        image,
        ax=ax,
        fraction=colorbar_width,
        pad=0.025,
        boundaries=np.linspace(0, 100, 21),
        ticks=[0, 20, 40, 60, 80, 100],
    )
    cbar.ax.set_title("Affected\n(%)", fontsize=8.0, pad=4)
    fig.tight_layout()
    save_panel(fig, "Fig2", stem)


def draw_module_burden() -> None:
    from fig2_module_burden import refresh

    refresh()
    data = pd.read_csv(DATA_DIR / "fig2/module_burden_distribution.csv")
    values = data.pivot(
        index="database", columns="affected_modules", values="model_pct"
    ).loc[DB_ORDER]
    colors = ["#DBDFE3", "#D7E9F2", "#A6CEDD", "#6BAEC6", "#367F9D", "#164B67"]
    fig, ax = plt.subplots(figsize=(6.6, 3.7))
    y = np.arange(len(DB_ORDER))
    left = np.zeros(len(DB_ORDER))
    for k, color in enumerate(colors):
        widths = values[k].to_numpy()
        ax.barh(
            y,
            widths,
            left=left,
            height=0.60,
            color=color,
            edgecolor="white",
            linewidth=0.7,
            label=str(k),
        )
        for row, width in enumerate(widths):
            if width >= 4:
                ax.text(
                    left[row] + width / 2,
                    row,
                    f"{width:.1f}",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="white" if k >= 4 else INK,
                    fontweight="bold",
                )
        left += widths
    ax.set_yticks(y, [f"{db}  ($n$={EXPECTED_COUNTS[db]:,})" for db in DB_ORDER])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Models (%)")
    centered_title(ax, "Number of affected QC modules per model")
    publication_axes(ax, grid_axis="x")
    means = data.groupby("database").mean_affected_modules.first()
    for row, db in enumerate(DB_ORDER):
        ax.text(
            1.025,
            row,
            f"{means[db]:.2f}",
            transform=ax.get_yaxis_transform(),
            va="center",
            fontsize=8.5,
            fontweight="bold",
            clip_on=False,
        )
    ax.text(1.025, 1.05, "Mean", transform=ax.transAxes, fontsize=8, clip_on=False)
    ax.legend(
        title="Affected modules (0 = all pass)",
        ncol=6,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.50, -0.20),
        fontsize=8,
        title_fontsize=8,
        handlelength=1.2,
        columnspacing=1.0,
    )
    fig.tight_layout()
    save_panel(fig, "Fig2", "Fig2d_module_burden_distribution")


def main() -> None:
    apply_style()
    draw_overall()
    draw_heatmap(
        "module_error_prevalence.csv",
        "Fig2b_module_error_prevalence",
        "Affected models by module",
        ("module", MODULES),
        ["#F7FAFC", "#C8E0E7", "#69A6B5", "#20566B"],
    )
    draw_heatmap(
        "biomass_issue_prevalence.csv",
        "Fig2c_biomass_issue_prevalence",
        "Biomass curation findings",
        ("failure_mode", BIOMASS_FAILURES),
        ["#FFF9F0", "#FCD7A1", "#F28A4B", "#B83D2E"],
        colorbar_width=0.040,
    )
    draw_module_burden()
    print("Fig2: wrote 4 standalone panels")


if __name__ == "__main__":
    main()
