#!/usr/bin/env python3
"""Replot archived Fig. 4c values; no metabolic model is optimized."""
from pathlib import Path
import csv
import math
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import font_manager

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "output"
# Retain normalized positions from the accepted 8-inch layout, but replot at
# 9.18 inches to use the space released by removing the main-figure workflow.
# Fonts, marker sizes and horizontal data coordinates are NOT stretched.
FIGSIZE_IN = (3.2, 9.18)
INK = "#20262E"
MUTED = "#68727D"
INITIAL = "#6F7882"
POST = "#009E73"
FONT = "Liberation Sans"


def read(name):
    with (ROOT / "data" / name).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def validate(no, matched):
    assert len(no) == 4 and len(matched) == 3
    assert [(r["objective"], float(r["initial"]), float(r["post_qc"])) for r in no] == [
        ("NADH", 665.802, 0.0),
        ("ATP/CTP/GTP/UTP/ITP", 96.25, 0.0),
        ("Urea", 5.0, 0.0),
        ("Allophanate", 5.0, 0.0),
    ]
    expected = [(520.555556, 18.333333), (20.833333, 14.166667), (0.251959, 0.169333)]
    assert [r["metric"] for r in matched] == [
        "Maximum NADH production",
        "Maximum ATP production",
        "Maximum biomass production",
    ]
    for row, values in zip(matched, expected):
        assert (float(row["initial"]), float(row["post_qc"])) == values
        assert row["carbon_source"] == "D-glucose" and float(row["uptake"]) == 10
        assert row["other_carbon_sources"] == "closed"


def setup():
    for p in (ROOT.parent / "assets" / "liberation-fonts-ttf-2.1.5").glob("*.ttf"):
        font_manager.fontManager.addfont(str(p))
    plt.rcParams.update(
        {
            "font.family": FONT,
            "font.size": 13.5,
            "text.color": INK,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "axes.edgecolor": INK,
            "axes.linewidth": 0.85,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "mathtext.fontset": "custom",
            "mathtext.rm": FONT,
            "mathtext.it": FONT + ":italic",
            "mathtext.bf": FONT + ":bold",
            "mathtext.sf": FONT,
            "mathtext.cal": FONT,
            "mathtext.tt": "Liberation Mono",
        }
    )


def dots(ax, pre, post, y):
    ax.plot([post, pre], [y, y], color="#AEB5BC", lw=1.2)
    ax.scatter(pre, y, s=30, color=INITIAL, edgecolor="white", linewidth=0.6, zorder=3)
    ax.scatter(post, y, s=38, marker="x", color=POST, linewidth=1.5, zorder=4)


def make(no, matched):
    fig = plt.figure(figsize=FIGSIZE_IN, facecolor="white")
    fig.text(
        0.5,
        0.987,
        "Post-MQC validation",
        ha="center",
        va="top",
        fontsize=15.5,
        weight="bold",
    )
    fig.legend(
        handles=[
            Line2D([], [], marker="o", color=INITIAL, ls="", label="Initial iCN718"),
            Line2D([], [], marker="x", color=POST, ls="", label="After MQC"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.52, 0.949),
        frameon=False,
        fontsize=13.5,
        labelspacing=0.35,
        handletextpad=0.35,
        borderaxespad=0,
    )
    fig.text(
        0.5,
        1 - 1.10 / 8,
        "No-carbon objectives",
        ha="center",
        va="top",
        fontsize=14.8,
        weight="bold",
    )
    ax = fig.add_axes([0.365, 1 - 3.10 / 8, 0.60, 1.70 / 8])
    names = ["NADH", "NTPs", "Urea", "Allophanate"]
    audit = []
    for yy, (row, label) in enumerate(zip(no, names)):
        y = 3 - yy
        pre = float(row["initial"])
        post = float(row["post_qc"])
        dots(ax, pre, post, y)
        ax.annotate(
            f"{pre:.3f}",
            (pre, y),
            xytext=(-1, 7),
            textcoords="offset points",
            ha="right",
            va="bottom",
            fontsize=13.5,
            color=INITIAL,
        )
        audit += [
            dict(
                section="no_carbon",
                metric=row["objective"],
                series=series,
                value=value,
                axis="symlog",
                xmin=-0.3,
                xmax=1600,
            )
            for series, value in [("initial", pre), ("post_mqc", post)]
        ]
    ax.set_xscale("symlog", linthresh=1)
    ax.set_xlim(-0.3, 1600)
    ax.set_ylim(-0.5, 3.62)
    ax.set_yticks([3, 2, 1, 0], names, fontsize=13.5)
    ax.set_xticks([0, 10, 1000], ["0", "10", "1000"], fontsize=13.5)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=5)
    ax.grid(axis="x", color="#D7DCE2", lw=0.55)
    ax.set_axisbelow(True)
    fig.text(0.5, 1 - 3.51 / 8, "Flux (symlog)", ha="center", va="top", fontsize=13.5)
    fig.text(
        0.5,
        1 - 3.86 / 8,
        "Matched-carbon re-test",
        ha="center",
        va="top",
        fontsize=14.8,
        weight="bold",
    )
    fig.text(
        0.5,
        1 - 4.18 / 8,
        "Glucose uptake: 10\nOther C sources closed",
        ha="center",
        va="top",
        fontsize=13.5,
        color=MUTED,
        linespacing=1.15,
    )
    configs = [
        ("NADH production", 4.81, 5.38, 600, [0, 300, 600]),
        ("ATP production", 5.89, 6.46, 25, [0, 10, 20]),
        ("Biomass production", 6.97, 7.54, 0.4, [0, 0.2, 0.4]),
    ]
    for row, (title, title_y, bottom, maxv, ticks) in zip(matched, configs):
        ax = fig.add_axes([0.12, 1 - bottom / 8, 0.82, 0.35 / 8])
        pre = float(row["initial"])
        post = float(row["post_qc"])
        dots(ax, pre, post, 0.42)
        ax.set_xlim(0, maxv)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.set_xticks(ticks)
        ax.tick_params(axis="x", labelsize=13.5)
        fig.text(
            0.53,
            1 - title_y / 8,
            title,
            fontsize=14.8,
            weight="bold",
            ha="center",
            va="top",
        )
        # Numeric labels use edge-aware placement; separate close ATP endpoints.
        ax.annotate(
            f"{post:.3f}",
            (post, 0.42),
            xytext=(0, 7),
            textcoords="offset points",
            ha="left" if post < 0.15 * maxv else "right",
            va="bottom",
            fontsize=13.5,
            color=INK,
        )
        ax.annotate(
            f"{pre:.3f}",
            (pre, 0.42),
            xytext=(0, 7),
            textcoords="offset points",
            ha="right" if pre > 0.75 * maxv else "left",
            va="bottom",
            fontsize=13.5,
            color=INITIAL,
        )
        if "Biomass" in title:
            ax.set_xlabel(r"h$^{-1}$", fontsize=13.5, labelpad=2)
        audit += [
            dict(
                section="matched_carbon",
                metric=row["metric"],
                series=series,
                value=value,
                axis="linear",
                xmin=0,
                xmax=maxv,
            )
            for series, value in [("initial", pre), ("post_mqc", post)]
        ]
    fig.canvas.draw()
    artists = [(a, collection) for a in fig.axes for collection in a.collections]
    assert len(artists) == len(audit) == 14
    for row, (a, collection) in zip(audit, artists):
        actual_x, actual_y = map(float, collection.get_offsets()[0])
        assert math.isclose(actual_x, row["value"], rel_tol=0, abs_tol=1e-12)
        normalized = a.transAxes.inverted().transform(
            a.transData.transform((actual_x, actual_y))
        )
        assert 0 <= normalized[0] <= 1 and 0 <= normalized[1] <= 1
        row.update(
            artist_x=actual_x,
            artist_y=actual_y,
            x_axis_fraction=float(normalized[0]),
            y_axis_fraction=float(normalized[1]),
        )
    return fig, audit


def main():
    no = read("no_carbon_objectives.csv")
    matched = read("matched_glucose_retest.csv")
    validate(no, matched)
    setup()
    fig, audit = make(no, matched)
    OUTPUT.mkdir(exist_ok=True)
    stem = OUTPUT / "Fig4_panel_c_post_QC_validation"
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(
            stem.with_suffix("." + ext),
            dpi=600,
            facecolor="white" if ext == "png" else "none",
            transparent=ext != "png",
        )
    with (OUTPUT / "Fig4_panel_c_coordinate_audit.csv").open(
        "w", encoding="utf-8", newline=""
    ) as f:
        writer = csv.DictWriter(f, fieldnames=list(audit[0]))
        writer.writeheader()
        writer.writerows(audit)
    plt.close(fig)
    print("Panel c rebuilt from unchanged CSVs; 3.20 x 9.18 in; 13.5–15.5 pt.")


if __name__ == "__main__":
    main()
