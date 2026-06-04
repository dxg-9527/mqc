#!/usr/bin/env python3
"""Redraw MQC Figure 1 as a reproducible publication-style vector figure."""

from __future__ import annotations

from pathlib import Path
import textwrap

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch


PAPER_DIR = Path("/home/dengxg/project/mqc/others/定量分析/result/论文")
OUT_STEM = PAPER_DIR / "MQC_Figure1_framework_redrawn"


mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.linewidth": 0.8,
    }
)


COLORS = {
    "ink": "#1F2933",
    "muted": "#65717F",
    "panel_edge": "#D7DEE8",
    "light_grid": "#EDF1F5",
    "energy": "#B9403D",
    "energy_light": "#F8DEDC",
    "biomass": "#D29A16",
    "biomass_light": "#FFF0C2",
    "yield": "#2F6FA3",
    "yield_light": "#DDECF8",
    "green": "#2F8C58",
    "green_light": "#DFF0E5",
    "gray": "#F6F8FA",
    "dark_gray": "#4A5561",
    "orange": "#DD7A45",
}


def add_panel(ax, x, y, w, h, label, title, title_size=14):
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.014",
        fc="white",
        ec=COLORS["panel_edge"],
        lw=1.4,
        zorder=0,
    )
    ax.add_patch(box)
    ax.text(x + 0.018 * w, y + h - 0.045 * h, label, fontsize=23, fontweight="bold", va="top", color="black")
    ax.text(
        x + 0.095 * w,
        y + h - 0.052 * h,
        title,
        fontsize=title_size,
        fontweight="bold",
        va="top",
        color="black",
    )
    return box


def rounded_box(ax, x, y, w, h, fc, ec, lw=1.5, radius=0.012, alpha=1.0, zorder=2):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        fc=fc,
        ec=ec,
        lw=lw,
        alpha=alpha,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def arrow(ax, start, end, color=COLORS["dark_gray"], lw=1.8, ms=15, style="-|>", alpha=1.0, connectionstyle="arc3"):
    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle=style,
        mutation_scale=ms,
        color=color,
        lw=lw,
        alpha=alpha,
        shrinkA=0,
        shrinkB=0,
        connectionstyle=connectionstyle,
        zorder=5,
    )
    ax.add_patch(patch)
    return patch


def wrap_text(text, width):
    return "\n".join(textwrap.wrap(text, width=width, break_long_words=False))


def draw_energy_icon(ax, cx, cy, scale=1.0):
    r = 0.038 * scale
    for angle1, angle2 in [(35, 170), (215, 350)]:
        path = MplPath.arc(angle1, angle2)
        verts = path.vertices * r + [cx, cy]
        patch = PathPatch(MplPath(verts, path.codes), fc="none", ec=COLORS["energy"], lw=2.2, zorder=4)
        ax.add_patch(patch)
    arrow(ax, (cx + r * 0.62, cy + r * 0.73), (cx + r * 0.98, cy + r * 0.52), color=COLORS["energy"], lw=0, ms=12)
    arrow(ax, (cx - r * 0.65, cy - r * 0.72), (cx - r * 1.0, cy - r * 0.48), color=COLORS["energy"], lw=0, ms=12)
    ax.text(cx, cy, "E", fontsize=17 * scale, fontweight="bold", color=COLORS["energy"], ha="center", va="center")


def draw_biomass_icon(ax, cx, cy, scale=1.0):
    end = (cx + 0.050 * scale, cy)
    for yy in [-0.035, 0.0, 0.035]:
        node = (cx - 0.06 * scale, cy + yy * scale)
        ax.add_patch(Circle(node, 0.008 * scale, fc="#D9E3D4", ec=COLORS["ink"], lw=1.0, zorder=4))
        arrow(ax, (node[0] + 0.009 * scale, node[1]), (end[0] - 0.014 * scale, end[1]), color="black", lw=1.2, ms=9)
    ax.add_patch(Circle(end, 0.015 * scale, fc="#D9E3D4", ec=COLORS["ink"], lw=1.2, zorder=4))
    ax.plot([cx - 0.012 * scale, cx + 0.020 * scale], [cy + 0.043 * scale, cy + 0.017 * scale], ls="--", lw=2.0, color=COLORS["energy"], zorder=4)
    ax.text(cx + 0.008 * scale, cy + 0.046 * scale, "missing\nprecursor", color=COLORS["energy"], fontsize=6.8 * scale)


def draw_yield_icon(ax, cx, cy, scale=1.0):
    left = (cx - 0.060 * scale, cy)
    right = (cx + 0.060 * scale, cy)
    ax.add_patch(Circle(left, 0.009 * scale, fc="#B9CFE1", ec=COLORS["ink"], lw=1.0, zorder=4))
    ax.add_patch(Circle(right, 0.013 * scale, fc="#B9CFE1", ec=COLORS["ink"], lw=1.0, zorder=4))
    arrow(ax, (left[0] + 0.010 * scale, cy), (right[0] - 0.015 * scale, cy), color=COLORS["yield"], lw=2.0, ms=10)
    ax.text(cx - 0.073 * scale, cy + 0.024 * scale, "Input", fontsize=7.0 * scale, color=COLORS["ink"])
    ax.text(cx + 0.032 * scale, cy + 0.024 * scale, "Output", fontsize=7.0 * scale, color=COLORS["ink"])


def draw_panel_a(ax):
    x, y, w, h = 0.035, 0.102, 0.318, 0.785
    add_panel(ax, x, y, w, h, "a", "Three recurrent quantitative infeasibilities", title_size=13.0)

    intro = "Raw GEMs can permit quantitatively impossible fluxes that distort phenotype prediction, flux analysis, and pathway evaluation."
    ax.text(x + 0.050 * w, y + h - 0.118 * h, wrap_text(intro, 62), fontsize=7.3, color=COLORS["muted"], va="top")

    cards = [
        (
            "Unbounded energy/redox production",
            "Closed uptake still permits ATP, XTP, or reducing equivalents to be produced from nothing.",
            COLORS["energy"],
            COLORS["energy_light"],
            draw_energy_icon,
            r"$v_{\mathrm{uptake}}=0;\quad v_{\mathrm{ATP/redox}}>0$",
            0.155,
            0.68,
        ),
        (
            "Biomass/growth infeasibility",
            "Growth or biomass components fail because required precursors or biomass definitions are inconsistent.",
            COLORS["biomass"],
            COLORS["biomass_light"],
            draw_biomass_icon,
            r"$v_{\mathrm{growth}}=0$ or biomass $\ne 1\,\mathrm{gDW}$",
            0.175,
            0.60,
        ),
        (
            "Yield beyond theoretical maximum",
            "Model-predicted product yield exceeds formula-based theoretical upper bounds.",
            COLORS["yield"],
            COLORS["yield_light"],
            draw_yield_icon,
            r"$Y_{\mathrm{model}}>Y_{\max}^{\mathrm{theoretical}}$",
            0.190,
            0.58,
        ),
    ]
    card_h = 0.153
    y0 = y + h - 0.392 * h
    for i, (title, body, color, light, icon_fn, formula, icon_x_frac, icon_scale) in enumerate(cards):
        cy = y0 - i * 0.190
        rounded_box(ax, x + 0.050 * w, cy, 0.900 * w, card_h, light, color, lw=1.5, radius=0.014)
        icon_fn(ax, x + icon_x_frac * w, cy + card_h * 0.55, scale=icon_scale)
        ax.text(x + 0.335 * w, cy + card_h * 0.77, title, fontsize=8.7, fontweight="bold", color=color, va="top")
        ax.text(x + 0.335 * w, cy + card_h * 0.55, wrap_text(body, 37), fontsize=6.8, color=COLORS["ink"], va="top")
        ax.text(x + 0.335 * w, cy + card_h * 0.13, formula, fontsize=7.5, color=color, va="bottom")

    ax.text(
        x + 0.050 * w,
        y + 0.035 * h,
        "MQC targets these errors as mechanistic, correctable constraints rather than cosmetic annotations.",
        fontsize=7.0,
        color=COLORS["muted"],
        va="bottom",
    )


def draw_panel_b(ax):
    x, y, w, h = 0.382, 0.585, 0.583, 0.302
    add_panel(ax, x, y, w, h, "b", "Orchestrated MQC workflow", title_size=14.0)
    formula = r"$M^* = A_{\mathrm{Yield}} \circ A_{\mathrm{Bio}} \circ A_{\mathrm{Cycle}} \circ A_{\mathrm{Carbon}} \circ A_{\mathrm{NS}}(M)$"
    ax.text(x + 0.50 * w, y + h - 0.160 * h, formula, fontsize=15.5, ha="center", va="center", color="black")

    steps = [
        ("Input\nGEM\n$M$", "", "#FFFFFF", COLORS["ink"]),
        ("Standardization\n$A_{NS}$", "Namespace,\nformula, charge", COLORS["gray"], "#7B8794"),
        ("Carbon/medium\nrefinement\n$A_{Carbon}$", "Minimal carbon\nsource set", "#F1F7F1", "#78A06A"),
        ("Cycle correction\n$A_{Cycle}$", "Energy/redox/net\nproduction removal", COLORS["energy_light"], COLORS["energy"]),
        ("Biomass curation\n$A_{Bio}$", "Precursor gap-filling\nand 1 gDW normalization", COLORS["biomass_light"], COLORS["biomass"]),
        ("Yield validation\n$A_{Yield}$", "Theoretical yield\nupper bounds", COLORS["yield_light"], COLORS["yield"]),
        ("Corrected\nGEM\n$M^*$", "QC report", COLORS["green_light"], COLORS["green"]),
    ]

    step_y = y + 0.132 * h
    step_h = 0.415 * h
    gap = 0.010 * w
    widths = [0.070, 0.126, 0.132, 0.136, 0.143, 0.132, 0.098]
    total = sum(widths) * w + gap * (len(steps) - 1)
    start_x = x + (w - total) / 2
    xs = []
    cursor = start_x
    for idx, ((head, body, fc, ec), frac) in enumerate(zip(steps, widths)):
        bw = frac * w
        xs.append((cursor, bw))
        rounded_box(ax, cursor, step_y, bw, step_h, fc, ec, lw=1.6 if idx in [3, 4, 5, 6] else 1.2, radius=0.010)
        ax.text(cursor + bw / 2, step_y + 0.280 * step_h, head, fontsize=7.2, fontweight="bold", ha="center", va="bottom", color="black")
        if body:
            ax.text(cursor + bw / 2, step_y + 0.085 * step_h, body, fontsize=6.1, ha="center", va="bottom", color=COLORS["ink"])
        cursor += bw + gap

    for (sx, sw), (ex, ew) in zip(xs[:-1], xs[1:]):
        arrow(ax, (sx + sw + 0.002, step_y + step_h / 2), (ex - 0.003, step_y + step_h / 2), color=COLORS["dark_gray"], lw=1.4, ms=11)

    ax.text(
        x + 0.050 * w,
        y + 0.045 * h,
        "Narrative operator order; implementation maps these operators to initial standardization, NADH/ATP/XTP, net-production, yield, and biomass modules.",
        fontsize=6.6,
        color=COLORS["muted"],
        va="bottom",
    )


def penalty_funnel(ax, x, y, w, h):
    colors = ["#E27373", "#F4A36B", "#F9D884", "#BFD7A4", "#9CBAD8"]
    labels = [r"$20\delta_C$", r"$3\delta_M$", r"$3\delta_{Rule}$", r"$\delta_{Annot}$", "low"]
    for i in range(5):
        top_w = w * (1 - i * 0.14)
        bot_w = w * (1 - (i + 1) * 0.14)
        yy_top = y + h - i * h / 5
        yy_bot = y + h - (i + 1) * h / 5
        poly = Polygon(
            [
                (x + (w - top_w) / 2, yy_top),
                (x + (w + top_w) / 2, yy_top),
                (x + (w + bot_w) / 2, yy_bot),
                (x + (w - bot_w) / 2, yy_bot),
            ],
            closed=True,
            fc=colors[i],
            ec="white",
            lw=0.8,
            zorder=4,
        )
        ax.add_patch(poly)
        ax.text(x + w / 2, (yy_top + yy_bot) / 2, labels[i], fontsize=6.8, ha="center", va="center", color="black")
    ax.text(x - 0.010, y + h * 0.50, "Penalty", fontsize=7.2, rotation=90, va="center", ha="right", color=COLORS["ink"])


def draw_panel_c(ax):
    x, y, w, h = 0.382, 0.102, 0.583, 0.422
    add_panel(ax, x, y, w, h, "c", "Core algorithm: penalty-based rollback validation", title_size=13.0)

    ax.text(
        x + 0.035 * w,
        y + h - 0.130 * h,
        "Detect infeasible fluxes, rank candidate reactions by biochemical penalty, suppress them iteratively, and validate causality by rollback.",
        fontsize=7.5,
        color=COLORS["muted"],
        va="top",
    )

    top_y = y + h - 0.350 * h
    steps = [
        ("1. Detect", "Closed-input tests reveal infeasible flux."),
        ("2. Penalize", r"$P_r=20\delta_C+3\delta_M+3\delta_{Rule}+\delta_{Annot}$"),
        ("3. Suppress", "Constrain high-penalty reactions until flux reaches zero."),
        ("4. Rollback", "Reintroduce reactions one at a time."),
    ]
    box_w = 0.220 * w
    box_h = 0.168 * h
    for i, (head, body) in enumerate(steps):
        bx = x + 0.035 * w + i * 0.240 * w
        fc = "white" if i != 1 else "#FFF8E6"
        ec = COLORS["dark_gray"] if i != 1 else COLORS["orange"]
        rounded_box(ax, bx, top_y, box_w, box_h, fc, ec, lw=1.2, radius=0.010)
        ax.text(bx + 0.012 * w, top_y + box_h - 0.022, head, fontsize=8.0, fontweight="bold", color=COLORS["ink"], va="top")
        ax.text(bx + 0.012 * w, top_y + box_h - 0.052, wrap_text(body, 24), fontsize=5.8, color=COLORS["ink"], va="top")
        if i < 3:
            arrow(ax, (bx + box_w + 0.002, top_y + box_h / 2), (bx + 0.240 * w - 0.002, top_y + box_h / 2), lw=1.5, ms=10)

    # Left mechanism mini-panel: anomalous flux and penalty ranking.
    left_x = x + 0.035 * w
    mid_y = y + 0.205 * h
    mini_h = 0.232 * h
    rounded_box(ax, left_x, mid_y, 0.420 * w, mini_h, "#FFFFFF", COLORS["panel_edge"], lw=1.0, radius=0.010)
    ax.text(left_x + 0.018 * w, mid_y + mini_h - 0.032, "Penalty-guided suppression", fontsize=8.2, fontweight="bold", color=COLORS["ink"])
    nodes = [(left_x + 0.062 * w, mid_y + 0.095 * h), (left_x + 0.145 * w, mid_y + 0.148 * h), (left_x + 0.145 * w, mid_y + 0.050 * h), (left_x + 0.245 * w, mid_y + 0.095 * h)]
    for n in nodes:
        ax.add_patch(Circle(n, 0.011, fc="#D6E4D2", ec=COLORS["ink"], lw=1.0, zorder=4))
    arrow(ax, nodes[0], nodes[1], color=COLORS["energy"], lw=1.9, ms=11)
    arrow(ax, nodes[0], nodes[2], color=COLORS["dark_gray"], lw=1.3, ms=10)
    arrow(ax, nodes[1], nodes[3], color=COLORS["dark_gray"], lw=1.3, ms=10)
    arrow(ax, nodes[2], nodes[3], color=COLORS["dark_gray"], lw=1.3, ms=10)
    ax.text(nodes[0][0] - 0.008, nodes[0][1] + 0.052, "P=20", fontsize=7.0, color=COLORS["energy"])
    ax.text(nodes[1][0] + 0.010, nodes[1][1] + 0.014, "P=3", fontsize=6.8, color=COLORS["muted"])
    penalty_funnel(ax, left_x + 0.295 * w, mid_y + 0.034 * h, 0.090 * w, 0.150 * h)

    # Right mechanism mini-panel: rollback validation.
    right_x = x + 0.500 * w
    rounded_box(ax, right_x, mid_y, 0.455 * w, mini_h, "#FFFFFF", COLORS["panel_edge"], lw=1.0, radius=0.010)
    ax.text(right_x + 0.018 * w, mid_y + mini_h - 0.025, "Rollback validation", fontsize=8.0, fontweight="bold", color=COLORS["ink"])
    decision = (right_x + 0.225 * w, mid_y + 0.070 * h)
    diamond_half_w = 0.040
    diamond_half_h = 0.035
    diamond = Polygon(
        [
            (decision[0], decision[1] + diamond_half_h),
            (decision[0] + diamond_half_w, decision[1]),
            (decision[0], decision[1] - diamond_half_h),
            (decision[0] - diamond_half_w, decision[1]),
        ],
        closed=True,
        fc="#F8FAFC",
        ec=COLORS["ink"],
        lw=1.1,
        zorder=3,
    )
    ax.add_patch(diamond)
    ax.text(decision[0], decision[1], "Does infeasibility\nreturn?", fontsize=6.5, ha="center", va="center")
    safe_box = rounded_box(ax, right_x + 0.015 * w, mid_y + 0.032 * h, 0.130 * w, 0.062 * h, COLORS["green_light"], COLORS["green"], lw=1.2, radius=0.008)
    lock_box = rounded_box(ax, right_x + 0.315 * w, mid_y + 0.032 * h, 0.130 * w, 0.062 * h, COLORS["energy_light"], COLORS["energy"], lw=1.2, radius=0.008)
    ax.text(safe_box.get_x() + safe_box.get_width() / 2, safe_box.get_y() + safe_box.get_height() / 2, "Retained\nsafe", fontsize=6.4, ha="center", va="center", color=COLORS["green"])
    ax.text(lock_box.get_x() + lock_box.get_width() / 2, lock_box.get_y() + lock_box.get_height() / 2, "Constrained\n$v=0$", fontsize=6.4, ha="center", va="center", color=COLORS["energy"])
    arrow(ax, (decision[0] - diamond_half_w, decision[1]), (safe_box.get_x() + safe_box.get_width(), safe_box.get_y() + safe_box.get_height()), color=COLORS["green"], lw=1.5, ms=10)
    arrow(ax, (decision[0] + diamond_half_w, decision[1]), (lock_box.get_x(), lock_box.get_y() + lock_box.get_height()), color=COLORS["energy"], lw=1.5, ms=10)
    ax.text(decision[0] - 0.046, decision[1] + 0.018, "No", fontsize=6.8, color=COLORS["green"], ha="center")
    ax.text(decision[0] + 0.046, decision[1] + 0.018, "Yes", fontsize=6.8, color=COLORS["energy"], ha="center")

    footer = (
        r"$\delta_C$: carbon imbalance; $\delta_M$: mass/charge imbalance; "
        r"$\delta_{Rule}$: biochemical direction-rule violation; $\delta_{Annot}$: missing annotation."
    )
    ax.text(x + 0.035 * w, y + 0.055 * h, footer, fontsize=6.7, color=COLORS["muted"], va="bottom")


def draw_figure():
    fig = plt.figure(figsize=(15.0, 8.4), dpi=600)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(
        0.035,
        0.955,
        "A constraint-based MQC framework for correcting quantitative infeasibilities in GEMs",
        fontsize=17.5,
        fontweight="bold",
        color="black",
        va="top",
    )
    ax.text(
        0.035,
        0.918,
        "From raw genome-scale metabolic models to corrected, validated, and quality-controlled GEMs",
        fontsize=9.2,
        color=COLORS["muted"],
        va="top",
    )

    draw_panel_a(ax)
    draw_panel_b(ax)
    draw_panel_c(ax)
    return fig


def main():
    fig = draw_figure()
    for suffix in [".png", ".pdf", ".svg"]:
        kwargs = {"bbox_inches": "tight", "pad_inches": 0.05}
        if suffix == ".png":
            kwargs["dpi"] = 600
        fig.savefig(OUT_STEM.with_suffix(suffix), **kwargs)
    plt.close(fig)
    print(OUT_STEM.with_suffix(".png"))
    print(OUT_STEM.with_suffix(".pdf"))
    print(OUT_STEM.with_suffix(".svg"))


if __name__ == "__main__":
    main()
