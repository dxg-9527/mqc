#!/usr/bin/env python3
"""Reproduce the recommended Fig. 6 from the archived September inputs.

Only KZ1544, JCM6485 and JCM9990 enter the experimental figures and summaries.
Existing nine-substrate relative values are preserved; no GEM run or curve fit.
"""
from pathlib import Path
import hashlib
import json
import platform

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import fitz

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "source_data"
DATA = ROOT / "source_data"
STRAINS = ["KZ1544", "JCM6485", "JCM9990"]
SUBSTRATES = [
    "Mannose",
    "Sucrose",
    "D-Glucose",
    "D-Galactose",
    "D-Xylose",
    "Trehalose",
    "Cellobiose",
    "D-Fructose",
    "Maltose",
]
P = "model_pre_raw_growth_norm_by_strain"
Q = "model_post_raw_growth_norm_by_strain"
Y = "specific_growth_rate_norm_by_strain"
PRE, POST = "#4979AA", "#DB8C4E"
INK, MUTED, GRID = "#273743", "#73838D", "#E1E8EC"
CHECKS = {}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_data():
    provenance = {}
    names = [
        "publication_model_experiment_raw_values.csv",
        "global_carbon_utilization_summary.csv",
    ]
    for name in names:
        path = ARCHIVE / name
        provenance[name] = {"path": str(path), "sha256_before": sha(path)}
    full = pd.read_csv(ARCHIVE / names[0])
    source_cols = [
        "strain_abbrev",
        "strain",
        "model_id",
        "substrate_name",
        "substrate",
        "substrate_id",
        "n_replicates",
        "replicate_ids",
        "specific_growth_rate",
        "specific_growth_r2",
        "specific_growth_window_start_h",
        "specific_growth_window_end_h",
        "model_pre_raw_growth",
        "model_post_raw_growth",
        P,
        Q,
        Y,
        "experimental_mu_missing",
        "pre_model_value_missing",
        "post_model_value_missing",
        "slope_fit_r2_below_0_8",
    ]
    t = full.loc[full.strain_abbrev.isin(STRAINS), source_cols].copy()
    assert len(t) == 27 and not t.duplicated(["strain_abbrev", "substrate"]).any()
    assert set(t.strain_abbrev) == set(STRAINS)
    for strain in STRAINS:
        assert set(t.loc[t.strain_abbrev.eq(strain), "substrate"]) == set(SUBSTRATES)
    assert np.isfinite(t[[P, Q, Y]].to_numpy()).all()
    assert t[[P, Q, Y]].ge(0).all().all() and t[[P, Q, Y]].le(1).all().all()
    for raw, relative in [
        ("model_pre_raw_growth", P),
        ("model_post_raw_growth", Q),
        ("specific_growth_rate", Y),
    ]:
        maximum = t.groupby("strain_abbrev")[raw].transform("max")
        # An all-zero model profile remains all zero in the archived table.
        expected = t[raw] / maximum.mask(maximum.eq(0), 1)
        assert np.allclose(expected, t[relative], atol=1e-12, rtol=0)
    t["signed_pre"] = t[P] - t[Y]
    t["signed_post"] = t[Q] - t[Y]
    t["absolute_error_pre"] = t.signed_pre.abs()
    t["absolute_error_post"] = t.signed_post.abs()
    t["agreement_pre"] = 1 - t.absolute_error_pre
    t["agreement_post"] = 1 - t.absolute_error_post
    t["delta_agreement"] = t.agreement_post - t.agreement_pre
    t["error_direction"] = np.select(
        [t.delta_agreement > 1e-6, t.delta_agreement < -1e-6],
        ["lower", "higher"],
        default="unchanged",
    )
    assert np.allclose(t.agreement_pre, full.loc[t.index, "agreement_pre"])
    assert np.allclose(t.agreement_post, full.loc[t.index, "agreement_post"])
    t = (
        t.set_index(["strain_abbrev", "substrate"])
        .loc[pd.MultiIndex.from_product([STRAINS, SUBSTRATES])]
        .reset_index()
    )
    t.to_csv(DATA / "three_strains_nine_substrates.csv", index=False)
    numeric = [
        "agreement_pre",
        "agreement_post",
        "delta_agreement",
        "absolute_error_pre",
        "absolute_error_post",
    ]
    summary = t.groupby("strain_abbrev")[numeric].mean().reindex(STRAINS)
    summary["n_substrate_comparisons"] = 9
    summary.loc["Mean (3 strains)", numeric] = t[numeric].mean()
    summary.loc["Mean (3 strains)", "n_substrate_comparisons"] = 27
    summary.to_csv(DATA / "three_strain_summary.csv")
    counts = pd.crosstab(t.strain_abbrev, t.error_direction).reindex(
        index=STRAINS, columns=["lower", "unchanged", "higher"], fill_value=0
    )
    counts.to_csv(DATA / "three_strain_change_counts.csv")
    t.loc[
        t.slope_fit_r2_below_0_8, ["strain_abbrev", "substrate", "specific_growth_r2"]
    ].to_csv(DATA / "experimental_fit_flags.csv", index=False)
    global_data = pd.read_csv(ARCHIVE / names[1])
    assert len(global_data) == 52
    for stage in ["pre", "post"]:
        assert np.allclose(
            global_data[f"{stage}_active_count"] / 213 * 100,
            global_data[f"{stage}_active_percent"],
        )
    # Exactly the archived ranking rule: larger of the two utilization fractions.
    idx = (
        global_data[["pre_active_percent", "post_active_percent"]]
        .max(axis=1)
        .sort_values(ascending=False)
        .index[:12]
    )
    top = global_data.loc[idx].copy()
    global_data.to_csv(DATA / "MAG_utilization_52_substrates.csv", index=False)
    top.to_csv(DATA / "panel_a_top12_substrates.csv", index=False)
    CHECKS.update(
        {
            "strains": STRAINS,
            "substrates": SUBSTRATES,
            "experimental_pairs": len(t),
            "lower_error": int((t.error_direction == "lower").sum()),
            "unchanged_error": int((t.error_direction == "unchanged").sum()),
            "higher_error": int((t.error_direction == "higher").sum()),
            "normalization_verified": True,
            "fit_flags": int(t.slope_fit_r2_below_0_8.sum()),
            "source_provenance": provenance,
            "exports": {},
        }
    )
    return t, summary, top


def style():
    plt.rcParams.update(
        {
            "font.family": "Liberation Sans",
            "font.size": 7.5,
            "text.color": INK,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "axes.edgecolor": MUTED,
            "axes.labelsize": 7.5,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7.3,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.65,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def canvas(w, h):
    return plt.figure(figsize=(w / 25.4, h / 25.4))


def axmm(fig, x, top, w, h):
    fw, fh = fig.get_size_inches() * 25.4
    return fig.add_axes([x / fw, 1 - (top + h) / fh, w / fw, h / fh])


def txt(fig, x, y, text, **kwargs):
    fw, fh = fig.get_size_inches() * 25.4
    return fig.text(x / fw, 1 - y / fh, text, va="center", **kwargs)


def heading(fig, x, y, title, subtitle=None):
    txt(fig, x, y, title, ha="center", fontsize=8.5, fontweight="bold")
    if subtitle:
        txt(fig, x, y + 6, subtitle, ha="center", fontsize=7, color=MUTED)


def letter(fig, x, y, value):
    txt(fig, x, y, value, fontsize=11, fontweight="bold")


def legend(fig, y, experiment=False):
    handles = [
        Line2D([], [], ls="", marker="o", mfc="white", mec=PRE, ms=4, label="Pre-MQC"),
        Line2D([], [], ls="", marker="o", color=POST, ms=4, label="Post-MQC"),
    ]
    if experiment:
        handles.append(
            Line2D([], [], ls="", marker="D", color=INK, ms=3.5, label="Experiment")
        )
    fig.legend(
        handles=handles,
        loc="center",
        bbox_to_anchor=(0.55, 1 - y / (fig.get_figheight() * 25.4)),
        ncol=len(handles),
        frameon=False,
        fontsize=7.5,
        columnspacing=1.8,
    )


def clean(ax):
    ax.set_axisbelow(True)
    ax.grid(axis="x", color=GRID, lw=0.55)
    ax.spines["left"].set_visible(False)
    ax.tick_params(length=2.5, pad=3)
    ax.tick_params(axis="y", length=0)


def global_panel(ax, t):
    yy = np.arange(len(t))
    for y, before, after in zip(yy, t.pre_active_percent, t.post_active_percent):
        ax.plot([before, after], [y - 0.075, y + 0.075], color="#ACBDC7", lw=1.2)
    # Tiny offsets on the category axis preserve overlapping stage symbols.
    ax.scatter(
        t.pre_active_percent, yy - 0.075, s=20, fc="white", ec=PRE, lw=1, zorder=3
    )
    ax.scatter(t.post_active_percent, yy + 0.075, s=17, color=POST, zorder=4)
    ax.set_yticks(yy, t.substrate)
    ax.set_ylim(len(t) - 0.45, -0.55)
    ax.set_xlim(-2, 102)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("MAGs with predicted growth (%)")
    clean(ax)


def agreement_panel(ax, summary):
    yy = np.array([0, 1, 2, 3.35])
    ax.axhspan(2.8, 3.9, color="#F0F4F6", zorder=0)
    for y, (_, row) in zip(yy, summary.iterrows()):
        ax.plot(
            [row.agreement_pre, row.agreement_post],
            [y - 0.055, y + 0.055],
            color="#ACBDC7",
            lw=1.2,
        )
        ax.text(
            1.08,
            y,
            f"{row.delta_agreement:+.3f}",
            transform=ax.get_yaxis_transform(),
            fontsize=7.1,
            ha="center",
            va="center",
            clip_on=False,
        )
    ax.text(
        1.08,
        1.05,
        "Change",
        transform=ax.transAxes,
        ha="center",
        fontsize=7,
        color=MUTED,
    )
    ax.scatter(
        summary.agreement_pre, yy - 0.055, s=23, fc="white", ec=PRE, lw=1.05, zorder=3
    )
    ax.scatter(summary.agreement_post, yy + 0.055, s=19, color=POST, zorder=4)
    ax.set_yticks(yy, summary.index)
    ax.set_ylim(3.9, -0.6)
    ax.set_xlim(-0.02, 1.02)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1], ["0", ".25", ".5", ".75", "1"])
    ax.set_xlabel("Mean agreement (1 - MAE)")
    clean(ax)


def strain_panel(ax, table, strain, labels=True, profiles=False):
    t = table[table.strain_abbrev.eq(strain)].set_index("substrate").loc[SUBSTRATES]
    yy = np.arange(len(t))
    for y in yy[::2]:
        ax.axhspan(y - 0.45, y + 0.45, fc="#F5F7F9", zorder=0)
    if profiles:
        pre, post = t[P], t[Q]
        ax.scatter(t[Y], yy, marker="D", s=14, color=INK, zorder=6)
        ax.set_xlim(-0.06, 1.06)
        ax.set_xticks([0, 0.5, 1], ["0", ".5", "1"])
    else:
        pre, post = t.signed_pre, t.signed_post
        ax.axvline(0, color=INK, lw=0.85, ls=(0, (3, 2)), zorder=2)
        ax.set_xlim(-1.06, 1.06)
        ax.set_xticks([-1, -0.5, 0, 0.5, 1], ["-1", "-.5", "0", ".5", "1"])
    for y, before, after in zip(yy, pre, post):
        ax.plot(
            [before, after], [y - 0.10, y + 0.10], color="#A9B8C1", lw=1.0, zorder=2
        )
    ax.scatter(pre, yy - 0.10, s=22, fc="white", ec=PRE, lw=1.05, zorder=4)
    ax.scatter(post, yy + 0.10, s=18, color=POST, zorder=5)
    for y, flagged in zip(yy, t.slope_fit_r2_below_0_8):
        if flagged:
            ax.text(
                1.03,
                y,
                "*",
                transform=ax.get_yaxis_transform(),
                ha="left",
                va="center",
                fontsize=9,
                clip_on=False,
            )
    ax.set_yticks(yy, SUBSTRATES if labels else [""] * len(SUBSTRATES))
    ax.set_ylim(len(t) - 0.45, -0.55)
    clean(ax)


def export(fig, name, folder=ROOT):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    outside = []
    for obj in fig.findobj(matplotlib.text.Text):
        if obj.get_visible() and obj.get_text().strip():
            b = obj.get_window_extent(renderer)
            if (
                b.x0 < -2
                or b.y0 < -2
                or b.x1 > fig.bbox.width + 2
                or b.y1 > fig.bbox.height + 2
            ):
                outside.append(obj.get_text())
    assert not outside, (name, outside)
    for ext in ["pdf", "svg", "png"]:
        fig.savefig(folder / f"{name}.{ext}", dpi=600)
    with fitz.open(folder / f"{name}.pdf") as doc:
        pdf_text = "\n".join(page.get_text() for page in doc)
        assert all(
            excluded not in pdf_text for excluded in ["LGMJ132", "8025", "CP1", "All 6"]
        )
        doc[0].get_pixmap(dpi=200, alpha=False).save(ROOT / "review" / f"{name}.png")
        fonts = doc[0].get_fonts(full=True)
        embedded = all(bool(doc.extract_font(font[0])[3]) for font in fonts)
        assert embedded
    CHECKS["exports"][name] = {
        "width_mm": float(fig.get_figwidth() * 25.4),
        "height_mm": float(fig.get_figheight() * 25.4),
        "png_dpi": 600,
        "outside_text": outside,
        "pdf_fonts_embedded": embedded,
    }
    plt.close(fig)


def main():
    for folder in [DATA, ROOT / "panels", ROOT / "review"]:
        folder.mkdir(parents=True, exist_ok=True)
    style()
    table, summary, top = load_data()
    fig = canvas(183, 179)
    letter(fig, 3, 8, "a")
    heading(
        fig,
        58,
        8,
        "Predicted carbon-source utilization",
        "Prediction survey across 213 MAGs",
    )
    global_panel(axmm(fig, 29, 23, 62, 52), top)
    letter(fig, 99, 8, "b")
    heading(
        fig, 143, 8, "Agreement with experiment", "3 strains; 9 substrates per strain"
    )
    agreement_panel(axmm(fig, 122, 23, 45, 52), summary)
    legend(fig, 91)
    for i, strain in enumerate(STRAINS):
        x = 29 + 52 * i
        letter(fig, [3, 77, 129][i], 103, "cde"[i])
        heading(fig, x + 21.5, 103, strain)
        strain_panel(axmm(fig, x, 113, 43, 42), table, strain, labels=i == 0)
    txt(
        fig,
        102,
        166,
        "Predicted - experimental relative growth",
        ha="center",
        fontsize=8,
    )
    txt(
        fig,
        102,
        174,
        "Dashed line: agreement     * Experimental growth fit R² < 0.8",
        ha="center",
        fontsize=7,
        color=MUTED,
    )
    export(fig, "Fig6_recommended_three_strains")

    # Companion chosen from the original C paired-profile encoding.
    fig = canvas(183, 91)
    for i, strain in enumerate(STRAINS):
        x = 29 + 52 * i
        letter(fig, [3, 77, 129][i], 8, "abc"[i])
        heading(fig, x + 21.5, 8, strain)
        strain_panel(
            axmm(fig, x, 17, 43, 45), table, strain, labels=i == 0, profiles=True
        )
    txt(fig, 102, 72, "Within-strain relative growth", ha="center", fontsize=8)
    legend(fig, 81, experiment=True)
    txt(
        fig,
        102,
        88,
        "* Experimental growth fit R² < 0.8",
        ha="center",
        fontsize=7,
        color=MUTED,
    )
    export(fig, "FigS5_three_strain_profiles")

    fig = canvas(91, 94)
    heading(
        fig,
        55,
        8,
        "Predicted carbon-source utilization",
        "Prediction survey across 213 MAGs",
    )
    global_panel(axmm(fig, 28, 22, 59, 52), top)
    legend(fig, 90)
    export(fig, "Fig6a_MAG_utilization", ROOT / "panels")
    fig = canvas(94, 94)
    heading(
        fig, 52, 8, "Agreement with experiment", "3 strains; 9 substrates per strain"
    )
    agreement_panel(axmm(fig, 29, 23, 50, 52), summary)
    legend(fig, 90)
    export(fig, "Fig6b_three_strain_agreement", ROOT / "panels")
    for panel, strain in zip("cde", STRAINS):
        fig = canvas(88, 83)
        heading(fig, 56, 8, strain)
        strain_panel(axmm(fig, 28, 17, 53, 44), table, strain)
        txt(
            fig,
            53,
            71,
            "Predicted - experimental relative growth",
            ha="center",
            fontsize=7.1,
        )
        legend(fig, 79)
        export(fig, f"Fig6{panel}_{strain}_deviation", ROOT / "panels")

    for item in CHECKS["source_provenance"].values():
        item["sha256_after"] = sha(Path(item["path"]))
        assert item["sha256_after"] == item["sha256_before"]
    CHECKS["environment"] = {
        "python": platform.python_version(),
        "matplotlib": matplotlib.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "PyMuPDF": fitz.VersionBind,
    }
    (ROOT / "review" / "validation.json").write_text(
        json.dumps(CHECKS, indent=2, ensure_ascii=False) + "\n"
    )
    print(
        json.dumps(
            {
                "output": str(ROOT),
                "pairs": len(table),
                "mean": summary.loc["Mean (3 strains)"].to_dict(),
                "exports": len(CHECKS["exports"]),
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
