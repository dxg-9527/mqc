#!/usr/bin/env python3
"""Author-approved Fig5 waterfall, rebuilt from packaged archived input tables."""
from matplotlib.patches import Patch
from plot_helpers import *

TOL = 1e-10
SHORT = {
    "glc": "Glc",
    "fru": "Fru",
    "6pgc": "6pgc",
    "succ&glu": "Succ+Glu",
    "glyc": "Glyc",
    "pyr": "Pyr",
    "mal": "Mal",
    "glc&mal": "Glc+Mal",
    "Acetate": "Acetate",
    "Succinate": "Succ",
    "L-Lactate": "L-Lac",
    "Pyruvate": "Pyr",
    "2-Oxoglutarate": "2-OG",
    "D-Glucose": "Glc",
    "L-Malate": "L-Mal",
    "Glycerol": "Glyc",
    "GlucoseUptake": "Glc",
    "O2uptake": "O$_2$",
    "CO2production": "CO$_2$",
    "NH4uptake": "NH$_4^+$",
}
LIMITS = {"bacillus": (-0.08, 0.78), "flux": (-2.4, 6.0), "ecoli": (-0.018, 0.075)}
TICKS = {
    "bacillus": [0, 0.2, 0.4, 0.6],
    "flux": [-2, 0, 2, 4, 6],
    "ecoli": [0, 0.02, 0.04, 0.06],
}
MANIFEST = []


def short(c, mu=None):
    return SHORT[c] + (f" {mu:.2f}" if mu is not None else "")


def summary(s, d):
    pre = np.abs(d[s["pre"]] - d[s["obs"]])
    post = np.abs(d[s["post"]] - d[s["obs"]])
    return float(pre.mean()), float(post.mean())


def tidy(s, d):
    records = []
    for db in DBS:
        block = d[d.database.eq(db)]
        for i, (c, mu) in enumerate(keys(s)):
            r = subset_row(block, s, c, mu)
            if r is None:
                vals = [np.nan] * 3
                status = "missing"
                gain = np.nan
            else:
                vals = [float(r[s[z]]) for z in ["obs", "pre", "post"]]
                status = r.initial_status
                gain = abs(vals[1] - vals[0]) - abs(vals[2] - vals[0])
            records.append(
                dict(
                    dataset=s["key"],
                    database=db,
                    condition=c,
                    growth_rate_h1=np.nan if mu is None else mu,
                    order=i,
                    present=r is not None,
                    observed=vals[0],
                    initial=vals[1],
                    post_mqc=vals[2],
                    initial_status=status,
                    error_reduction=gain,
                    direction=(
                        "missing"
                        if r is None
                        else (
                            "improved"
                            if gain > TOL
                            else "worsened" if gain < -TOL else "unchanged"
                        )
                    ),
                )
            )
    return pd.DataFrame(records)


def waterfall(ax, s, records, labels=True):
    t = (
        records[records.present]
        .sort_values(
            ["error_reduction", "order"], ascending=[False, True], kind="stable"
        )
        .copy()
    )
    t["rank"] = np.arange(1, len(t) + 1)
    for k, r in enumerate(t.itertuples()):
        hatch = "////" if r.direction == "worsened" else None
        ax.barh(
            k,
            r.error_reduction,
            height=0.65,
            color=color(r.condition),
            edgecolor=INK if hatch else "none",
            linewidth=0.25,
            hatch=hatch,
            zorder=3,
        )
        if abs(r.error_reduction) <= TOL:
            ax.scatter([0], [k], s=12, color=color(r.condition), marker="|", zorder=4)
    lo, hi = LIMITS[s["key"]]
    assert t.error_reduction.min() >= lo and t.error_reduction.max() <= hi
    ax.set(
        xlim=(lo, hi),
        ylim=(len(t) - 0.5, -0.5),
        yticks=np.arange(len(t)),
        xticks=TICKS[s["key"]],
    )
    ax.set_yticklabels(
        [
            short(r.condition, None if pd.isna(r.growth_rate_h1) else r.growth_rate_h1)
            for r in t.itertuples()
        ],
        fontsize=5.8,
    )
    for tick, r in zip(ax.get_yticklabels(), t.itertuples()):
        tick.set_color(color(r.condition))
    ax.axvline(0, color=INK, lw=0.75, zorder=2)
    ax.grid(axis="x", color=GRID, lw=0.45)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=2)
    ax.tick_params(axis="x", length=2, pad=2, labelsize=5.8)
    return t


def save(fig, name, kind, s=None, db=None):
    folder = ROOT if s is None else ROOT / "panels"
    export(fig, name, folder)
    if s:
        row = dict(
            panel=name,
            kind=kind,
            dataset=s["key"],
            database=db,
            source="source_data/" + s["file"],
            code="code/draw_fig5_waterfall.py",
        )
        for ext in ["png", "pdf", "svg"]:
            p = folder / (name + "." + ext)
            row[ext] = str(p.relative_to(ROOT))
            row[ext + "_sha256"] = digest(p)
        MANIFEST.append(row)


def main():
    style()
    tables = {s["key"]: pd.read_csv(DATA / s["file"]) for s in SPECS}
    records = {s["key"]: tidy(s, tables[s["key"]]) for s in SPECS}
    long = pd.concat(records.values(), ignore_index=True)
    assert (
        len(long) == 160
        and long.present.sum() == 155
        and long[["observed", "initial", "post_mqc"]].notna().sum().sum() == 465
    )
    long.to_csv(DATA / "condition_values_and_changes.csv", index=False)
    mae = []
    ranks = []
    fig = canvas(183, 266)
    for i, s in enumerate(SPECS):
        top, h = [(30, 30), (91, 83), (205, 30)][i]
        text(fig, 3, top - 22, "abc"[i], fontsize=11, fontweight="bold")
        text(
            fig,
            17,
            top - 22,
            s["title"] + " (" + s["unit"] + ")",
            fontsize=8.6,
            fontweight="bold",
        )
        for j, db in enumerate(DBS):
            d = tables[s["key"]].query("database == @db")
            r = records[s["key"]].query("database == @db")
            ax = axmm(fig, 17 + 44 * j, top, 28, h)
            rank = waterfall(ax, s, r)
            ranks.append(rank)
            text(
                fig,
                31 + 44 * j,
                top - 11,
                db + ("†" if tagged(d) else ""),
                ha="center",
                fontsize=7.8,
                fontweight="bold",
            )
            pos = int((r.direction == "improved").sum())
            neg = int((r.direction == "worsened").sum())
            eq = int((r.direction == "unchanged").sum())
            text(
                fig,
                31 + 44 * j,
                top - 5,
                f"{pos} better / {neg} worse / {eq} equal",
                ha="center",
                fontsize=5.1,
                color=MUTED,
            )
            a, b = summary(s, d)
            mae.append(
                dict(
                    dataset=s["key"],
                    database=db,
                    n=len(d),
                    mae_initial=a,
                    mae_post_mqc=b,
                    mae_reduction=a - b,
                )
            )
        if i == 2:
            text(
                fig,
                97,
                top + h + 10,
                "Absolute-error reduction:  better →",
                ha="center",
                fontsize=6.6,
            )
    text(
        fig,
        94,
        252,
        "Reduction = |Initial − Experiment| − |Post-MQC − Experiment|. Bars sorted within each model.",
        ha="center",
        fontsize=5.9,
    )
    text(
        fig,
        94,
        259,
        "Hatched: worsened. Zero mark: unchanged. Flux labels include growth rate (h$^{-1}$). † See source provenance.",
        ha="center",
        fontsize=5.8,
        color=MUTED,
    )
    save(fig, "Fig5_final_waterfall", "waterfall")
    pd.concat(ranks, ignore_index=True).to_csv(
        DATA / "waterfall_ranked_coordinates.csv", index=False
    )
    pd.DataFrame(mae).to_csv(DATA / "mae_summary.csv", index=False)
    for s in SPECS:
        for db in DBS:
            h = 108 if s["key"] == "flux" else 77
            fig = canvas(88, h)
            ax = axmm(fig, 23, 17, 59, h - 32)
            waterfall(ax, s, records[s["key"]].query("database == @db"))
            ax.set_xlabel("Absolute-error reduction (" + s["unit"] + ")", fontsize=6.6)
            text(fig, 50, 7, db, ha="center", fontsize=8)
            save(
                fig,
                "Fig5"
                + {"bacillus": "a", "flux": "b", "ecoli": "c"}[s["key"]]
                + "_"
                + db,
                "waterfall",
                s,
                db,
            )

    pd.DataFrame(MANIFEST).to_csv(ROOT / "panel_manifest.csv", index=False)
    report = dict(
        input_records=155,
        display_positions=160,
        missing_positions=5,
        source_raw_values=465,
        displayed_changes=155,
        source_derived_initial_records=int(
            (
                long.present & long.initial_status.ne("exact_from_uploaded_workbook")
            ).sum()
        ),
        input_hashes={s["file"]: digest(DATA / s["file"]) for s in SPECS},
        exports=EXPORTS,
    )
    (ROOT / "review/validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(dict(exports=len(EXPORTS), panels=len(MANIFEST)), indent=2))


if __name__ == "__main__":
    main()
