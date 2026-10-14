#!/usr/bin/env python3
"""Generate 14 panels: ATP, NADH and biomass for manuscript Fig. 3."""

from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from panel_style import (
    DATA_DIR,
    DB_COLORS,
    DB_ORDER,
    GRID,
    INK,
    MUTED,
    apply_style,
    centered_title,
    publication_axes,
    save_panel,
)


def plot_profile_ecdf(ax, values, color, linestyle):
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    vals = np.sort(np.clip(vals, 1e-3, None))
    fractions = np.arange(1, len(vals) + 1) / len(vals)
    ax.plot(vals, fractions, color=color, linewidth=2.0, linestyle=linestyle)


def format_value(value: float) -> str:
    if abs(value) >= 100:
        return f"{value:.0f}"
    if abs(value) >= 10:
        return f"{value:.1f}"
    if abs(value) >= 1:
        return f"{value:.2f}"
    return f"{value:.2g}"


def draw_profiles(metric: str, manuscript_panel: str, source: str) -> None:
    profiles = pd.read_csv(DATA_DIR / "fig3" / source)
    profiles["value"] = pd.to_numeric(profiles["value"], errors="coerce")
    profiles = profiles.dropna(subset=["value"])
    for database in DB_ORDER:
        subset = profiles[profiles["database"].eq(database)]
        initial = subset[subset["stage"].eq("Initial")]["value"].to_numpy(dtype=float)
        final = subset[subset["stage"].eq("Final")]["value"].to_numpy(dtype=float)
        if len(initial) != len(final):
            raise ValueError(f"{metric}/{database}: unequal initial and final counts")
        color = DB_COLORS[database]
        fig, ax = plt.subplots(figsize=(4.4, 3.6))
        plot_profile_ecdf(ax, initial, color, "-")
        plot_profile_ecdf(ax, final, color, "--")
        ax.set_xscale("log")
        ax.set_xlim(1e-3, 1500)
        ax.set_ylim(0, 1.02)
        ax.set_xlabel("Recorded value (log scale)")
        ax.set_ylabel("Cumulative fraction")
        centered_title(ax, f"{database}: {metric} production profile")
        publication_axes(ax)
        ax.text(
            0.045,
            0.94,
            f"$n$ = {len(initial):,}\nInitial median  {format_value(float(np.median(initial)))}\n"
            f"After MQC median*  {format_value(float(np.median(final)))}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.0,
            bbox=dict(
                boxstyle="round,pad=0.28", facecolor="white", edgecolor=GRID, alpha=0.93
            ),
        )
        handles = [
            Line2D(
                [0],
                [0],
                color=color,
                lw=2.0,
                ls="-",
                label="Initial maximum-production flux",
            ),
            Line2D(
                [0],
                [0],
                color=color,
                lw=2.0,
                ls="--",
                label="Post-MQC carbon-normalized value*",
            ),
        ]
        ax.legend(
            handles=handles,
            frameon=True,
            facecolor="white",
            edgecolor=GRID,
            framealpha=0.96,
            fontsize=7.3,
            loc="lower right",
        )
        fig.tight_layout()
        save_panel(fig, "Fig3", f"Fig3{manuscript_panel}_{metric}_{database}")


def pct_in_range(series: pd.Series, lo: float, hi: float) -> float:
    values = pd.to_numeric(series, errors="coerce").dropna()
    return float(((values >= lo) & (values <= hi)).mean() * 100)


def violin_with_points(ax, values_by_db, y_max):
    positions = np.arange(1, len(DB_ORDER) + 1)
    arrays = []
    for db in DB_ORDER:
        values = np.asarray(values_by_db.get(db, []), dtype=float)
        values = values[np.isfinite(values)]
        arrays.append(np.clip(values, 0, y_max))
    violin = ax.violinplot(
        arrays,
        positions=positions,
        widths=0.78,
        showmeans=False,
        showmedians=False,
        showextrema=False,
    )
    for body, db in zip(violin["bodies"], DB_ORDER):
        body.set_facecolor(DB_COLORS[db])
        body.set_edgecolor(DB_COLORS[db])
        body.set_alpha(0.25)
        body.set_linewidth(0.8)
    rng = np.random.default_rng(7)
    for position, db, values in zip(positions, DB_ORDER, arrays):
        sample = (
            values
            if len(values) <= 700
            else rng.choice(values, size=700, replace=False)
        )
        jitter = rng.normal(0, 0.045, len(sample))
        ax.scatter(
            np.full(len(sample), position) + jitter,
            sample,
            s=5.0,
            color=DB_COLORS[db],
            alpha=0.42,
            linewidths=0,
            rasterized=True,
        )
        q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
        ax.add_patch(
            plt.Rectangle(
                (position - 0.14, q1),
                0.28,
                q3 - q1,
                edgecolor=INK,
                facecolor="white",
                linewidth=0.8,
                alpha=0.72,
            )
        )
        ax.plot(
            [position - 0.16, position + 0.16],
            [median, median],
            color=INK,
            linewidth=1.1,
        )
    ax.set_xticks(positions, DB_ORDER)


def draw_stage_distributions(growth: pd.DataFrame) -> None:
    settings = {
        "Initial": (230, "Initial-model biomass output"),
        "Middle": (2.0, "Intermediate biomass output"),
        "Final": (2.2, "Final normalized biomass output"),
    }
    for stage, (y_max, title) in settings.items():
        fig, ax = plt.subplots(figsize=(4.8, 3.8))
        values = {
            db: growth[(growth["stage"] == stage) & (growth["database"] == db)][
                "growth_rate_h1"
            ]
            .dropna()
            .to_numpy()
            for db in DB_ORDER
        }
        if stage == "Final":
            ax.axhspan(0, 1, color="#DDEFE2", alpha=0.65, zorder=0)
        violin_with_points(ax, values, y_max)
        ax.set_ylim(0, y_max)
        ax.set_ylabel(r"Biomass output (h$^{-1}$)")
        centered_title(ax, title)
        publication_axes(ax, grid_axis="y")
        if stage == "Middle":
            ax.text(
                0.98,
                0.95,
                r"Values >2 h$^{-1}$ clipped",
                transform=ax.transAxes,
                fontsize=7.5,
                color=MUTED,
                va="top",
                ha="right",
            )
        if stage == "Final":
            ax.text(
                0.98,
                0.96,
                "Shaded: 0–1 h$^{-1}$\nreference interval",
                transform=ax.transAxes,
                fontsize=7.5,
                color="#2E8B43",
                va="top",
                ha="right",
            )
        fig.tight_layout()
        save_panel(fig, "Fig3", f"Fig3c_{stage}_growth_distribution")


def plot_ecdf(ax, values, color, label, linestyle):
    vals = np.asarray(values, dtype=float)
    vals = np.sort(np.clip(vals[np.isfinite(vals)], 1e-3, None))
    fractions = np.arange(1, len(vals) + 1) / len(vals)
    ax.plot(
        vals, fractions, color=color, linewidth=2.0, linestyle=linestyle, label=label
    )


def draw_database_ecdf(database: str, data: pd.DataFrame) -> None:
    color = DB_COLORS[database]
    fig, ax = plt.subplots(figsize=(4.5, 3.6))
    ax.axvspan(1e-3, 1, color="#DDEFE2", alpha=0.55, zorder=0)
    plot_ecdf(ax, data["initial_growth_h1"].to_numpy(), color, "Initial", "-")
    plot_ecdf(ax, data["final_growth_h1"].to_numpy(), color, "Final", "--")
    ax.axvline(1, color="#2E8B43", linewidth=1.0, linestyle=":")
    ax.set_xscale("log")
    ax.set_xlim(1e-3, 250)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel(r"Growth rate (h$^{-1}$, log scale)")
    ax.set_ylabel("Cumulative fraction")
    centered_title(ax, f"{database} growth-rate distribution")
    publication_axes(ax)
    percentage = pct_in_range(data["final_growth_h1"], 0, 1)
    initial_median = pd.to_numeric(data["initial_growth_h1"], errors="coerce").median()
    final_median = pd.to_numeric(data["final_growth_h1"], errors="coerce").median()
    ax.text(
        0.03,
        0.18,
        f"Final in 0–1 h$^{{-1}}$: {percentage:.1f}%\nMedian: {initial_median:.1f} → {final_median:.2f}",
        transform=ax.transAxes,
        fontsize=8.0,
        bbox=dict(
            boxstyle="round,pad=0.25", facecolor="white", edgecolor=GRID, alpha=0.92
        ),
    )
    ax.legend(
        frameon=True,
        facecolor="white",
        edgecolor=GRID,
        framealpha=0.96,
        loc="lower right",
    )
    fig.tight_layout()
    save_panel(fig, "Fig3", f"Fig3d_{database}_growth_ECDF")


def draw_final_summary(growth: pd.DataFrame) -> None:
    final = growth[growth["stage"] == "Final"].copy()
    rows = []
    for db in DB_ORDER:
        values = final[final["database"] == db]["growth_rate_h1"]
        rows.append(
            (
                db,
                pct_in_range(values, 0, 1),
                float(pd.to_numeric(values, errors="coerce").median()),
            )
        )
    y = np.arange(len(rows))
    percentages = [row[1] for row in rows]
    fig, ax = plt.subplots(figsize=(4.7, 3.5))
    ax.axvspan(80, 100, color="#DDEFE2", alpha=0.55, zorder=0)
    bars = ax.barh(
        y,
        percentages,
        color=[DB_COLORS[row[0]] for row in rows],
        alpha=0.84,
        height=0.55,
    )
    for bar, (_, percentage, median) in zip(bars, rows):
        ax.text(
            min(percentage + 1.2, 99.2),
            bar.get_y() + bar.get_height() / 2,
            f"{percentage:.1f}%",
            va="center",
            ha="left",
        )
        ax.text(
            2.0,
            bar.get_y() + bar.get_height() / 2,
            f"median {median:.2f}",
            va="center",
            ha="left",
            fontsize=7.3,
            color="white",
            fontweight="bold",
        )
    ax.set_yticks(y, [row[0] for row in rows])
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel(r"Final models within 0–1 h$^{-1}$ (%)")
    centered_title(ax, "Final biomass-output range")
    publication_axes(ax, grid_axis="x")
    fig.tight_layout()
    save_panel(fig, "Fig3", "Fig3d_final_range_summary")


def main() -> None:
    apply_style()
    draw_profiles("ATP", "a", "atp_profiles.csv")
    draw_profiles("NADH", "b", "nadh_profiles.csv")
    growth = pd.read_csv(DATA_DIR / "fig3" / "growth_distributions.csv")
    growth["growth_rate_h1"] = pd.to_numeric(growth["growth_rate_h1"], errors="coerce")
    draw_stage_distributions(growth)
    for database, source in [
        ("CarveMe", "carveme_growth_ecdf.csv"),
        ("VMH", "vmh_growth_ecdf.csv"),
    ]:
        draw_database_ecdf(database, pd.read_csv(DATA_DIR / "fig3" / source))
    draw_final_summary(growth)
    print("Fig3: wrote 14 standalone panels")


if __name__ == "__main__":
    main()
