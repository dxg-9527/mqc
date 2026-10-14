"""Shared publication style and export helpers for all MQC panels."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt


PACKAGE = Path(__file__).resolve().parents[1]
DATA_DIR = PACKAGE / "source_data"
PANEL_DIR = PACKAGE / "panels"

FONT = "Liberation Sans"
INK = "#20262E"
MUTED = "#68727D"
GRID = "#D7DCE2"
DB_ORDER = ["CarveMe", "BiGG", "ModelSEED", "VMH"]
DB_COLORS = {
    "CarveMe": "#E6862E",
    "BiGG": "#2F73B7",
    "ModelSEED": "#54A65B",
    "VMH": "#D85863",
}


def apply_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": FONT,
            "font.size": 9.0,
            "axes.titlesize": 11.0,
            "axes.titleweight": "bold",
            "axes.labelsize": 9.0,
            "axes.labelcolor": INK,
            "axes.edgecolor": INK,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "text.color": INK,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "mathtext.fontset": "custom",
            "mathtext.rm": FONT,
            "mathtext.it": f"{FONT}:italic",
            "mathtext.bf": f"{FONT}:bold",
            "mathtext.cal": FONT,
            "mathtext.sf": FONT,
            "mathtext.tt": "Liberation Mono",
            "savefig.facecolor": "white",
        }
    )


def centered_title(ax: mpl.axes.Axes, title: str, *, pad: float = 10) -> None:
    ax.set_title(title, loc="center", fontsize=11.0, fontweight="bold", pad=pad)


def save_panel(
    fig: mpl.figure.Figure, figure: str, stem: str, *, tight: bool = True
) -> None:
    out_dir = PANEL_DIR / figure
    out_dir.mkdir(parents=True, exist_ok=True)
    kwargs = {"facecolor": "white"}
    if tight:
        kwargs.update({"bbox_inches": "tight", "pad_inches": 0.05})
    fig.savefig(out_dir / f"{stem}.png", dpi=600, **kwargs)
    fig.savefig(out_dir / f"{stem}.pdf", **kwargs)
    fig.savefig(out_dir / f"{stem}.svg", **kwargs)
    plt.close(fig)


def publication_axes(ax: mpl.axes.Axes, *, grid_axis: str = "both") -> None:
    ax.set_facecolor("#FBFCFD")
    ax.grid(axis=grid_axis, color=GRID, linewidth=0.6, linestyle="--", alpha=0.75)
    ax.set_axisbelow(True)
