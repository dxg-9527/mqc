"""Three condition-resolved Fig5 alternatives; archived numeric values are immutable."""

from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
import fitz

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "source_data"
DBS = ["CarveMe", "BiGG", "ModelSEED", "VMH"]
INK = "#273743"
MUTED = "#77828B"
GRID = "#E6EBEF"
UNIT = "h$^{-1}$"
FLUX = "mmol gDCW$^{-1}$ h$^{-1}$"
COLORS = {
    "glucose": "#2F73B7",
    "fructose": "#E6862E",
    "6pgc": "#8E67C6",
    "succ_glu": "#CB5262",
    "glycerol": "#139D9C",
    "pyruvate": "#B68A12",
    "malate": "#4D9B58",
    "glc_mal": "#805738",
    "acetate": "#BD5796",
    "succinate": "#6D68B0",
    "lactate": "#D65B4A",
    "2oxoglutarate": "#697886",
    "GlucoseUptake": "#2F73B7",
    "O2uptake": "#D85B61",
    "CO2production": "#E18A32",
    "NH4uptake": "#4C9B60",
}
MAP = {
    "glc": ("glucose", "Glucose"),
    "fru": ("fructose", "Fructose"),
    "6pgc": ("6pgc", "6pgc"),
    "succ&glu": ("succ_glu", "Succinate + Glu"),
    "glyc": ("glycerol", "Glycerol"),
    "pyr": ("pyruvate", "Pyruvate"),
    "mal": ("malate", "Malate"),
    "glc&mal": ("glc_mal", "Glucose + malate"),
    "Acetate": ("acetate", "Acetate"),
    "Succinate": ("succinate", "Succinate"),
    "L-Lactate": ("lactate", "L-Lactate"),
    "Pyruvate": ("pyruvate", "Pyruvate"),
    "2-Oxoglutarate": ("2oxoglutarate", "2-Oxoglutarate"),
    "D-Glucose": ("glucose", "D-Glucose"),
    "L-Malate": ("malate", "L-Malate"),
    "Glycerol": ("glycerol", "Glycerol"),
    "GlucoseUptake": ("GlucoseUptake", "Glucose uptake"),
    "O2uptake": ("O2uptake", "O$_2$ uptake"),
    "CO2production": ("CO2production", "CO$_2$ production"),
    "NH4uptake": ("NH4uptake", "NH$_4^+$ uptake"),
}
SPECS = [
    dict(
        key="bacillus",
        file="bacillus_substrate_growth.csv",
        title="$\\it{B.\\ subtilis}$ growth",
        cat="condition",
        obs="experimental_growth_h1",
        pre="initial_model_growth_h1",
        post="control_model_growth_h1",
        order=["glc", "fru", "6pgc", "succ&glu", "glyc", "pyr", "mal", "glc&mal"],
        unit=UNIT,
        lim=1.08,
        ticks=[0, 0.5, 1],
        reslim=(-0.85, 0.45),
        resticks=[-0.8, -0.4, 0, 0.4],
        heatlim=0.8,
    ),
    dict(
        key="flux",
        file="ecoli_exchange_fluxes.csv",
        title="$\\it{E.\\ coli}$ exchange fluxes",
        cat="measurement",
        obs="experimental_flux",
        pre="initial_model_flux",
        post="control_model_flux",
        order=["GlucoseUptake", "O2uptake", "CO2production", "NH4uptake"],
        unit=FLUX,
        lim=10.5,
        ticks=[0, 5, 10],
        reslim=(-9.5, 3.5),
        resticks=[-8, -4, 0],
        heatlim=9,
    ),
    dict(
        key="ecoli",
        file="ecoli_substrate_growth.csv",
        title="$\\it{E.\\ coli}$ growth",
        cat="substrate",
        obs="experimental_growth_h1",
        pre="initial_model_growth_h1",
        post="control_model_growth_h1",
        order=[
            "Acetate",
            "Succinate",
            "L-Lactate",
            "Pyruvate",
            "2-Oxoglutarate",
            "D-Glucose",
            "L-Malate",
            "Glycerol",
        ],
        unit=UNIT,
        lim=0.15,
        ticks=[0, 0.05, 0.1, 0.15],
        reslim=(-0.09, 0.065),
        resticks=[-0.08, -0.04, 0, 0.04],
        heatlim=0.09,
    ),
]
EXPORTS = {}
CMAP = LinearSegmentedColormap.from_list(
    "signed_error", ["#35669B", "#F9FAFB", "#C55A56"]
)
CMAP.set_bad("#ECEFF1")


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def style():
    plt.rcParams.update(
        {
            "font.family": "Liberation Sans",
            "font.size": 7,
            "text.color": INK,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "axes.edgecolor": MUTED,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.65,
            "axes.labelsize": 7,
            "xtick.labelsize": 6.7,
            "ytick.labelsize": 6.7,
            "mathtext.fontset": "dejavusans",
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def canvas(w, h):
    return plt.figure(figsize=(w / 25.4, h / 25.4))


def axmm(fig, x, y, w, h):
    fw, fh = fig.get_size_inches() * 25.4
    return fig.add_axes([x / fw, 1 - (y + h) / fh, w / fw, h / fh])


def text(fig, x, y, s, **kw):
    fw, fh = fig.get_size_inches() * 25.4
    return fig.text(x / fw, 1 - y / fh, s, va="center", **kw)


def color(cat):
    return COLORS[MAP[cat][0]]


def label(cat):
    return MAP[cat][1]


def tagged(d):
    return bool(d.initial_status.ne("exact_from_uploaded_workbook").any())


def keys(s):
    return (
        [(c, mu) for c in s["order"] for mu in [0.15, 0.2, 0.25, 0.3, 0.35, 0.4]]
        if s["key"] == "flux"
        else [(c, None) for c in s["order"]]
    )


def subset_row(d, s, c, mu):
    t = d[d[s["cat"]].eq(c)]
    if mu is not None:
        t = t[np.isclose(t.growth_rate_h1, mu)]
    assert len(t) <= 1
    return None if t.empty else t.iloc[0]


def export(fig, name, folder=ROOT):
    folder.mkdir(parents=True, exist_ok=True)
    (ROOT / "review").mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    outside = []
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_visible() and t.get_text().strip():
            b = t.get_window_extent(r)
            if (
                b.x0 < -2
                or b.y0 < -2
                or b.x1 > fig.bbox.width + 2
                or (b.y1 > fig.bbox.height + 2)
            ):
                outside.append(t.get_text())
    assert not outside, (name, outside)
    for ext in ["pdf", "svg", "png"]:
        fig.savefig(folder / f"{name}.{ext}", dpi=600)
    with fitz.open(folder / f"{name}.pdf") as doc:
        assert all((doc.extract_font(f[0])[3] for f in doc[0].get_fonts(full=True)))
        if folder == ROOT:
            doc[0].get_pixmap(dpi=200, alpha=False).save(
                ROOT / "review" / f"{name}.png"
            )
    EXPORTS[name] = {
        "width_mm": round(fig.get_figwidth() * 25.4, 3),
        "height_mm": round(fig.get_figheight() * 25.4, 3),
        "outside_text": outside,
    }
    plt.close(fig)
