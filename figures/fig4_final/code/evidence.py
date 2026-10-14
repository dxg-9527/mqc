"""Archived values used in native PowerPoint panels; never run optimization."""

import csv
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "source_data" / "fig4"


def rows(name):
    with (DATA / name).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


DIAG = {
    r["target"]: float(r["initial"]) for r in rows("diagnostic_quantitative_data.csv")
}
ROUTES = {
    (r["path"], r["reaction_id"]): float(r["flux"])
    for r in rows("causal_path_fluxes.csv")
}


def initial(target):
    return f"{DIAG[target]:.3f}"


def flux(case, reaction, signed=False):
    value = ROUTES[(case, reaction)]
    return f"{value:+.3f}" if signed else f"{value:.3f}"


assert len({DIAG[t] for t in ["ATP", "CTP", "GTP", "UTP", "ITP"]}) == 1
