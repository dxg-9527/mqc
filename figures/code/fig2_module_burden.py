"""Count affected initial-QC modules per model, retaining the full core cohort."""

from pathlib import Path
import csv
import hashlib

import numpy as np
import pandas as pd

PACKAGE = Path(__file__).resolve().parents[1]
DATA = PACKAGE / "source_data/fig2"
COUNTS = {"CarveMe": 5585, "BiGG": 73, "ModelSEED": 79, "VMH": 806}
FLAGS = {
    "Reducing equivalents": "reducing_power_pass_initial",
    "Energy": "energy_pass_initial",
    "Net metabolite production": "metabolite_pass_initial",
    "Yield": "yield_pass_initial",
    "Biomass": "biomass_pass_initial",
}


def derive() -> pd.DataFrame:
    """Reconcile row-level flags with both existing prevalence summaries."""
    core = pd.read_csv(DATA / "core_model_qc_flags.csv")
    assert len(core) == 6543 and core.stable_model_id.is_unique
    assert core.groupby("database").size().to_dict() == COUNTS
    assert core[list(FLAGS.values())].isin([0, 1]).all().all()
    burden = core[list(FLAGS.values())].eq(0).sum(axis=1)
    assert burden.gt(0).equals(core.all_modules_pass_initial.eq(0))
    module = pd.read_csv(DATA / "module_error_prevalence.csv").set_index(
        ["database", "module"]
    )
    overall = pd.read_csv(DATA / "overall_error_prevalence.csv").set_index("database")
    rows = []
    for db, total in COUNTS.items():
        selected = core.database.eq(db)
        values = burden[selected]
        assert int(values.gt(0).sum()) == int(overall.loc[db, "error_models"])
        for label, flag in FLAGS.items():
            assert int(core.loc[selected, flag].eq(0).sum()) == int(
                module.loc[(db, label), "error_models"]
            )
        for count in range(6):
            affected = int(values.eq(count).sum())
            rows.append(
                dict(
                    database=db,
                    affected_modules=count,
                    models=affected,
                    total_models=total,
                    model_pct=100 * affected / total,
                    mean_affected_modules=float(values.mean()),
                    source_file="source_data/fig2/core_model_qc_flags.csv",
                    definition="Number of failed initial-QC modules (0=pass all; 5=fail all)",
                )
            )
    result = pd.DataFrame(rows)
    np.testing.assert_allclose(result.groupby("database").model_pct.sum(), 100)
    return result


def refresh() -> None:
    result = derive()
    result.to_csv(
        DATA / "module_burden_distribution.csv", index=False, encoding="utf-8-sig"
    )
    manifest = PACKAGE / "source_data_manifest.csv"
    rows = []
    if manifest.exists():
        with manifest.open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
    names = ["core_model_qc_flags.csv", "module_burden_distribution.csv"]
    rows = [
        r
        for r in rows
        if r["package_file"] not in [f"source_data/fig2/{n}" for n in names]
    ]
    for name, n, upstream, element, purpose in [
        (
            names[0],
            6543,
            "core_model_roster_6543.csv",
            "whole_file",
            "Fig2 model-level initial-QC flags; 0=fail, 1=pass",
        ),
        (
            names[1],
            len(result),
            "source_data/fig2/core_model_qc_flags.csv",
            "code/fig2_module_burden.py:derive",
            "Fig2d distribution of the number of affected QC modules",
        ),
    ]:
        path = DATA / name
        # Upstream roster identity is recorded even in a relocated, self-contained package.
        rows.append(
            dict(
                package_file=str(path.relative_to(PACKAGE)),
                upstream_source=upstream,
                upstream_element=element,
                rows=n,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                purpose=purpose,
            )
        )
    rows.sort(key=lambda r: r["package_file"])
    with manifest.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    refresh()
