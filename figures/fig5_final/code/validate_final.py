#!/usr/bin/env python3
"""Verify raw coordinates, ranking, aggregate errors and submission exports."""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
import fitz
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "source_data"
CONFIG = {
    "bacillus": ("bacillus_substrate_growth.csv", "condition", "growth_h1"),
    "flux": ("ecoli_exchange_fluxes.csv", "measurement", "flux"),
    "ecoli": ("ecoli_substrate_growth.csv", "substrate", "growth_h1"),
}
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    report = json.loads((ROOT / "review/validation.json").read_text())
    coords = pd.read_csv(DATA / "condition_values_and_changes.csv")
    ranks = pd.read_csv(DATA / "waterfall_ranked_coordinates.csv")
    mae = pd.read_csv(DATA / "mae_summary.csv")
    n = 0
    maxerr = 0.0
    source_counts = {}
    originals = 0
    for key, (filename, category, suffix) in CONFIG.items():
        table = pd.read_csv(DATA / filename)
        source_counts[key] = len(table)
        assert sha(DATA / filename) == report["input_hashes"][filename]
        original = ROOT.parent / "source_data/fig5" / filename
        if original.exists():
            assert sha(original) == sha(DATA / filename)
            originals += 1
        fields = [
            "experimental_" + suffix,
            "initial_model_" + suffix,
            "control_model_" + suffix,
        ]
        for row in coords[coords.dataset.eq(key)].itertuples():
            t = table[
                (table.database == row.database) & (table[category] == row.condition)
            ]
            if key == "flux":
                t = t[np.isclose(t.growth_rate_h1, row.growth_rate_h1)]
            if not row.present:
                assert (
                    t.empty
                    and np.isnan(
                        [row.observed, row.initial, row.post_mqc, row.error_reduction]
                    ).all()
                )
                continue
            assert len(t) == 1
            raw = t.iloc[0][fields].to_numpy(dtype=float)
            shown = np.array([row.observed, row.initial, row.post_mqc])
            maxerr = max(maxerr, float(np.max(abs(raw - shown))))
            assert np.allclose(raw, shown, rtol=1e-13, atol=1e-14)
            assert row.initial_status == t.iloc[0].initial_status
            delta = abs(raw[1] - raw[0]) - abs(raw[2] - raw[0])
            assert np.isclose(delta, row.error_reduction, atol=1e-14)
            direction = (
                "improved"
                if delta > 1e-10
                else "worsened" if delta < -1e-10 else "unchanged"
            )
            assert row.direction == direction
            rr = ranks[
                (ranks.dataset == key)
                & (ranks.database == row.database)
                & (ranks.condition == row.condition)
            ]
            if key == "flux":
                rr = rr[np.isclose(rr.growth_rate_h1, row.growth_rate_h1)]
            assert len(rr) == 1 and np.isclose(rr.iloc[0].error_reduction, delta)
            n += 1
        for db, t in table.groupby("database"):
            r = mae[(mae.dataset == key) & (mae.database == db)].iloc[0]
            a = abs(t[fields[1]] - t[fields[0]]).mean()
            b = abs(t[fields[2]] - t[fields[0]]).mean()
            assert np.allclose(
                [r.mae_initial, r.mae_post_mqc, r.mae_reduction], [a, b, a - b]
            )
            ranked = ranks[(ranks.dataset == key) & (ranks.database == db)].sort_values(
                "rank"
            )
            assert len(ranked) == len(t) and np.all(
                np.diff(ranked.error_reduction) <= 1e-14
            )
            assert np.allclose(ranked.error_reduction.mean(), a - b)
    assert (
        n == 155
        and len(ranks) == 155
        and len(coords) == 160
        and coords.present.sum() == 155
    )
    assert maxerr == 0
    assert len(mae) == 12 and int((mae.mae_reduction > 0).sum()) == 10
    manifest = pd.read_csv(ROOT / "panel_manifest.csv")
    assert len(manifest) == 12
    for row in manifest.itertuples():
        for ext in ["png", "pdf", "svg"]:
            assert sha(ROOT / getattr(row, ext)) == getattr(row, ext + "_sha256")
        assert (ROOT / row.source).is_file() and (ROOT / row.code).is_file()
    for name, size in report["exports"].items():
        folder = ROOT if name.startswith("Fig5_") else ROOT / "panels"
        with Image.open(folder / (name + ".png")) as im:
            assert all(abs(d - 600) < 0.1 for d in im.info["dpi"])
        with fitz.open(folder / (name + ".pdf")) as doc:
            assert len(doc) == 1 and all(
                doc.extract_font(f[0])[3] for f in doc[0].get_fonts(full=True)
            )
            assert abs(doc[0].rect.width / 72 * 25.4 - size["width_mm"]) < 0.02
        assert size["outside_text"] == []
    report["independent_validation"] = dict(
        passed=True,
        raw_values_checked=465,
        max_raw_difference=maxerr,
        ranked_changes_checked=155,
        mae_pairs_checked=12,
        source_counts=source_counts,
        original_input_tables_checked=originals,
        direction_counts=coords.direction.value_counts().to_dict(),
        panel_files_checked=36,
        exported_figures_checked=13,
    )
    (ROOT / "review/validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["independent_validation"], indent=2))


if __name__ == "__main__":
    main()
