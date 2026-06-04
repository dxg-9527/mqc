import ast
import math
from pathlib import Path

import pandas as pd


BASE_DIR = Path("/hpcfs/fhome/yuan_qq/jupyter/MQC结果分析/MQC 文章图/图 2 中的 8 个原始表")
OUTPUT_DIR = BASE_DIR / "analysis_outputs"
SOURCE_CSV = OUTPUT_DIR / "new_mechanistic_feature_matrix_source.csv"

NEW_FILES = {
    "BiGG": BASE_DIR / "bigg_new.xlsx",
    "VMH": BASE_DIR / "virtual_new.xlsx",
    "ModelSEED": BASE_DIR / "seed_new.xlsx",
    "CarveMe": BASE_DIR / "carveme_new.xlsx",
}

COFACTOR_COLS = ["NADH", "NADPH", "FADH2", "FMNH2", "Q8H2", "MQL8", "DMMQL8", "ATP", "CTP", "GTP", "UTP", "ITP"]


def parse_second_value(value):
    if pd.isna(value):
        return 0.0
    try:
        parsed = ast.literal_eval(str(value))
    except Exception:
        return 0.0
    if isinstance(parsed, tuple) and len(parsed) >= 2:
        try:
            second = float(parsed[1])
        except Exception:
            return 0.0
        if math.isnan(second) or math.isinf(second):
            return 0.0
        return second
    return 0.0


def yes_no_to_int(value):
    if pd.isna(value):
        return 0
    return 1 if str(value).strip().lower() in {"yes", "y", "true", "1"} else 0


def safe_float(value, default=0.0):
    try:
        out = float(value)
    except Exception:
        return default
    if math.isnan(out) or math.isinf(out):
        return default
    return out


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for database, path in NEW_FILES.items():
        df = pd.read_excel(path)
        for idx, row in df.iterrows():
            item = {
                "database": database,
                "row_id": f"{database}_{idx}",
                "model": row.get("model"),
            }
            for col in COFACTOR_COLS:
                item[col] = parse_second_value(row.get(col))
            item["metabolite_violation"] = yes_no_to_int(row.get("metabolite"))
            item["yield_violation"] = yes_no_to_int(row.get("yield"))
            item["biomass_gap"] = yes_no_to_int(row.get("is_gap"))
            item["gap_num_log1p"] = math.log1p(max(safe_float(row.get("gap_num")), 0.0))
            item["bio_norxn_fail"] = 1 - int(safe_float(row.get("bio_norxn"), 1))
            item["bio_nogrow_fail"] = 1 - int(safe_float(row.get("bio_nogrow"), 1))
            item["bio_coupling_fail"] = 1 - int(safe_float(row.get("bio_coupling"), 1))
            item["bio_1g_fail"] = 1 - int(safe_float(row.get("bio_1g"), 1))
            item["ini_weight_deviation"] = abs(safe_float(row.get("ini_weight(g)"), 1.0) - 1.0)
            item["carbon_biomass_ratio"] = safe_float(row.get("c_bio"), 0.0)
            rows.append(item)
    df = pd.DataFrame(rows)
    df.to_csv(SOURCE_CSV, index=False)
    print(SOURCE_CSV)


if __name__ == "__main__":
    main()
