from __future__ import annotations

import ast
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook


ROOT = Path("/home/dengxg/project/mqc/others/定量分析")
TMP_DIR = Path("/home/dengxg/project/mqc/tmp")
PAPER_DIR = ROOT / "result" / "论文"
FIG5_DIR = ROOT / "result" / "定量数据"
ARTICLE_WORKBOOK = ROOT / "定量生长数据文章 (version 1).xlsx"

OUTPUT_XLSX = PAPER_DIR / "MQC图2-5原始数据整理.xlsx"
OUTPUT_NOTE = PAPER_DIR / "MQC图2-5原始数据说明.md"
SOURCE_SCRIPT = ROOT / "figure_audit_workspace" / "remediation_scripts" / "定量分析.py"

DATABASE_FILES = {
    "CarveMe": TMP_DIR / "carveme_new.xlsx",
    "BiGG": TMP_DIR / "bigg_new.xlsx",
    "ModelSEED": TMP_DIR / "seed_new.xlsx",
    "VMH": TMP_DIR / "virtual_new.xlsx",
}


def parse_tuple(value: object, min_len: int) -> list[float] | None:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text == "无" or "nan" in text.lower():
        return None
    try:
        parsed = ast.literal_eval(text)
    except (SyntaxError, ValueError):
        return None
    if not isinstance(parsed, (list, tuple)) or len(parsed) < min_len:
        return None
    values: list[float] = []
    for item in parsed[:min_len]:
        try:
            values.append(float(item))
        except (TypeError, ValueError):
            return None
    return values


def article_raw(sheet_name: str) -> pd.DataFrame:
    return pd.read_excel(ARTICLE_WORKBOOK, sheet_name=sheet_name, header=None)


def as_float(value: object) -> float | None:
    if pd.isna(value):
        return None
    return float(value)


def build_overview_sheet() -> pd.DataFrame:
    rows = [
        {
            "figure_panel": "MQC图2-a",
            "data_status": "exact",
            "source_type": "excel",
            "source_path": str(TMP_DIR / "web_CARVEME_COMMEN_check.xlsx")
            + " | "
            + str(TMP_DIR / "bigg_analysis_check.xlsx")
            + " | "
            + str(TMP_DIR / "seed_analysis_check.xlsx")
            + " | "
            + str(TMP_DIR / "virtual_analysis_check.xlsx"),
            "note": "数据库总体错误率，可精确追溯。",
        },
        {
            "figure_panel": "MQC图2-b",
            "data_status": "not_applicable",
            "source_type": "conceptual",
            "source_path": "",
            "note": "示意图，无独立数值原始表。",
        },
        {
            "figure_panel": "MQC图2-c",
            "data_status": "exact",
            "source_type": "excel",
            "source_path": str(TMP_DIR / "carveme_check.xlsx")
            + " | "
            + str(TMP_DIR / "bigg_check.xlsx")
            + " | "
            + str(TMP_DIR / "seed_check.xlsx")
            + " | "
            + str(TMP_DIR / "virtual_check.xlsx"),
            "note": "QC 模块错误分布，可精确追溯。",
        },
        {
            "figure_panel": "MQC图2-d",
            "data_status": "exact",
            "source_type": "excel",
            "source_path": " | ".join(str(path) for path in DATABASE_FILES.values()),
            "note": "biomass 错误类型分布，可精确追溯。",
        },
        {
            "figure_panel": "MQC图3-a",
            "data_status": "not_applicable",
            "source_type": "conceptual",
            "source_path": "",
            "note": "目标函数/流程说明图，无独立数值原始表。",
        },
        {
            "figure_panel": "MQC图3-b",
            "data_status": "exact",
            "source_type": "excel",
            "source_path": " | ".join(str(path) for path in DATABASE_FILES.values()),
            "note": "ATP/GTP/CTP 箱线图原始值，可精确追溯。",
        },
        {
            "figure_panel": "MQC图3-c",
            "data_status": "exact",
            "source_type": "excel",
            "source_path": " | ".join(str(path) for path in DATABASE_FILES.values()),
            "note": "NADH/NADPH/Q8H2 箱线图原始值，可精确追溯。",
        },
        {
            "figure_panel": "MQC图3-d",
            "data_status": "not_applicable",
            "source_type": "conceptual",
            "source_path": "",
            "note": "流程示意图，无独立数值原始表。",
        },
        {
            "figure_panel": "MQC图4-a",
            "data_status": "exact",
            "source_type": "excel",
            "source_path": " | ".join(str(path) for path in DATABASE_FILES.values()),
            "note": "Initial/Middle/Final 生长率分布，可精确追溯。",
        },
        {
            "figure_panel": "MQC图4-b",
            "data_status": "exact",
            "source_type": "excel",
            "source_path": " | ".join(str(path) for path in DATABASE_FILES.values()),
            "note": "CarveMe 与 VMH 逐模型折线图数据，可精确追溯。",
        },
        {
            "figure_panel": "MQC图5-a",
            "data_status": "mostly_exact",
            "source_type": "uploaded_excel + legacy_script",
            "source_path": str(ARTICLE_WORKBOOK) + " | " + str(SOURCE_SCRIPT),
            "note": "上传工作簿补齐了 Bacillus 的 experimental、CarveMe final、ModelSEED final、BiGG initial/final、VMH final；缺失的初始列按单独 status/source 标注。",
        },
        {
            "figure_panel": "MQC图5-b",
            "data_status": "exact",
            "source_type": "uploaded_excel",
            "source_path": str(ARTICLE_WORKBOOK),
            "note": "E. coli 吸收/外排散点图原始值已切换为上传工作簿中的精确表。",
        },
        {
            "figure_panel": "MQC图5-c",
            "data_status": "exact",
            "source_type": "uploaded_excel + legacy_script",
            "source_path": str(ARTICLE_WORKBOOK) + " | " + str(SOURCE_SCRIPT),
            "note": "E. coli 生长率散点图已切换为上传工作簿中的精确表；ModelSEED 初始全 0 由旧脚本补足。",
        },
    ]
    return pd.DataFrame(rows)


def build_fig2a() -> pd.DataFrame:
    files = {
        "CarveMe": TMP_DIR / "web_CARVEME_COMMEN_check.xlsx",
        "BiGG": TMP_DIR / "bigg_analysis_check.xlsx",
        "ModelSEED": TMP_DIR / "seed_analysis_check.xlsx",
        "VMH": TMP_DIR / "virtual_analysis_check.xlsx",
    }
    rows: list[dict[str, object]] = []
    for database, path in files.items():
        df = pd.read_excel(path)
        error_count = int(df["all"].eq(0).sum())
        total_models = int(len(df))
        rows.append(
            {
                "database": database,
                "total_models": total_models,
                "error_models": error_count,
                "error_rate_pct": round(error_count / total_models * 100, 1),
                "source_file": str(path),
            }
        )
    return pd.DataFrame(rows)


def build_fig2c() -> pd.DataFrame:
    files = {
        "CarveMe": TMP_DIR / "carveme_check.xlsx",
        "BiGG": TMP_DIR / "bigg_check.xlsx",
        "ModelSEED": TMP_DIR / "seed_check.xlsx",
        "VMH": TMP_DIR / "virtual_check.xlsx",
    }
    module_map = {
        "reducing_power": "Reducing power",
        "energy": "Energy",
        "metabolite": "Metabolite",
        "yield": "Yield",
        "biomass": "Biomass",
    }
    rows: list[dict[str, object]] = []
    for database, path in files.items():
        df = pd.read_excel(path)
        total_models = int(len(df))
        counts = {name: int(df[name].eq(0).sum()) for name in module_map}
        module_total = sum(counts.values())
        for col, label in module_map.items():
            count = counts[col]
            rows.append(
                {
                    "database": database,
                    "module": label,
                    "error_count": count,
                    "total_models": total_models,
                    "error_rate_pct": round(count / total_models * 100, 1),
                    "module_total_errors": module_total,
                    "normalized_share_pct": round(count / module_total * 100, 1) if module_total else 0.0,
                    "source_file": str(path),
                }
            )
    return pd.DataFrame(rows)


def build_fig2d() -> pd.DataFrame:
    error_types = {
        "bio_norxn": "does_not_exist",
        "bio_nogrow": "cannot_calculate_growth",
        "bio_coupling": "coupled_macromolecules",
        "bio_1g": "is_not_1g",
    }
    rows: list[dict[str, object]] = []
    for database, path in DATABASE_FILES.items():
        df = pd.read_excel(path)
        counts = {col: int(df[col].eq(0).sum()) for col in error_types}
        total_error_events = sum(counts.values())
        for col, label in error_types.items():
            count = counts[col]
            rows.append(
                {
                    "database": database,
                    "error_type": label,
                    "error_count": count,
                    "total_error_events": total_error_events,
                    "share_pct": round(count / total_error_events * 100, 1) if total_error_events else 0.0,
                    "source_file": str(path),
                }
            )
    return pd.DataFrame(rows)


def build_fig3(metrics: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for database, path in DATABASE_FILES.items():
        df = pd.read_excel(path)
        for _, row in df.iterrows():
            model = row.get("model", "")
            for metric in metrics:
                parsed = parse_tuple(row.get(metric), 2)
                if not parsed:
                    continue
                for stage, value in zip(("Initial", "Final"), (parsed[0], parsed[-1])):
                    rows.append(
                        {
                            "database": database,
                            "model": model,
                            "metric": metric,
                            "stage": stage,
                            "value": value,
                            "source_file": str(path),
                        }
                    )
    return pd.DataFrame(rows)


def build_fig4_distribution() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for database, path in DATABASE_FILES.items():
        df = pd.read_excel(path)
        for _, row in df.iterrows():
            model = row.get("model", "")
            parsed = parse_tuple(row.get("biomass"), 3)
            if not parsed:
                continue
            if database == "CarveMe" and parsed[2] >= 2:
                continue
            if database == "BiGG" and not parsed[0]:
                continue
            for stage, value in zip(("Initial", "Middle", "Final"), parsed[:3]):
                rows.append(
                    {
                        "database": database,
                        "model": model,
                        "stage": stage,
                        "growth_rate_h1": value,
                        "source_file": str(path),
                    }
                )
    return pd.DataFrame(rows)


def build_fig4_lines(database: str) -> pd.DataFrame:
    path = DATABASE_FILES[database]
    df = pd.read_excel(path)
    rows: list[dict[str, object]] = []
    order_index = 0
    for _, row in df.iterrows():
        model = row.get("model", "")
        parsed = parse_tuple(row.get("biomass"), 3)
        if not parsed:
            continue
        if database == "CarveMe" and parsed[2] >= 2:
            continue
        if database == "BiGG" and not parsed[0]:
            continue
        order_index += 1
        rows.append(
            {
                "order_index": order_index,
                "database": database,
                "model": model,
                "initial_growth_h1": parsed[0],
                "final_growth_h1": parsed[2],
                "source_file": str(path),
            }
        )
    return pd.DataFrame(rows)


def build_fig5a() -> pd.DataFrame:
    sheet = article_raw("Bacillus subtilis不同底物下的生长速率")
    condition_map = {
        "glc=7.63": "glc",
        "fru=5.72": "fru",
        "6pgc=5.13": "6pgc",
        "succ=3.35, glu__L=2.21": "succ&glu",
        "glyc=6.22": "glyc",
        "pyr=8.26": "pyr",
        "mal=26.51": "mal",
        "glc=5.95, mal=14.6": "glc&mal",
    }
    rows: list[dict[str, object]] = []

    for row_idx in range(3, 11):
        raw_condition = str(sheet.iat[row_idx, 0]).strip()
        condition = condition_map[raw_condition]
        experimental = as_float(sheet.iat[row_idx, 1])
        published_ibsu1147 = as_float(sheet.iat[row_idx, 2])
        published_ecbsu1 = as_float(sheet.iat[row_idx, 3])

        if pd.notna(sheet.iat[row_idx, 4]):
            rows.append(
                {
                    "database": "CarveMe",
                    "condition": condition,
                    "raw_condition": raw_condition,
                    "experimental_growth_h1": experimental,
                    "published_iBsu1147R_growth_h1": published_ibsu1147,
                    "published_ecBSU1_growth_h1": published_ecbsu1,
                    "control_model_growth_h1": as_float(sheet.iat[row_idx, 4]),
                    "initial_model_growth_h1": 0.0,
                    "experimental_status": "exact_from_uploaded_workbook",
                    "control_status": "exact_from_uploaded_workbook",
                    "initial_status": "inferred_zero_from_uploaded_workbook_layout",
                    "experimental_source": str(ARTICLE_WORKBOOK),
                    "control_source": str(ARTICLE_WORKBOOK),
                    "initial_source": str(ARTICLE_WORKBOOK),
                    "note": "Uploaded workbook only provides the final CarveMe column; initial values are filled as 0 based on the omitted initial column layout.",
                }
            )

        if pd.notna(sheet.iat[row_idx, 6]):
            rows.append(
                {
                    "database": "ModelSEED",
                    "condition": condition,
                    "raw_condition": raw_condition,
                    "experimental_growth_h1": experimental,
                    "published_iBsu1147R_growth_h1": published_ibsu1147,
                    "published_ecBSU1_growth_h1": published_ecbsu1,
                    "control_model_growth_h1": as_float(sheet.iat[row_idx, 6]),
                    "initial_model_growth_h1": 0.0,
                    "experimental_status": "exact_from_uploaded_workbook",
                    "control_status": "exact_from_uploaded_workbook",
                    "initial_status": "exact_from_legacy_script",
                    "experimental_source": str(ARTICLE_WORKBOOK),
                    "control_source": str(ARTICLE_WORKBOOK),
                    "initial_source": str(SOURCE_SCRIPT),
                    "note": "The uploaded workbook omits a ModelSEED initial column; zeros are backfilled from the legacy plotting script and remain consistent with the uploaded layout.",
                }
            )

        if pd.notna(sheet.iat[row_idx, 7]) and pd.notna(sheet.iat[row_idx, 8]):
            rows.append(
                {
                    "database": "BiGG",
                    "condition": condition,
                    "raw_condition": raw_condition,
                    "experimental_growth_h1": experimental,
                    "published_iBsu1147R_growth_h1": published_ibsu1147,
                    "published_ecBSU1_growth_h1": published_ecbsu1,
                    "control_model_growth_h1": as_float(sheet.iat[row_idx, 8]),
                    "initial_model_growth_h1": as_float(sheet.iat[row_idx, 7]),
                    "experimental_status": "exact_from_uploaded_workbook",
                    "control_status": "exact_from_uploaded_workbook",
                    "initial_status": "exact_from_uploaded_workbook",
                    "experimental_source": str(ARTICLE_WORKBOOK),
                    "control_source": str(ARTICLE_WORKBOOK),
                    "initial_source": str(ARTICLE_WORKBOOK),
                    "note": "",
                }
            )

        if pd.notna(sheet.iat[row_idx, 9]):
            rows.append(
                {
                    "database": "VMH",
                    "condition": condition,
                    "raw_condition": raw_condition,
                    "experimental_growth_h1": experimental,
                    "published_iBsu1147R_growth_h1": published_ibsu1147,
                    "published_ecBSU1_growth_h1": published_ecbsu1,
                    "control_model_growth_h1": as_float(sheet.iat[row_idx, 9]),
                    "initial_model_growth_h1": 0.0,
                    "experimental_status": "exact_from_uploaded_workbook",
                    "control_status": "exact_from_uploaded_workbook",
                    "initial_status": "inferred_zero_from_uploaded_workbook_layout",
                    "experimental_source": str(ARTICLE_WORKBOOK),
                    "control_source": str(ARTICLE_WORKBOOK),
                    "initial_source": str(ARTICLE_WORKBOOK),
                    "note": "Uploaded workbook only provides the final VMH column; initial values are filled as 0 based on the omitted initial column layout.",
                }
            )

    return pd.DataFrame(rows)


def build_fig5b() -> pd.DataFrame:
    sheet = article_raw("E. coli不同生长速率下的吸收和外排")
    measurements = ["GlucoseUptake", "O2uptake", "CO2production", "NH4uptake"]
    label_map = {
        "carveme": "CarveMe",
        "modelseed": "ModelSEED",
        "iECIAI1_1343": "BiGG",
        "virtual": "VMH",
    }
    block_offsets = {
        "carveme": 0,
        "modelseed": 5,
        "iECIAI1_1343": 10,
        "virtual": 15,
    }
    experimental_lookup: dict[tuple[float, str], float] = {}
    for row_idx in range(2, 8):
        growth_rate = as_float(sheet.iat[row_idx, 0])
        for metric_idx, measurement in enumerate(measurements):
            experimental_lookup[(growth_rate, measurement)] = as_float(sheet.iat[row_idx, metric_idx + 1])

    rows: list[dict[str, object]] = []
    for label, offset in block_offsets.items():
        database = label_map[label]
        for row_idx in range(16, 22):
            growth_rate = as_float(sheet.iat[row_idx, offset])
            initial_row_idx = row_idx + 9
            for metric_idx, measurement in enumerate(measurements):
                rows.append(
                    {
                        "database": database,
                        "model_label": label,
                        "growth_rate_h1": growth_rate,
                        "measurement": measurement,
                        "experimental_flux": experimental_lookup[(growth_rate, measurement)],
                        "control_model_flux": as_float(sheet.iat[row_idx, offset + metric_idx + 1]),
                        "initial_model_flux": as_float(sheet.iat[initial_row_idx, offset + metric_idx + 1]),
                        "experimental_status": "exact_from_uploaded_workbook",
                        "control_status": "exact_from_uploaded_workbook",
                        "initial_status": "exact_from_uploaded_workbook",
                        "experimental_source": str(ARTICLE_WORKBOOK),
                        "control_source": str(ARTICLE_WORKBOOK),
                        "initial_source": str(ARTICLE_WORKBOOK),
                        "source_file": str(ARTICLE_WORKBOOK),
                    }
                )
    return pd.DataFrame(rows)


def build_fig5c() -> pd.DataFrame:
    sheet = article_raw("E. coli 不同底物下的生长速率")
    metadata_rows = list(range(1, 9))
    simulation_rows = list(range(19, 27))

    rows: list[dict[str, object]] = []
    for meta_idx, sim_idx in zip(metadata_rows, simulation_rows):
        substrate = str(sheet.iat[meta_idx, 0]).strip()
        publication = str(sheet.iat[meta_idx, 1]).strip()
        experimental_uptake = as_float(sheet.iat[meta_idx, 2])
        experimental_growth = as_float(sheet.iat[meta_idx, 3])
        me_model_uptake = as_float(sheet.iat[meta_idx, 4])
        me_model_growth = as_float(sheet.iat[meta_idx, 5])

        rows.append(
            {
                "database": "CarveMe",
                "model_label": "carveme",
                "substrate": substrate,
                "publication": publication,
                "experimental_substrate_uptake_rate_mmol_gDW_h": experimental_uptake,
                "experimental_growth_h1": experimental_growth,
                "published_me_model_substrate_uptake_rate_mmol_gDW_h": me_model_uptake,
                "published_me_model_growth_h1": me_model_growth,
                "control_model_growth_h1": as_float(sheet.iat[sim_idx, 3]),
                "initial_model_growth_h1": as_float(sheet.iat[sim_idx, 2]),
                "experimental_status": "exact_from_uploaded_workbook",
                "control_status": "exact_from_uploaded_workbook",
                "initial_status": "exact_from_uploaded_workbook",
                "experimental_source": str(ARTICLE_WORKBOOK),
                "control_source": str(ARTICLE_WORKBOOK),
                "initial_source": str(ARTICLE_WORKBOOK),
                "source_file": str(ARTICLE_WORKBOOK),
                "note": "",
            }
        )
        rows.append(
            {
                "database": "BiGG",
                "model_label": "iECIAI1_1343",
                "substrate": substrate,
                "publication": publication,
                "experimental_substrate_uptake_rate_mmol_gDW_h": experimental_uptake,
                "experimental_growth_h1": experimental_growth,
                "published_me_model_substrate_uptake_rate_mmol_gDW_h": me_model_uptake,
                "published_me_model_growth_h1": me_model_growth,
                "control_model_growth_h1": as_float(sheet.iat[sim_idx, 5]),
                "initial_model_growth_h1": as_float(sheet.iat[sim_idx, 4]),
                "experimental_status": "exact_from_uploaded_workbook",
                "control_status": "exact_from_uploaded_workbook",
                "initial_status": "exact_from_uploaded_workbook",
                "experimental_source": str(ARTICLE_WORKBOOK),
                "control_source": str(ARTICLE_WORKBOOK),
                "initial_source": str(ARTICLE_WORKBOOK),
                "source_file": str(ARTICLE_WORKBOOK),
                "note": "",
            }
        )
        rows.append(
            {
                "database": "VMH",
                "model_label": "virtual",
                "substrate": substrate,
                "publication": publication,
                "experimental_substrate_uptake_rate_mmol_gDW_h": experimental_uptake,
                "experimental_growth_h1": experimental_growth,
                "published_me_model_substrate_uptake_rate_mmol_gDW_h": me_model_uptake,
                "published_me_model_growth_h1": me_model_growth,
                "control_model_growth_h1": as_float(sheet.iat[sim_idx, 7]),
                "initial_model_growth_h1": as_float(sheet.iat[sim_idx, 6]),
                "experimental_status": "exact_from_uploaded_workbook",
                "control_status": "exact_from_uploaded_workbook",
                "initial_status": "exact_from_uploaded_workbook",
                "experimental_source": str(ARTICLE_WORKBOOK),
                "control_source": str(ARTICLE_WORKBOOK),
                "initial_source": str(ARTICLE_WORKBOOK),
                "source_file": str(ARTICLE_WORKBOOK),
                "note": "",
            }
        )
        rows.append(
            {
                "database": "ModelSEED",
                "model_label": "modelseed",
                "substrate": substrate,
                "publication": publication,
                "experimental_substrate_uptake_rate_mmol_gDW_h": experimental_uptake,
                "experimental_growth_h1": experimental_growth,
                "published_me_model_substrate_uptake_rate_mmol_gDW_h": me_model_uptake,
                "published_me_model_growth_h1": me_model_growth,
                "control_model_growth_h1": as_float(sheet.iat[sim_idx, 9]),
                "initial_model_growth_h1": 0.0,
                "experimental_status": "exact_from_uploaded_workbook",
                "control_status": "exact_from_uploaded_workbook",
                "initial_status": "exact_from_legacy_script",
                "experimental_source": str(ARTICLE_WORKBOOK),
                "control_source": str(ARTICLE_WORKBOOK),
                "initial_source": str(SOURCE_SCRIPT),
                "source_file": str(ARTICLE_WORKBOOK),
                "note": "The uploaded workbook omits a ModelSEED initial column; zeros are backfilled from the legacy plotting script.",
            }
        )
    return pd.DataFrame(rows)


def autofit_workbook(path: Path) -> None:
    workbook = load_workbook(path)
    for sheet in workbook.worksheets:
        sheet.freeze_panes = "A2"
        sheet.auto_filter.ref = sheet.dimensions
        for column_cells in sheet.columns:
            values = ["" if cell.value is None else str(cell.value) for cell in column_cells]
            width = min(max(len(value) for value in values) + 2, 60)
            sheet.column_dimensions[column_cells[0].column_letter].width = width
    workbook.save(path)


def write_note() -> None:
    note = f"""# MQC图2-5原始数据整理说明

- 输出文件：`{OUTPUT_XLSX}`
- 生成脚本：`{PAPER_DIR / 'organize_mqc_figure_data.py'}`
- 补充来源：`{ARTICLE_WORKBOOK}`

## 说明

- `图2`、`图3`、`图4` 的数值表均可从 `tmp` 目录中的中间 Excel 结果精确恢复。
- `图2-b`、`图3-a`、`图3-d` 为示意/流程面板，没有独立的数值原始表。
- `图5a_Bacillus` 已优先切换到上传工作簿中的原始表：experimental、CarveMe final、ModelSEED final、BiGG initial/final、VMH final 均来自 `{ARTICLE_WORKBOOK}`。
- `图5a_Bacillus` 中 CarveMe/VMH 的 initial 列在上传工作簿里没有单独给出，我保留为 0，并在表内用 `inferred_zero_from_uploaded_workbook_layout` 单独标注来源级别。
- `图5a_Bacillus` 和 `图5c_Ecoli生长` 中 ModelSEED initial 全 0 由 `{SOURCE_SCRIPT}` 补足，并在表内单独标注。
- `图5b_Ecoli通量` 已切换为 `{ARTICLE_WORKBOOK}` 中的精确实验值、质控前值和质控后值。
- `图5c_Ecoli生长` 已切换为 `{ARTICLE_WORKBOOK}` 中的 experimental / CarveMe / BiGG / VMH / ModelSEED final 数据，并额外保留论文里的 publication 与 ME-model 对照列。
"""
    OUTPUT_NOTE.write_text(note, encoding="utf-8")


def main() -> None:
    sheets = {
        "说明": build_overview_sheet(),
        "图2a_总错误率": build_fig2a(),
        "图2c_模块分布": build_fig2c(),
        "图2d_biomass错误": build_fig2d(),
        "图3b_能量原始值": build_fig3(["ATP", "GTP", "CTP"]),
        "图3c_还原力原始值": build_fig3(["NADH", "NADPH", "Q8H2"]),
        "图4a_生长分布": build_fig4_distribution(),
        "图4b_CarveMe折线": build_fig4_lines("CarveMe"),
        "图4b_VMH折线": build_fig4_lines("VMH"),
        "图5a_Bacillus": build_fig5a(),
        "图5b_Ecoli通量": build_fig5b(),
        "图5c_Ecoli生长": build_fig5c(),
    }

    with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as writer:
        for sheet_name, dataframe in sheets.items():
            dataframe.to_excel(writer, sheet_name=sheet_name, index=False)

    autofit_workbook(OUTPUT_XLSX)
    write_note()
    print(f"wrote {OUTPUT_XLSX}")
    print(f"wrote {OUTPUT_NOTE}")


if __name__ == "__main__":
    main()
