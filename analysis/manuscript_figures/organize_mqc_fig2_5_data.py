#!/usr/bin/env python3
"""Organize source data for MQC figures 2-5 into a workbook."""

from __future__ import annotations

import ast
from pathlib import Path

import pandas as pd


BASE = Path("/home/dengxg/project/mqc/others/定量分析")
TMP = Path("/home/dengxg/project/mqc/tmp")
OUT_DIR = BASE / "result" / "论文"
OUT_XLSX = OUT_DIR / "MQC图2-5原始数据整理.xlsx"
OUT_MD = OUT_DIR / "MQC图2-5原始数据说明.md"


def parse_tuple(value):
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text == "无" or "nan" in text.lower():
        return None
    try:
        parsed = ast.literal_eval(text)
    except Exception:
        return None
    if isinstance(parsed, (list, tuple)):
        cleaned = []
        for item in parsed:
            if item in ("", None):
                cleaned.append(None)
            else:
                try:
                    cleaned.append(float(item))
                except Exception:
                    cleaned.append(None)
        return cleaned
    return None


def fig2a_overall_error():
    files = {
        "BiGG": TMP / "bigg_analysis_check.xlsx",
        "ModelSEED": TMP / "seed_analysis_check.xlsx",
        "VMH": TMP / "virtual_analysis_check.xlsx",
        "CarveMe": TMP / "web_CARVEME_COMMEN_check.xlsx",
    }
    rows = []
    for db, path in files.items():
        df = pd.read_excel(path)
        total = len(df)
        pass_count = int(df["all"].sum())
        error_count = total - pass_count
        rows.append(
            {
                "database": db,
                "total_models": total,
                "pass_models": pass_count,
                "error_models": error_count,
                "error_rate_pct": round(error_count / total * 100, 1),
                "source_file": str(path),
                "note": "对应 MQC图2 panel a",
            }
        )
    return pd.DataFrame(rows)


def fig2c_qc_module_distribution():
    files = {
        "BiGG": TMP / "bigg_check.xlsx",
        "ModelSEED": TMP / "seed_check.xlsx",
        "VMH": TMP / "virtual_check.xlsx",
        "CarveMe": TMP / "carveme_check.xlsx",
    }
    modules = [
        ("reducing_power", "Reducing Power"),
        ("energy", "Energy"),
        ("metabolite", "Metabolite"),
        ("yield", "Yield"),
        ("biomass", "Biomass"),
    ]
    rows = []
    for db, path in files.items():
        df = pd.read_excel(path)
        total = len(df)
        module_errors = {}
        for col, label in modules:
            error_count = int(total - df[col].sum())
            module_errors[label] = error_count
        total_error_pct_sum = sum((count / total * 100) for count in module_errors.values())
        for label in [m[1] for m in modules]:
            error_count = module_errors[label]
            raw_pct = round(error_count / total * 100, 1)
            normalized_pct = round(raw_pct / total_error_pct_sum * 100, 1)
            rows.append(
                {
                    "database": db,
                    "module": label,
                    "error_models": error_count,
                    "error_rate_pct_of_models": raw_pct,
                    "normalized_share_pct_for_plot": normalized_pct,
                    "source_file": str(path),
                    "note": "normalized_share_pct_for_plot 对应 MQC图2 panel c",
                }
            )
    return pd.DataFrame(rows)


def fig2d_biomass_error_breakdown():
    files = {
        "BiGG": TMP / "bigg_new.xlsx",
        "ModelSEED": TMP / "seed_new.xlsx",
        "VMH": TMP / "virtual_new.xlsx",
        "CarveMe": TMP / "carveme_new.xlsx",
    }
    cols = [
        ("bio_norxn", "does not exist"),
        ("bio_nogrow", "cannot calculate growth"),
        ("bio_coupling", "coupled macromolecules"),
        ("bio_1g", "is not 1g"),
    ]
    rows = []
    for db, path in files.items():
        df = pd.read_excel(path)
        counts = {label: int((df[col] == 0).sum()) for col, label in cols}
        total = sum(counts.values())
        for label in [c[1] for c in cols]:
            count = counts[label]
            rows.append(
                {
                    "database": db,
                    "error_type": label,
                    "error_models": count,
                    "normalized_share_pct_for_plot": round(count / total * 100, 1),
                    "source_file": str(path),
                    "note": "对应 MQC图2 panel d",
                }
            )
    return pd.DataFrame(rows)


def _metric_long_dataframe(metrics):
    files = {
        "CarveMe": TMP / "carveme_new.xlsx",
        "BiGG": TMP / "bigg_new.xlsx",
        "ModelSEED": TMP / "seed_new.xlsx",
        "VMH": TMP / "virtual_new.xlsx",
    }
    rows = []
    for db, path in files.items():
        df = pd.read_excel(path)
        for idx, row in df.iterrows():
            model_name = row.get("model")
            if pd.isna(model_name):
                model_name = f"{db}_model_{idx + 1}"
            for metric in metrics:
                parsed = parse_tuple(row.get(metric))
                if not parsed or len(parsed) < 2:
                    continue
                initial = parsed[0]
                final = parsed[-1]
                if initial is not None:
                    rows.append(
                        {
                            "database": db,
                            "model": model_name,
                            "metric": metric,
                            "stage": "Initial",
                            "value": initial,
                            "source_file": str(path),
                        }
                    )
                if final is not None:
                    rows.append(
                        {
                            "database": db,
                            "model": model_name,
                            "metric": metric,
                            "stage": "Final",
                            "value": final,
                            "source_file": str(path),
                        }
                    )
    return pd.DataFrame(rows)


def fig3_energy_raw():
    return _metric_long_dataframe(["ATP", "GTP", "CTP"])


def fig3_redox_raw():
    return _metric_long_dataframe(["NADH", "NADPH", "Q8H2"])


def fig4_growth_raw():
    files = {
        "CarveMe": TMP / "carveme_new.xlsx",
        "BiGG": TMP / "bigg_new.xlsx",
        "ModelSEED": TMP / "seed_new.xlsx",
        "VMH": TMP / "virtual_new.xlsx",
    }
    rows = []
    for db, path in files.items():
        df = pd.read_excel(path)
        for idx, row in df.iterrows():
            model_name = row.get("model")
            if pd.isna(model_name):
                model_name = f"{db}_model_{idx + 1}"
            parsed = parse_tuple(row.get("biomass"))
            if not parsed or len(parsed) < 3:
                continue
            initial, middle, final = parsed[0], parsed[1], parsed[2]
            if db == "CarveMe" and (final is None or final >= 2):
                continue
            if db == "BiGG" and not initial:
                continue
            for stage, value in [("Initial", initial), ("Middle", middle), ("Final", final)]:
                rows.append(
                    {
                        "database": db,
                        "model_order": idx + 1,
                        "model": model_name,
                        "stage": stage,
                        "value": value,
                        "source_file": str(path),
                        "note": "筛选逻辑与 tips_analysis.xlsx / MQC图4 一致",
                    }
                )
    return pd.DataFrame(rows)


def fig4_line_series(growth_df, database):
    wide = (
        growth_df[growth_df["database"] == database]
        .pivot_table(index=["model_order", "model"], columns="stage", values="value", aggfunc="first")
        .reset_index()
        .sort_values("model_order")
    )
    wide.columns.name = None
    wide["database"] = database
    wide["note"] = "对应 MQC图4 中每个模型的折线/序列图"
    return wide[["database", "model_order", "model", "Initial", "Final", "note"]]


def fig5_bacillus_data():
    rows = []

    def add_rows(database, labels, experimental, control, initial, source_file, source_status):
        for label, exp, ctrl, ini in zip(labels, experimental, control, initial):
            rows.append(
                {
                    "database": database,
                    "substrate": label,
                    "experimental_growth_rate_h-1": exp,
                    "control_model_growth_rate_h-1": ctrl,
                    "initial_model_growth_rate_h-1": ini,
                    "source_status": source_status,
                    "source_file": source_file,
                    "note": "对应 MQC图5 panel a（Bacillus）",
                }
            )

    add_rows(
        "ModelSEED",
        ["glc", "fru", "succ&glu", "glyc", "mal", "glc&mal"],
        [0.59, 0.53, 0.22, 0.40, 0.57, 0.75],
        [0.33, 0.24, 0.15, 0.14, 0.56, 0.82],
        [0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
        str(BASE / "figure_audit_workspace" / "remediation_scripts" / "定量分析.py"),
        "exact_from_notebook",
    )

    add_rows(
        "CarveMe",
        ["glc", "fru", "succ&glu", "glyc", "pyr", "mal", "glc&mal"],
        [0.59, 0.53, 0.22, 0.40, 0.17, 0.57, 0.75],
        [0.60, 0.44, 0.25, 0.29, 0.27, 0.93, 1.00],
        [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
        str(BASE / "result" / "定量数据" / "Bacillus_carveme.png"),
        "approx_recovered_from_panel_png",
    )

    add_rows(
        "BiGG",
        ["glc", "fru", "6pgc", "succ&glu", "glyc", "pyr", "mal", "glc&mal"],
        [0.59, 0.53, 0.42, 0.22, 0.40, 0.17, 0.57, 0.75],
        [0.65, 0.52, 0.45, 0.30, 0.32, 0.28, 0.65, 0.65],
        [0.62, 0.47, 0.40, 0.26, 0.28, 0.24, 0.62, 0.62],
        str(BASE / "result" / "定量数据" / "Bacillus_bigg.png"),
        "approx_recovered_from_panel_png",
    )

    add_rows(
        "VMH",
        ["glc", "fru", "succ&glu", "glyc", "mal", "glc&mal"],
        [0.59, 0.53, 0.22, 0.40, 0.57, 0.75],
        [0.43, 0.32, 0.10, 0.20, 0.69, 0.71],
        [0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
        str(BASE / "result" / "定量数据" / "Bacillus_virtual.png"),
        "approx_recovered_from_panel_png",
    )

    return pd.DataFrame(rows)


def fig5_ecoli_flux_data():
    growth_conditions = [0.15, 0.20, 0.25, 0.30, 0.35, 0.40]
    measurements = ["GlucoseUptake", "O2uptake", "CO2production", "NH4uptake"]
    experimental = [
        1.59, 3.67, 3.98, 1.26,
        2.08, 5.49, 5.36, 1.70,
        2.79, 6.18, 6.61, 2.26,
        3.01, 7.41, 7.38, 2.51,
        3.36, 7.57, 7.73, 2.85,
        4.01, 8.23, 8.68, 3.39,
    ]
    control_data = {
        "CarveMe": [
            1.52, 2.68, 2.84, 1.63,
            2.02, 3.57, 3.79, 2.17,
            2.53, 4.47, 4.74, 2.71,
            3.04, 5.36, 5.68, 3.25,
            3.54, 6.25, 6.63, 3.80,
            4.05, 7.15, 7.58, 4.34,
        ],
        "ModelSEED": [
            2.49, 0.00, 3.69, 0.00,
            3.23, 0.00, 4.93, 0.00,
            4.15, 0.01, 6.16, 0.00,
            4.98, 0.01, 7.39, 0.00,
            5.81, 0.01, 8.62, 0.00,
            6.64, 0.01, 9.85, 0.00,
        ],
        "VMH": [
            1.42, 2.07, 2.35, 1.62,
            1.89, 2.76, 3.13, 2.16,
            2.36, 3.46, 3.91, 2.70,
            2.83, 4.15, 4.70, 3.23,
            3.30, 4.84, 5.48, 3.77,
            3.78, 5.53, 6.26, 4.31,
        ],
        "BiGG": [
            1.53, 2.54, 2.86, 1.62,
            2.04, 3.39, 3.81, 2.16,
            2.55, 4.23, 4.77, 2.70,
            3.06, 5.08, 5.72, 3.24,
            3.57, 5.93, 6.67, 3.78,
            4.09, 6.77, 7.63, 4.32,
        ],
    }
    initial_data = {
        "CarveMe": [
            0.00, 2.08, 1.07, 0.73,
            0.00, 2.78, 1.43, 0.97,
            0.00, 3.47, 1.79, 1.22,
            0.00, 4.16, 2.14, 1.46,
            0.00, 4.86, 2.50, 1.71,
            0.00, 5.55, 2.86, 1.95,
        ],
        "ModelSEED": [
            0.00, 0.00, 0.93, 0.00,
            0.00, 0.00, 1.24, 0.00,
            0.00, 0.00, 1.56, 0.00,
            0.00, 0.00, 1.87, 0.00,
            0.00, 0.00, 2.18, 0.00,
            0.00, 0.00, 2.49, 0.00,
        ],
        "VMH": [
            0.00, 1.22, 1.05, 0.00,
            0.00, 1.62, 1.40, 0.00,
            0.00, 2.03, 1.76, 0.00,
            0.00, 2.44, 2.11, 0.00,
            0.00, 2.84, 2.46, 0.00,
            0.00, 3.25, 2.81, 0.00,
        ],
        "BiGG": [
            0.50, 4.02, 4.34, 1.62,
            0.69, 5.06, 5.49, 2.16,
            0.88, 6.10, 6.63, 2.70,
            1.08, 7.14, 7.78, 3.24,
            1.27, 8.18, 8.93, 3.78,
            1.46, 9.22, 10.07, 4.32,
        ],
    }
    rows = []
    source_file = str(BASE / "figure_audit_workspace" / "remediation_scripts" / "定量分析.py")
    for database in ["CarveMe", "BiGG", "ModelSEED", "VMH"]:
        idx = 0
        for growth_condition in growth_conditions:
            for measurement in measurements:
                rows.append(
                    {
                        "database": database,
                        "growth_condition_h-1": growth_condition,
                        "measurement": measurement,
                        "experimental_value": experimental[idx],
                        "control_model_value": control_data[database][idx],
                        "initial_model_value": initial_data[database][idx],
                        "source_file": source_file,
                        "note": "对应 MQC图5 panel b（E. coli uptake/secretion）",
                    }
                )
                idx += 1
    return pd.DataFrame(rows)


def fig5_ecoli_growth_data():
    substrates = [
        "Acetate",
        "Succinate",
        "L-Lactate",
        "Pyruvate",
        "2-Oxoglutarate",
        "D-Glucose",
        "L-Malate",
        "Glycerol",
    ]
    experimental = [
        0.019677419,
        0.034836066,
        0.041618497,
        0.032110092,
        0.040540541,
        0.080703625,
        0.042608696,
        0.051301685,
    ]
    control_data = {
        "CarveMe": [
            0.0260382591074756,
            0.0512294450120816,
            0.0419465690665474,
            0.0,
            0.0641661735504868,
            0.1015828665533370,
            0.0470896918797927,
            0.0584871049852542,
        ],
        "BiGG": [
            0.0260108452194865,
            0.0527081724424213,
            0.0411109890250687,
            0.0371956567369668,
            0.0630007652282947,
            0.0994522862169264,
            0.0462541105106966,
            0.0575691236474373,
        ],
        "VMH": [
            0.0285088476588867,
            0.0586314036758236,
            0.0486593960333508,
            0.0,
            0.0670146008850796,
            0.1069296585056840,
            0.0543281813876899,
            0.0613209176059072,
        ],
        "ModelSEED": [
            0.0253492866188815,
            0.0547314142907671,
            0.0426328911317554,
            0.0362955694770347,
            0.0655750704045048,
            0.0949482342107753,
            0.0483940926360467,
            0.0535791739899091,
        ],
    }
    initial_data = {
        "CarveMe": [
            0.0329734051643870,
            0.0642113679516513,
            0.0555341560662965,
            0.0,
            0.0885075612306545,
            0.1301581782803730,
            0.0642113679516513,
            0.0711531374599405,
        ],
        "BiGG": [
            0.0,
            0.0205545991234181,
            0.0102772995617073,
            0.0102772995617073,
            0.0308318986851254,
            0.0719410969319120,
            0.0205545991234181,
            0.0411091982467866,
        ],
        "VMH": [
            0.0472170099561749,
            0.0708255149342583,
            0.0708255149342585,
            0.0,
            0.0826297674232919,
            0.1416510298685000,
            0.0708255149342698,
            0.0826297674232922,
        ],
        "ModelSEED": [0.0] * 8,
    }
    rows = []
    source_file = str(BASE / "figure_audit_workspace" / "remediation_scripts" / "定量分析.py")
    for database in ["CarveMe", "BiGG", "ModelSEED", "VMH"]:
        for substrate, exp, ctrl, ini in zip(
            substrates,
            experimental,
            control_data[database],
            initial_data[database],
        ):
            rows.append(
                {
                    "database": database,
                    "substrate": substrate,
                    "experimental_growth_rate_h-1": exp,
                    "control_model_growth_rate_h-1": ctrl,
                    "initial_model_growth_rate_h-1": ini,
                    "source_file": source_file,
                    "note": "对应 MQC图5 panel c（E. coli growth）",
                }
            )
    return pd.DataFrame(rows)


def build_readme():
    rows = [
        {
            "figure": "MQC图2",
            "panel": "a",
            "content": "Overall Error Rate by Database",
            "status": "exact_from_table",
            "source": str(TMP / "bigg_analysis_check.xlsx"),
            "details": "BiGG/ModelSEED/VMH/CarveMe 分别来自 bigg_analysis_check.xlsx、seed_analysis_check.xlsx、virtual_analysis_check.xlsx、web_CARVEME_COMMEN_check.xlsx 的 all 列。",
        },
        {
            "figure": "MQC图2",
            "panel": "b",
            "content": "Common Errors in Genome-scale network models",
            "status": "no_numeric_source",
            "source": str(BASE / "result" / "其他" / "架构图3.png"),
            "details": "示意图/概念图，没有对应数值型原始表。",
        },
        {
            "figure": "MQC图2",
            "panel": "c",
            "content": "Error Distribution by QC Module",
            "status": "exact_from_table",
            "source": str(TMP / "bigg_check.xlsx"),
            "details": "按 *_check.xlsx 各模块错误率重新归一化后得到，和图中百分比一致。",
        },
        {
            "figure": "MQC图2",
            "panel": "d",
            "content": "Biomass Error Types Distribution",
            "status": "exact_from_table",
            "source": str(TMP / "bigg_new.xlsx"),
            "details": "按 *_new.xlsx 中 bio_norxn / bio_nogrow / bio_coupling / bio_1g 统计 0 值并归一化。",
        },
        {
            "figure": "MQC图3",
            "panel": "a,d",
            "content": "cycpFBA / reaction-penalty schematic",
            "status": "no_numeric_source",
            "source": str(BASE / "result" / "其他" / "循环罚分法.png"),
            "details": "流程/规则示意，没有对应数值型原始表。",
        },
        {
            "figure": "MQC图3",
            "panel": "b",
            "content": "ATP/GTP/CTP 分布箱线图",
            "status": "exact_from_table",
            "source": str(TMP / "carveme_new.xlsx"),
            "details": "从四个 *_new.xlsx 的 ATP/GTP/CTP 元组列解析 Initial/Final。",
        },
        {
            "figure": "MQC图3",
            "panel": "c",
            "content": "NADH/NADPH/Q8H2 分布箱线图",
            "status": "exact_from_table",
            "source": str(TMP / "carveme_new.xlsx"),
            "details": "从四个 *_new.xlsx 的 NADH/NADPH/Q8H2 元组列解析 Initial/Final。",
        },
        {
            "figure": "MQC图4",
            "panel": "a",
            "content": "Initial/Middle/Final 生长分布",
            "status": "exact_from_table",
            "source": str(TMP / "tips_analysis.xlsx"),
            "details": "底层来自四个 *_new.xlsx 的 biomass 元组列；已按 notebook 的筛选逻辑整理。",
        },
        {
            "figure": "MQC图4",
            "panel": "b",
            "content": "CarveMe / VMH 单模型序列图",
            "status": "exact_from_table",
            "source": str(TMP / "carveme_new.xlsx"),
            "details": "按 biomass 元组列提取 Initial/Final，保留原文件顺序。",
        },
        {
            "figure": "MQC图5",
            "panel": "a",
            "content": "Bacillus 不同底物生长率散点图",
            "status": "mixed_exact_and_approx",
            "source": str(BASE / "result" / "定量数据" / "Bacillus_modelseed.png"),
            "details": "ModelSEED 数组能从 notebook 精确恢复；CarveMe/BiGG/VMH 当前仓库只保留单 panel PNG，因此做了近似转录并单独标注 source_status。",
        },
        {
            "figure": "MQC图5",
            "panel": "b",
            "content": "E. coli uptake/secretion 对比",
            "status": "exact_from_notebook",
            "source": str(BASE / "figure_audit_workspace" / "remediation_scripts" / "定量分析.py"),
            "details": "24 个点按 6 个 growth condition × 4 个 flux type 组织。",
        },
        {
            "figure": "MQC图5",
            "panel": "c",
            "content": "E. coli 不同底物生长率对比",
            "status": "exact_from_notebook",
            "source": str(BASE / "figure_audit_workspace" / "remediation_scripts" / "定量分析.py"),
            "details": "8 个底物的 experimental / control / initial 数组可直接从 notebook 恢复。",
        },
    ]
    return pd.DataFrame(rows)


def write_markdown(readme_df):
    lines = [
        "# MQC图2-5原始数据说明",
        "",
        f"- 输出工作簿: `{OUT_XLSX}`",
        "- 说明: 图2/图3中的概念示意 panel 没有对应数值型原始数据；图5 panel a 中除 ModelSEED 外，其余数据库数值为从单 panel PNG 近似转录。",
        "",
        "## 面板来源总览",
        "",
        "| Figure | Panel | 内容 | 状态 | 说明 |",
        "| --- | --- | --- | --- | --- |",
    ]
    for _, row in readme_df.iterrows():
        lines.append(
            f"| {row['figure']} | {row['panel']} | {row['content']} | {row['status']} | {row['details']} |"
        )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")


def main():
    readme = build_readme()
    fig2a = fig2a_overall_error()
    fig2c = fig2c_qc_module_distribution()
    fig2d = fig2d_biomass_error_breakdown()
    fig3b = fig3_energy_raw()
    fig3c = fig3_redox_raw()
    fig4a = fig4_growth_raw()
    fig4b_carveme = fig4_line_series(fig4a, "CarveMe")
    fig4b_vmh = fig4_line_series(fig4a, "VMH")
    fig5a = fig5_bacillus_data()
    fig5b = fig5_ecoli_flux_data()
    fig5c = fig5_ecoli_growth_data()

    with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as writer:
        readme.to_excel(writer, sheet_name="说明", index=False)
        fig2a.to_excel(writer, sheet_name="图2a_总错误率", index=False)
        fig2c.to_excel(writer, sheet_name="图2c_模块分布", index=False)
        fig2d.to_excel(writer, sheet_name="图2d_biomass错误", index=False)
        fig3b.to_excel(writer, sheet_name="图3b_能量原始值", index=False)
        fig3c.to_excel(writer, sheet_name="图3c_还原力原始值", index=False)
        fig4a.to_excel(writer, sheet_name="图4a_生长分布", index=False)
        fig4b_carveme.to_excel(writer, sheet_name="图4b_CarveMe序列", index=False)
        fig4b_vmh.to_excel(writer, sheet_name="图4b_VMH序列", index=False)
        fig5a.to_excel(writer, sheet_name="图5a_Bacillus", index=False)
        fig5b.to_excel(writer, sheet_name="图5b_Ecoli通量", index=False)
        fig5c.to_excel(writer, sheet_name="图5c_Ecoli生长", index=False)

    write_markdown(readme)
    print(f"written: {OUT_XLSX}")
    print(f"written: {OUT_MD}")


if __name__ == "__main__":
    main()
