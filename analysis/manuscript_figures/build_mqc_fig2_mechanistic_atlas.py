from __future__ import annotations

import ast
import math
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Ellipse
from openpyxl import load_workbook


PAPER_DIR = Path("/home/dengxg/project/mqc/others/定量分析/result/论文")
TMP_DIR = Path("/home/dengxg/project/mqc/tmp")

ORGANIZED_XLSX = PAPER_DIR / "MQC图2-5原始数据整理.xlsx"
ORGANIZED_NOTE = PAPER_DIR / "MQC图2-5原始数据说明.md"
HOME_XLSX = Path("/home/dengxg/MQC图2-5原始数据整理.xlsx")
HOME_NOTE = Path("/home/dengxg/MQC图2-5原始数据说明.md")

INPUT_FILES = {
    "BiGG": TMP_DIR / "bigg_new.xlsx",
    "VMH": TMP_DIR / "virtual_new.xlsx",
    "ModelSEED": TMP_DIR / "seed_new.xlsx",
    "CarveMe": TMP_DIR / "carveme_new.xlsx",
}
DB_ORDER = ["BiGG", "VMH", "ModelSEED", "CarveMe"]
EXPECTED_COUNTS = {"BiGG": 73, "VMH": 806, "ModelSEED": 79, "CarveMe": 5585}

ERROR_PREVALENCE_FILES = {
    "BiGG": TMP_DIR / "bigg_analysis_check.xlsx",
    "VMH": TMP_DIR / "virtual_analysis_check.xlsx",
    "ModelSEED": TMP_DIR / "seed_analysis_check.xlsx",
    "CarveMe": TMP_DIR / "web_CARVEME_COMMEN_check.xlsx",
}

DB_COLORS = {
    "BiGG": "#2F6FA3",
    "VMH": "#C65B66",
    "ModelSEED": "#4F9D69",
    "CarveMe": "#D9892B",
}
CLUSTER_PALETTE = [
    "#355C7D",
    "#D96C4B",
    "#4F9D69",
    "#8A6F9E",
    "#E2A72E",
    "#3E8C8C",
    "#C65B84",
    "#6E6A5E",
]

COFACTOR_COLS = [
    "NADH",
    "NADPH",
    "FADH2",
    "FMNH2",
    "Q8H2",
    "MQL8",
    "DMMQL8",
    "ATP",
    "CTP",
    "GTP",
    "UTP",
    "ITP",
]
MATERIAL_COLS = ["metabolite_violation", "yield_violation"]
BIOMASS_COLS = [
    "biomass_gap",
    "gap_num_log1p",
    "bio_norxn_fail",
    "bio_nogrow_fail",
    "bio_coupling_fail",
    "bio_1g_fail",
    "ini_weight_deviation",
    "carbon_biomass_ratio",
]
FEATURE_COLS = COFACTOR_COLS + MATERIAL_COLS + BIOMASS_COLS
BINARY_COLS = {
    "metabolite_violation",
    "yield_violation",
    "biomass_gap",
    "bio_norxn_fail",
    "bio_nogrow_fail",
    "bio_coupling_fail",
    "bio_1g_fail",
}
CONTINUOUS_COLS = [col for col in FEATURE_COLS if col not in BINARY_COLS]

FEATURE_LABELS = {
    "NADH": "NADH",
    "NADPH": "NADPH",
    "FADH2": "FADH2",
    "FMNH2": "FMNH2",
    "Q8H2": "Q8H2",
    "MQL8": "MQL8",
    "DMMQL8": "DMMQL8",
    "ATP": "ATP",
    "CTP": "CTP",
    "GTP": "GTP",
    "UTP": "UTP",
    "ITP": "ITP",
    "metabolite_violation": "Metabolite",
    "yield_violation": "Yield",
    "biomass_gap": "Gap-filled precursor",
    "gap_num_log1p": "Gap count",
    "bio_norxn_fail": "No biomass rxn",
    "bio_nogrow_fail": "No biomass growth",
    "bio_coupling_fail": "Biomass uncoupled",
    "bio_1g_fail": "Not 1 gDW",
    "ini_weight_deviation": "Biomass wt dev.",
    "carbon_biomass_ratio": "C-bio ratio",
}


def parse_second_value(value: object) -> float:
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


def yes_no_to_int(value: object) -> int:
    if pd.isna(value):
        return 0
    return int(str(value).strip().lower() in {"yes", "y", "true", "1"})


def safe_float(value: object, default: float = 0.0) -> float:
    try:
        out = float(value)
    except Exception:
        return default
    if math.isnan(out) or math.isinf(out):
        return default
    return out


def build_feature_matrix() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for database in DB_ORDER:
        path = INPUT_FILES[database]
        df = pd.read_excel(path)
        for idx, row in df.iterrows():
            item: dict[str, object] = {
                "database": database,
                "row_id": f"{database}_{idx}",
                "model": row.get("model"),
                "source_file": str(path),
            }
            for col in COFACTOR_COLS:
                item[col] = parse_second_value(row.get(col))
            item["metabolite_violation"] = yes_no_to_int(row.get("metabolite"))
            item["yield_violation"] = yes_no_to_int(row.get("yield"))
            item["biomass_gap"] = yes_no_to_int(row.get("is_gap"))
            item["gap_num_log1p"] = math.log1p(max(safe_float(row.get("gap_num")), 0.0))
            item["bio_norxn_fail"] = 1 - int(safe_float(row.get("bio_norxn"), 1.0))
            item["bio_nogrow_fail"] = 1 - int(safe_float(row.get("bio_nogrow"), 1.0))
            item["bio_coupling_fail"] = 1 - int(safe_float(row.get("bio_coupling"), 1.0))
            item["bio_1g_fail"] = 1 - int(safe_float(row.get("bio_1g"), 1.0))
            item["ini_weight_deviation"] = abs(safe_float(row.get("ini_weight(g)"), 1.0) - 1.0)
            item["carbon_biomass_ratio"] = safe_float(row.get("c_bio"), 0.0)
            rows.append(item)
    feature_df = pd.DataFrame(rows)
    feature_df[FEATURE_COLS] = feature_df[FEATURE_COLS].fillna(0.0)
    return feature_df


def load_error_prevalence() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for database in DB_ORDER:
        path = ERROR_PREVALENCE_FILES[database]
        df = pd.read_excel(path)
        total = int(len(df))
        failing = int(df["all"].eq(0).sum())
        rows.append(
            {
                "database": database,
                "total_models": total,
                "failing_models": failing,
                "error_rate_pct": failing / total * 100,
                "source_file": str(path),
            }
        )
    return pd.DataFrame(rows)


def standardize_matrix(values: np.ndarray) -> tuple[np.ndarray, pd.DataFrame]:
    mean = values.mean(axis=0)
    std = values.std(axis=0, ddof=0)
    std[std == 0] = 1.0
    scaled = (values - mean) / std
    stats = pd.DataFrame(
        {
            "feature": FEATURE_COLS,
            "transform": "zscore",
            "raw_min": values.min(axis=0),
            "raw_median": np.median(values, axis=0),
            "raw_max": values.max(axis=0),
            "clip_lower": np.nan,
            "clip_upper": np.nan,
            "center": mean,
            "scale": std,
            "scaled_min": scaled.min(axis=0),
            "scaled_max": scaled.max(axis=0),
        }
    )
    return scaled, stats


def refined_scale(feature_df: pd.DataFrame) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
    raw = feature_df[FEATURE_COLS].astype(float)
    scaled = pd.DataFrame(index=raw.index)
    stats_rows: list[dict[str, object]] = []

    for col in FEATURE_COLS:
        series = raw[col]
        if col in CONTINUOUS_COLS:
            lower = float(series.quantile(0.01))
            upper = float(series.quantile(0.99))
            clipped = series.clip(lower=lower, upper=upper)
            center = float(clipped.median())
            q25 = float(clipped.quantile(0.25))
            q75 = float(clipped.quantile(0.75))
            scale = q75 - q25
            if scale <= 1e-12:
                scale = float(clipped.std(ddof=0))
            if scale <= 1e-12:
                scale = 1.0
            transformed = ((clipped - center) / scale).clip(lower=-3.0, upper=3.0)
            transform = "winsorized_robust"
        else:
            lower = 0.0
            upper = 1.0
            center = float(series.mean())
            scale = 1.0
            transformed = series - center
            transform = "mean_centered_binary"
        scaled[col] = transformed
        stats_rows.append(
            {
                "feature": col,
                "feature_label": FEATURE_LABELS[col],
                "feature_type": "binary" if col in BINARY_COLS else "continuous",
                "transform": transform,
                "raw_min": float(series.min()),
                "raw_p01": float(series.quantile(0.01)),
                "raw_median": float(series.median()),
                "raw_p99": float(series.quantile(0.99)),
                "raw_max": float(series.max()),
                "clip_lower": lower,
                "clip_upper": upper,
                "center": center,
                "scale": scale,
                "scaled_min": float(transformed.min()),
                "scaled_max": float(transformed.max()),
            }
        )
    return scaled.to_numpy(dtype=float), pd.DataFrame(stats_rows), scaled


def pca_project(x: np.ndarray, n_components: int = 2) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    mean = x.mean(axis=0, keepdims=True)
    x_centered = x - mean
    _, s, vt = np.linalg.svd(x_centered, full_matrices=False)
    coords = x_centered @ vt[:n_components].T
    eigenvalues = (s**2) / max(len(x) - 1, 1)
    explained = eigenvalues / eigenvalues.sum()
    return coords, vt[:n_components], explained[:n_components], mean


def init_kmeans_pp(x: np.ndarray, k: int, rng: np.random.Generator) -> np.ndarray:
    centers = np.empty((k, x.shape[1]), dtype=float)
    first = rng.integers(0, len(x))
    centers[0] = x[first]
    closest_sq = np.sum((x - centers[0]) ** 2, axis=1)
    for i in range(1, k):
        total = closest_sq.sum()
        if total <= 1e-12:
            idx = rng.integers(0, len(x))
        else:
            idx = rng.choice(len(x), p=closest_sq / total)
        centers[i] = x[idx]
        dist_sq = np.sum((x - centers[i]) ** 2, axis=1)
        closest_sq = np.minimum(closest_sq, dist_sq)
    return centers


def assign_clusters(x: np.ndarray, centers: np.ndarray) -> tuple[np.ndarray, float]:
    dist_sq = ((x[:, None, :] - centers[None, :, :]) ** 2).sum(axis=2)
    labels = dist_sq.argmin(axis=1)
    inertia = float(dist_sq[np.arange(len(x)), labels].sum())
    return labels, inertia


def fit_kmeans(
    x: np.ndarray,
    k: int,
    n_init: int = 12,
    max_iter: int = 100,
    rng_seed: int = 42,
) -> tuple[np.ndarray, np.ndarray, float]:
    best_labels: np.ndarray | None = None
    best_centers: np.ndarray | None = None
    best_inertia = np.inf
    rng = np.random.default_rng(rng_seed)
    for init_idx in range(n_init):
        centers = init_kmeans_pp(x, k, np.random.default_rng(rng.integers(1_000_000_000) + init_idx))
        for _ in range(max_iter):
            labels, _ = assign_clusters(x, centers)
            new_centers = centers.copy()
            for cluster_idx in range(k):
                members = x[labels == cluster_idx]
                if len(members) == 0:
                    new_centers[cluster_idx] = x[rng.integers(0, len(x))]
                else:
                    new_centers[cluster_idx] = members.mean(axis=0)
            if np.allclose(new_centers, centers, atol=1e-6):
                centers = new_centers
                break
            centers = new_centers
        labels, inertia = assign_clusters(x, centers)
        if inertia < best_inertia:
            best_labels = labels.copy()
            best_centers = centers.copy()
            best_inertia = inertia
    if best_labels is None or best_centers is None:
        raise RuntimeError("k-means failed to initialize")
    return best_labels, best_centers, best_inertia


def silhouette_score_sample(
    x: np.ndarray,
    labels: np.ndarray,
    sample_n: int = 1800,
    rng_seed: int = 42,
) -> float:
    unique = np.unique(labels)
    if len(unique) < 2:
        return -1.0
    rng = np.random.default_rng(rng_seed)
    if len(x) > sample_n:
        sample_idx = rng.choice(len(x), size=sample_n, replace=False)
        sample = x[sample_idx]
        sample_labels = labels[sample_idx]
    else:
        sample = x
        sample_labels = labels
    dist = np.sqrt(np.maximum(((sample[:, None, :] - sample[None, :, :]) ** 2).sum(axis=2), 0.0))
    silhouettes: list[float] = []
    for i in range(len(sample)):
        same = sample_labels == sample_labels[i]
        same[i] = False
        a = float(dist[i, same].mean()) if same.any() else 0.0
        b = np.inf
        for lab in unique:
            if lab == sample_labels[i]:
                continue
            other = sample_labels == lab
            if other.any():
                b = min(b, float(dist[i, other].mean()))
        denom = max(a, b)
        silhouettes.append((b - a) / denom if denom > 0 else 0.0)
    return float(np.mean(silhouettes))


def choose_k(
    x: np.ndarray,
    mode: str,
    min_cluster_threshold: int | None,
    fallback_k: int = 5,
) -> tuple[int, pd.DataFrame, str]:
    records: list[dict[str, object]] = []
    for k in range(4, 9):
        labels, _, inertia = fit_kmeans(x, k, n_init=8, max_iter=80, rng_seed=3100 + k)
        sizes = pd.Series(labels).value_counts().sort_index()
        min_size = int(sizes.min())
        score = silhouette_score_sample(x, labels, rng_seed=4100 + k)
        eligible = min_cluster_threshold is None or min_size >= min_cluster_threshold
        records.append(
            {
                "mode": mode,
                "k": k,
                "silhouette": score,
                "inertia": inertia,
                "min_cluster_size": min_size,
                "eligible": eligible,
                "chosen": False,
            }
        )
    selection_df = pd.DataFrame(records)
    if min_cluster_threshold is None:
        best_score = float(selection_df["silhouette"].max())
        near_best = selection_df[selection_df["silhouette"] >= best_score - 0.01]
        chosen_k = int(near_best.sort_values("k").iloc[0]["k"])
        reason = "highest silhouette; smaller k preferred within 0.01"
    else:
        eligible_df = selection_df[selection_df["eligible"]]
        if eligible_df.empty:
            chosen_k = fallback_k
            reason = f"no k passed min cluster size threshold; fallback to k={fallback_k}"
        else:
            best_score = float(eligible_df["silhouette"].max())
            near_best = eligible_df[eligible_df["silhouette"] >= best_score - 0.01]
            chosen_k = int(near_best.sort_values("k").iloc[0]["k"])
            reason = "eligible k with highest silhouette; smaller k preferred within 0.01"
    selection_df.loc[selection_df["k"] == chosen_k, "chosen"] = True
    selection_df["selection_reason"] = reason
    selection_df["min_cluster_threshold"] = min_cluster_threshold if min_cluster_threshold is not None else np.nan
    return chosen_k, selection_df, reason


def describe_cluster(center: pd.Series) -> str:
    energy = float(center[["ATP", "CTP", "GTP", "UTP", "ITP"]].mean())
    redox = float(center[["NADH", "NADPH", "FADH2", "FMNH2", "Q8H2", "MQL8", "DMMQL8"]].mean())
    cofactor = max(energy, redox)
    gap_burden = float(center[["biomass_gap", "gap_num_log1p"]].mean())
    growth = float(center[["biomass_gap", "gap_num_log1p", "bio_nogrow_fail", "bio_1g_fail"]].mean())
    biomass_structure = float(
        center[["bio_norxn_fail", "bio_coupling_fail", "ini_weight_deviation", "carbon_biomass_ratio"]].max()
    )
    material = float(center[["metabolite_violation", "yield_violation"]].mean())
    if cofactor > 0.65 and growth > 0.18:
        return "Cofactor + growth defect"
    if cofactor > 0.65:
        return "Cofactor-dominant"
    if growth > 0.35 or biomass_structure > 0.65:
        return "Biomass-defect enriched"
    if gap_burden > 0.25:
        return "Gap-burden enriched"
    if cofactor > 0.22:
        return "Mild cofactor burden"
    if growth > 0.35:
        return "Biomass-defect enriched"
    if biomass_structure > 0.35:
        return "Biomass-structure enriched"
    if material > 0.12:
        return "Material/yield enriched"
    if float(center.abs().max()) < 0.35:
        return "Low-intensity baseline"
    return "Mixed biomass/gap burden"


def make_unique_descriptions(centroids: pd.DataFrame) -> dict[str, str]:
    descriptions: dict[str, str] = {}
    seen: dict[str, int] = {}
    for cluster_name, row in centroids.iterrows():
        desc = describe_cluster(row)
        if desc in seen:
            top_feature = row.sort_values(ascending=False).index[0]
            desc = f"{desc}: {FEATURE_LABELS[top_feature]}"
        seen[desc] = seen.get(desc, 0) + 1
        descriptions[cluster_name] = desc
    return descriptions


def build_cluster_outputs(
    feature_df: pd.DataFrame,
    x: np.ndarray,
    mode: str,
    min_cluster_threshold: int | None,
) -> dict[str, object]:
    chosen_k, k_selection, selection_reason = choose_k(x, mode, min_cluster_threshold)
    labels_raw, centers_raw, _ = fit_kmeans(x, chosen_k, n_init=22, max_iter=130, rng_seed=2026)

    coords, pca_components, explained, pca_mean = pca_project(x, n_components=2)
    center_coords = (centers_raw - pca_mean) @ pca_components[:2].T
    order = np.argsort(center_coords[:, 0])

    raw_to_cluster = {int(raw_idx): f"C{new_idx + 1}" for new_idx, raw_idx in enumerate(order)}
    ordered_centers = centers_raw[order]
    centroids = pd.DataFrame(ordered_centers, columns=FEATURE_COLS, index=[f"C{i + 1}" for i in range(chosen_k)])
    descriptions = make_unique_descriptions(centroids)

    assignment = feature_df[["database", "row_id", "model", "source_file"]].copy()
    assignment["cluster_raw"] = labels_raw + 1
    assignment["cluster"] = [raw_to_cluster[int(label)] for label in labels_raw]
    assignment["cluster_description"] = assignment["cluster"].map(descriptions)
    assignment["PC1"] = coords[:, 0]
    assignment["PC2"] = coords[:, 1]

    ordered_cluster_names = centroids.index.tolist()
    size_map = assignment["cluster"].value_counts().reindex(ordered_cluster_names).fillna(0).astype(int).to_dict()
    centroid_records = []
    for cluster_name, row in centroids.iterrows():
        top = row.sort_values(ascending=False).head(6)
        centroid_records.append(
            {
                "mode": mode,
                "cluster": cluster_name,
                "cluster_description": descriptions[cluster_name],
                "cluster_size": size_map[cluster_name],
                "top_positive_features": "; ".join(f"{FEATURE_LABELS[col]}={value:.2f}" for col, value in top.items()),
                **{FEATURE_LABELS[col]: float(row[col]) for col in FEATURE_COLS},
            }
        )
    centroid_df = pd.DataFrame(centroid_records)

    comp_rows = []
    for database in DB_ORDER:
        sub = assignment[assignment["database"] == database]
        total = len(sub)
        for cluster_name in ordered_cluster_names:
            count = int((sub["cluster"] == cluster_name).sum())
            comp_rows.append(
                {
                    "mode": mode,
                    "database": database,
                    "cluster": cluster_name,
                    "cluster_description": descriptions[cluster_name],
                    "count": count,
                    "pct": count / total * 100 if total else 0.0,
                }
            )
    composition = pd.DataFrame(comp_rows)

    return {
        "mode": mode,
        "x": x,
        "k": chosen_k,
        "k_selection": k_selection,
        "selection_reason": selection_reason,
        "assignment": assignment,
        "centroids": centroids,
        "centroid_df": centroid_df,
        "composition": composition,
        "descriptions": descriptions,
        "cluster_sizes": size_map,
        "pca_explained": explained,
        "ordered_cluster_names": ordered_cluster_names,
    }


def cluster_color_map(cluster_names: list[str]) -> dict[str, str]:
    return {cluster: CLUSTER_PALETTE[i % len(CLUSTER_PALETTE)] for i, cluster in enumerate(cluster_names)}


def add_cov_ellipse(ax: plt.Axes, x: np.ndarray, y: np.ndarray, color: str) -> None:
    if len(x) < 4:
        return
    cov = np.cov(x, y)
    vals, vecs = np.linalg.eigh(cov)
    if np.any(vals <= 0):
        return
    order = vals.argsort()[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    width, height = 2.0 * np.sqrt(vals) * 2.0
    ellipse = Ellipse(
        xy=(float(np.mean(x)), float(np.mean(y))),
        width=float(width),
        height=float(height),
        angle=float(angle),
        facecolor=color,
        edgecolor=color,
        alpha=0.10,
        linewidth=1.0,
        zorder=1,
    )
    ax.add_patch(ellipse)


def plot_pca_panel(ax: plt.Axes, assignment: pd.DataFrame, explained: np.ndarray, title: str) -> None:
    rng = np.random.default_rng(21)
    max_points = {"BiGG": 73, "VMH": 360, "ModelSEED": 79, "CarveMe": 520}
    for database in DB_ORDER:
        full = assignment[assignment["database"] == database]
        plot_df = full
        if len(full) > max_points[database]:
            plot_df = full.iloc[rng.choice(len(full), size=max_points[database], replace=False)]
        ax.scatter(
            plot_df["PC1"],
            plot_df["PC2"],
            s=15,
            alpha=0.30,
            c=DB_COLORS[database],
            label=f"{database} (n={len(full):,})",
            linewidths=0,
            zorder=2,
        )
        add_cov_ellipse(ax, full["PC1"].to_numpy(), full["PC2"].to_numpy(), DB_COLORS[database])
        centroid = full[["PC1", "PC2"]].mean()
        ax.scatter(
            [centroid["PC1"]],
            [centroid["PC2"]],
            s=95,
            c=DB_COLORS[database],
            edgecolor="white",
            linewidth=1.3,
            zorder=5,
        )
    ax.set_title(title, loc="left", fontsize=11.2, fontweight="bold", pad=8)
    ax.set_xlabel(f"PC1 ({explained[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({explained[1] * 100:.1f}%)")
    ax.grid(True, linestyle="--", alpha=0.22)
    ax.legend(frameon=False, fontsize=7.3, loc="upper right")


def plot_heatmap_panel(
    fig: plt.Figure,
    ax: plt.Axes,
    centroids: pd.DataFrame,
    centroid_df: pd.DataFrame,
    title: str,
) -> None:
    heat = centroids[FEATURE_COLS].to_numpy(dtype=float)
    vmax = float(np.nanpercentile(np.abs(heat), 95))
    vmax = max(vmax, 0.5)
    cmap = LinearSegmentedColormap.from_list("mqc_balance", ["#355C7D", "#F7F5EF", "#B94748"])
    im = ax.imshow(heat, aspect="auto", cmap=cmap, vmin=-vmax, vmax=vmax)
    ylabels = [
        f"{row.cluster}  {shorten_label(row.cluster_description)}  (n={int(row.cluster_size)})"
        for row in centroid_df.itertuples(index=False)
    ]
    ax.set_yticks(np.arange(len(centroids)))
    ax.set_yticklabels(ylabels, fontsize=7.1)
    ax.set_xticks(np.arange(len(FEATURE_COLS)))
    ax.set_xticklabels([FEATURE_LABELS[col] for col in FEATURE_COLS], rotation=52, ha="right", fontsize=6.7)
    ax.set_title(title, loc="left", fontsize=11.2, fontweight="bold", pad=24)
    ax.tick_params(length=0)
    for boundary in [len(COFACTOR_COLS) - 0.5, len(COFACTOR_COLS) + len(MATERIAL_COLS) - 0.5]:
        ax.axvline(boundary, color="white", linewidth=1.5)
    group_centers = [
        (len(COFACTOR_COLS) / 2 - 0.5, "Cofactor infeasibility"),
        (len(COFACTOR_COLS) + len(MATERIAL_COLS) / 2 - 0.5, "Material"),
        (len(COFACTOR_COLS) + len(MATERIAL_COLS) + len(BIOMASS_COLS) / 2 - 0.5, "Biomass/growth"),
    ]
    for xpos, label in group_centers:
        ax.text(xpos, 1.035, label, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6.9, color="#58606A")
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, ax=ax, fraction=0.028, pad=0.012)
    cbar.set_label("Standardized intensity", fontsize=7.2)
    cbar.ax.tick_params(labelsize=7)


def plot_composition_panel(
    ax: plt.Axes,
    composition: pd.DataFrame,
    cluster_names: list[str],
    colors: dict[str, str],
    title: str,
    show_legend: bool = True,
) -> None:
    left = np.zeros(len(DB_ORDER))
    y_positions = np.arange(len(DB_ORDER))
    for cluster in cluster_names:
        vals = (
            composition[composition["cluster"] == cluster]
            .set_index("database")
            .reindex(DB_ORDER)["pct"]
            .fillna(0.0)
            .to_numpy()
        )
        desc = composition[composition["cluster"] == cluster]["cluster_description"].iloc[0]
        ax.barh(
            y_positions,
            vals,
            left=left,
            height=0.58,
            color=colors[cluster],
            edgecolor="white",
            linewidth=0.8,
            label=f"{cluster}  {desc}",
        )
        left += vals
    ax.set_yticks(y_positions)
    ax.set_yticklabels(DB_ORDER)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Cluster composition within database (%)")
    ax.set_title(title, loc="left", fontsize=11.2, fontweight="bold", pad=8)
    ax.grid(True, axis="x", linestyle="--", alpha=0.22)
    ax.set_axisbelow(True)
    ax.invert_yaxis()
    if show_legend:
        ax.legend(frameon=False, fontsize=7.3, ncol=1, loc="lower right", bbox_to_anchor=(1.02, -0.03))


def shorten_label(label: str, max_chars: int = 34) -> str:
    replacements = {
        "Biomass-defect enriched: ": "Bio-defect: ",
        "Cofactor-dominant: ": "Cofactor: ",
    }
    for old, new in replacements.items():
        label = label.replace(old, new)
    return label if len(label) <= max_chars else f"{label[: max_chars - 1]}..."


def make_reproduction_figure(results: dict[str, object]) -> None:
    assignment = results["assignment"]
    explained = results["pca_explained"]
    centroids = results["centroids"]
    centroid_df = results["centroid_df"]
    composition = results["composition"]
    k_selection = results["k_selection"]
    cluster_names = results["ordered_cluster_names"]
    colors = cluster_color_map(cluster_names)

    fig = plt.figure(figsize=(14.2, 9.2), dpi=300)
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.0, 1.12],
        height_ratios=[1.0, 0.92],
        left=0.075,
        right=0.982,
        top=0.88,
        bottom=0.155,
        wspace=0.40,
        hspace=0.55,
    )
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    plot_pca_panel(ax1, assignment, explained, "a  PCA of model-level mechanistic error profiles")
    plot_heatmap_panel(fig, ax2, centroids, centroid_df, f"b  Cluster centroids, original scaling (k = {results['k']})")
    plot_composition_panel(ax3, composition, cluster_names, colors, "c  Cluster mixture across GEM resources", show_legend=False)
    handles, labels = ax3.get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize=6.8, ncol=3, loc="lower center", bbox_to_anchor=(0.43, 0.04))

    ax4.plot(k_selection["k"], k_selection["silhouette"], marker="o", color="#2F6FA3", linewidth=2.0)
    chosen = k_selection[k_selection["chosen"]].iloc[0]
    ax4.scatter([chosen["k"]], [chosen["silhouette"]], s=75, color="#D9892B", zorder=5)
    ax4.set_title("d  Cluster-number selection", loc="left", fontsize=12.5, fontweight="bold")
    ax4.set_xlabel("Number of clusters (k)")
    ax4.set_ylabel("Silhouette score")
    ax4.set_xticks(k_selection["k"])
    ax4.grid(True, linestyle="--", alpha=0.22)

    fig.suptitle("Reproduction of mechanistic clustering from MQC *_new.xlsx tables", x=0.06, y=0.985, ha="left", fontsize=17, fontweight="bold")
    fig.text(0.06, 0.952, "Original z-score scaling is retained here as a traceable baseline for comparison.", fontsize=9.5, color="#5F6B7A")
    save_figure(fig, "MQC图2_mechanistic_reproduction")


def make_refined_figure(results: dict[str, object], prevalence: pd.DataFrame) -> None:
    assignment = results["assignment"]
    explained = results["pca_explained"]
    centroids = results["centroids"]
    centroid_df = results["centroid_df"]
    composition = results["composition"]
    cluster_names = results["ordered_cluster_names"]
    colors = cluster_color_map(cluster_names)

    fig = plt.figure(figsize=(15.8, 10.2), dpi=300)
    gs = fig.add_gridspec(
        2,
        3,
        width_ratios=[0.76, 1.16, 1.10],
        height_ratios=[0.88, 1.12],
        left=0.075,
        right=0.982,
        top=0.905,
        bottom=0.165,
        wspace=0.48,
        hspace=0.60,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1:])
    ax_c = fig.add_subplot(gs[1, :2])
    ax_d = fig.add_subplot(gs[1, 2])

    plot_prevalence_panel(ax_a, prevalence)
    plot_pca_panel(ax_b, assignment, explained, "b  Source-specific mechanistic error space")
    plot_heatmap_panel(fig, ax_c, centroids, centroid_df, f"c  Robust cluster signatures (k = {results['k']})")
    plot_composition_panel(ax_d, composition, cluster_names, colors, "d  Cluster composition by GEM resource", show_legend=False)

    handles, labels = ax_d.get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize=7.1, ncol=2, loc="lower center", bbox_to_anchor=(0.61, 0.055))

    fig.suptitle("Mechanistic error fingerprints captured by MQC", x=0.06, y=0.987, ha="left", fontsize=18, fontweight="bold")
    fig.text(
        0.06,
        0.955,
        "MQC errors are widespread, form reproducible mechanism-level combinations, and remain source-specific across GEM resources.",
        fontsize=9.8,
        color="#5F6B7A",
    )
    fig.text(
        0.06,
        0.025,
        "Continuous features were winsorized and robust-scaled; binary defect indicators were mean-centered. ModelSEED rows include many missing model names, so row_id is retained as the stable identifier.",
        fontsize=7.8,
        color="#7A8591",
    )
    save_figure(fig, "MQC图2_mechanistic_refined")


def plot_prevalence_panel(ax: plt.Axes, prevalence: pd.DataFrame) -> None:
    prevalence = prevalence.set_index("database").reindex(DB_ORDER).reset_index()
    x = np.arange(len(prevalence))
    bars = ax.bar(
        x,
        prevalence["error_rate_pct"],
        color=[DB_COLORS[db] for db in prevalence["database"]],
        edgecolor="white",
        linewidth=1.0,
        width=0.68,
    )
    for bar, row in zip(bars, prevalence.itertuples(index=False)):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            row.error_rate_pct + 2.0,
            f"{row.error_rate_pct:.1f}%",
            ha="center",
            va="bottom",
            fontsize=8,
            fontweight="bold",
        )
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            3,
            f"n={row.total_models:,}",
            ha="center",
            va="bottom",
            fontsize=7.2,
            color="white",
        )
    ax.set_ylim(0, 108)
    ax.set_xticks(x)
    ax.set_xticklabels(prevalence["database"], rotation=35, ha="right")
    ax.set_ylabel("Models failing MQC (%)")
    ax.set_title("a  MQC error prevalence", loc="left", fontsize=11.2, fontweight="bold", pad=8)
    ax.grid(True, axis="y", linestyle="--", alpha=0.20)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def save_figure(fig: plt.Figure, stem: str) -> None:
    for suffix in [".png", ".pdf", ".svg"]:
        fig.savefig(PAPER_DIR / f"{stem}{suffix}", dpi=300 if suffix == ".png" else None, bbox_inches="tight")
    plt.close(fig)


def write_outputs(
    feature_df: pd.DataFrame,
    refined_results: dict[str, object],
    reproduction_results: dict[str, object],
    refined_scaling: pd.DataFrame,
    z_scaling: pd.DataFrame,
    prevalence: pd.DataFrame,
) -> None:
    feature_df.to_csv(PAPER_DIR / "fig2_mechanistic_feature_matrix.csv", index=False)
    refined_results["assignment"].to_csv(PAPER_DIR / "fig2_mechanistic_cluster_assignments.csv", index=False)
    refined_results["centroid_df"].to_csv(PAPER_DIR / "fig2_mechanistic_cluster_centroids.csv", index=False)
    refined_results["composition"].to_csv(PAPER_DIR / "fig2_mechanistic_database_cluster_composition.csv", index=False)
    k_selection = pd.concat([reproduction_results["k_selection"], refined_results["k_selection"]], ignore_index=True)
    k_selection.to_csv(PAPER_DIR / "fig2_mechanistic_k_selection.csv", index=False)
    scaling = pd.concat(
        [
            z_scaling.assign(mode="reproduction"),
            refined_scaling.assign(mode="refined"),
        ],
        ignore_index=True,
        sort=False,
    )
    scaling.to_csv(PAPER_DIR / "fig2_mechanistic_feature_scaling.csv", index=False)
    reproduction_results["assignment"].to_csv(PAPER_DIR / "fig2_mechanistic_reproduction_cluster_assignments.csv", index=False)
    reproduction_results["centroid_df"].to_csv(PAPER_DIR / "fig2_mechanistic_reproduction_cluster_centroids.csv", index=False)
    prevalence.to_csv(PAPER_DIR / "fig2_mechanistic_error_prevalence.csv", index=False)

    replace_workbook_sheets(
        {
            "图2_new_特征矩阵": feature_df,
            "图2_new_聚类分配": refined_results["assignment"],
            "图2_new_聚类中心": refined_results["centroid_df"],
            "图2_new_数据库构成": refined_results["composition"],
            "图2_new_k选择": k_selection,
            "图2_new_特征缩放": scaling,
        }
    )
    update_note(refined_results)
    shutil.copy2(ORGANIZED_XLSX, HOME_XLSX)
    shutil.copy2(ORGANIZED_NOTE, HOME_NOTE)


def replace_workbook_sheets(sheets: dict[str, pd.DataFrame]) -> None:
    with pd.ExcelWriter(ORGANIZED_XLSX, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
        for sheet_name, dataframe in sheets.items():
            dataframe.to_excel(writer, sheet_name=sheet_name, index=False)

    workbook = load_workbook(ORGANIZED_XLSX)
    for sheet_name in sheets:
        ws = workbook[sheet_name]
        ws.freeze_panes = "A2"
        ws.auto_filter.ref = ws.dimensions
        for column_cells in ws.columns:
            values = ["" if cell.value is None else str(cell.value) for cell in column_cells]
            width = min(max(len(value) for value in values) + 2, 55)
            ws.column_dimensions[column_cells[0].column_letter].width = width
    workbook.save(ORGANIZED_XLSX)


def update_note(refined_results: dict[str, object]) -> None:
    start = "<!-- FIG2_MECHANISTIC_START -->"
    end = "<!-- FIG2_MECHANISTIC_END -->"
    if ORGANIZED_NOTE.exists():
        text = ORGANIZED_NOTE.read_text(encoding="utf-8")
    else:
        text = "# MQC图2-5原始数据整理说明\n"
    if start in text and end in text:
        before = text.split(start)[0].rstrip()
        after = text.split(end, 1)[1].lstrip()
        text = before + "\n\n" + after
    section = f"""{start}
## 图2机制聚类补充

- 新增脚本：`{PAPER_DIR / 'build_mqc_fig2_mechanistic_atlas.py'}`
- 复现图：`{PAPER_DIR / 'MQC图2_mechanistic_reproduction.png'}`
- refined 主图：`{PAPER_DIR / 'MQC图2_mechanistic_refined.png'}`
- 输入表：`{INPUT_FILES['BiGG']}`、`{INPUT_FILES['VMH']}`、`{INPUT_FILES['ModelSEED']}`、`{INPUT_FILES['CarveMe']}`
- refined 聚类口径：连续特征 winsorize 后 robust scaling，二元缺陷特征 mean-centering；候选 k=4..8，要求最小簇规模不低于 `max(20, 0.5% total)`，并在 silhouette 差距 0.01 内优先较小 k。
- 本次 refined 选择：`k={refined_results['k']}`；新增 sheet 包括 `图2_new_特征矩阵`、`图2_new_聚类分配`、`图2_new_聚类中心`、`图2_new_数据库构成`、`图2_new_k选择`、`图2_new_特征缩放`。
{end}
"""
    ORGANIZED_NOTE.write_text(text.rstrip() + "\n\n" + section + "\n", encoding="utf-8")


def validate_outputs(
    feature_df: pd.DataFrame,
    refined_results: dict[str, object],
    refined_x: np.ndarray,
    refined_scaling: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    counts = feature_df["database"].value_counts().to_dict()
    for database, expected in EXPECTED_COUNTS.items():
        actual = int(counts.get(database, 0))
        rows.append({"check": f"input_count_{database}", "passed": actual == expected, "detail": f"{actual} / {expected}"})
    total_expected = sum(EXPECTED_COUNTS.values())
    rows.append({"check": "input_total", "passed": len(feature_df) == total_expected, "detail": f"{len(feature_df)} / {total_expected}"})
    rows.append({"check": "feature_no_nan", "passed": not feature_df[FEATURE_COLS].isna().any().any(), "detail": ""})
    rows.append({"check": "refined_matrix_finite", "passed": bool(np.isfinite(refined_x).all()), "detail": str(refined_x.shape)})
    binary_ok = all(set(feature_df[col].dropna().unique()).issubset({0, 1}) for col in BINARY_COLS)
    rows.append({"check": "binary_features_are_0_1", "passed": binary_ok, "detail": ", ".join(sorted(BINARY_COLS))})
    assignment = refined_results["assignment"]
    rows.append({"check": "one_cluster_per_model", "passed": len(assignment) == len(feature_df) and assignment["cluster"].notna().all(), "detail": str(len(assignment))})
    min_cluster = int(assignment["cluster"].value_counts().min())
    threshold = max(20, math.ceil(0.005 * len(feature_df)))
    rows.append({"check": "refined_min_cluster_size", "passed": min_cluster >= threshold, "detail": f"{min_cluster} >= {threshold}"})
    comp_sum = refined_results["composition"].groupby("database")["pct"].sum().round(6)
    rows.append({"check": "composition_sums_to_100", "passed": bool(np.allclose(comp_sum.to_numpy(), 100.0)), "detail": comp_sum.to_dict()})
    rows.append({"check": "scaling_stats_written", "passed": len(refined_scaling) == len(FEATURE_COLS), "detail": str(len(refined_scaling))})
    for stem in ["MQC图2_mechanistic_refined", "MQC图2_mechanistic_reproduction"]:
        for suffix in [".png", ".pdf", ".svg"]:
            path = PAPER_DIR / f"{stem}{suffix}"
            rows.append({"check": f"exists_{path.name}", "passed": path.exists() and path.stat().st_size > 0, "detail": str(path.stat().st_size if path.exists() else 0)})
    validation = pd.DataFrame(rows)
    validation.to_csv(PAPER_DIR / "fig2_mechanistic_validation_summary.csv", index=False)
    return validation


def main() -> None:
    feature_df = build_feature_matrix()
    prevalence = load_error_prevalence()

    raw_values = feature_df[FEATURE_COLS].to_numpy(dtype=float)
    z_x, z_scaling = standardize_matrix(raw_values)
    min_cluster_threshold = max(20, math.ceil(0.005 * len(feature_df)))
    refined_x, refined_scaling, _ = refined_scale(feature_df)

    reproduction_results = build_cluster_outputs(feature_df, z_x, "reproduction", min_cluster_threshold=None)
    refined_results = build_cluster_outputs(feature_df, refined_x, "refined", min_cluster_threshold=min_cluster_threshold)

    make_reproduction_figure(reproduction_results)
    make_refined_figure(refined_results, prevalence)
    validation = validate_outputs(feature_df, refined_results, refined_x, refined_scaling)
    write_outputs(feature_df, refined_results, reproduction_results, refined_scaling, z_scaling, prevalence)

    print("refined_k:", refined_results["k"])
    print("reproduction_k:", reproduction_results["k"])
    print("validation_passed:", bool(validation["passed"].all()))
    for path in [
        PAPER_DIR / "MQC图2_mechanistic_refined.png",
        PAPER_DIR / "MQC图2_mechanistic_reproduction.png",
        PAPER_DIR / "fig2_mechanistic_cluster_assignments.csv",
        ORGANIZED_XLSX,
    ]:
        print(path)


if __name__ == "__main__":
    main()
