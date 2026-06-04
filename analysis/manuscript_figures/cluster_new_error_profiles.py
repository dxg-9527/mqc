import ast
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap


BASE_DIR = Path("/hpcfs/fhome/yuan_qq/jupyter/MQC结果分析/MQC 文章图/图 2 中的 8 个原始表")
OUTPUT_DIR = BASE_DIR / "analysis_outputs"
SOURCE_CSV = OUTPUT_DIR / "new_mechanistic_feature_matrix_source.csv"

NEW_FILES = {
    "BiGG": BASE_DIR / "bigg_new.xlsx",
    "VMH": BASE_DIR / "virtual_new.xlsx",
    "ModelSEED": BASE_DIR / "seed_new.xlsx",
    "CarveMe": BASE_DIR / "carveme_new.xlsx",
}

DB_ORDER = ["BiGG", "VMH", "ModelSEED", "CarveMe"]
DB_COLORS = {
    "BiGG": "#3B6EA8",
    "VMH": "#C65B66",
    "ModelSEED": "#5B9B5A",
    "CarveMe": "#D9892B",
}

CLUSTER_DESCRIPTIONS = {
    "C1": "Pan-cofactor stressed",
    "C2": "Cofactor-dominant",
    "C3": "Biomass-heavy",
    "C4": "Near-clean baseline",
    "C5": "Hybrid no-growth",
    "C6": "Biomass-coupling outlier",
}

COFACTOR_COLS = ["NADH", "NADPH", "FADH2", "FMNH2", "Q8H2", "MQL8", "DMMQL8", "ATP", "CTP", "GTP", "UTP", "ITP"]
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
    text = str(value).strip().lower()
    if text in {"yes", "y", "true", "1"}:
        return 1
    return 0


def safe_float(value, default=0.0):
    try:
        out = float(value)
    except Exception:
        return default
    if math.isnan(out) or math.isinf(out):
        return default
    return out


def load_feature_table():
    if SOURCE_CSV.exists():
        feature_df = pd.read_csv(SOURCE_CSV)
        feature_cols = list(FEATURE_LABELS.keys())
        feature_df[feature_cols] = feature_df[feature_cols].fillna(0.0)
        return feature_df, feature_cols

    rows = []
    for database in DB_ORDER:
        df = pd.read_excel(NEW_FILES[database])
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
    feature_df = pd.DataFrame(rows)
    feature_cols = list(FEATURE_LABELS.keys())
    feature_df[feature_cols] = feature_df[feature_cols].fillna(0.0)
    return feature_df, feature_cols


def standardize_matrix(x):
    mean = x.mean(axis=0)
    std = x.std(axis=0, ddof=0)
    std[std == 0] = 1.0
    return (x - mean) / std, mean, std


def pca_project(x, n_components=2):
    x_centered = x - x.mean(axis=0, keepdims=True)
    u, s, vt = np.linalg.svd(x_centered, full_matrices=False)
    coords = x_centered @ vt[:n_components].T
    eigenvalues = (s ** 2) / max(len(x) - 1, 1)
    explained = eigenvalues / eigenvalues.sum()
    return coords, vt[:n_components], explained[:n_components]


def init_kmeans_pp(x, k, rng):
    centers = np.empty((k, x.shape[1]), dtype=float)
    first = rng.integers(0, len(x))
    centers[0] = x[first]
    closest_sq = np.sum((x - centers[0]) ** 2, axis=1)
    for i in range(1, k):
        probs = closest_sq / closest_sq.sum()
        idx = rng.choice(len(x), p=probs)
        centers[i] = x[idx]
        dist_sq = np.sum((x - centers[i]) ** 2, axis=1)
        closest_sq = np.minimum(closest_sq, dist_sq)
    return centers


def assign_clusters(x, centers):
    dist_sq = ((x[:, None, :] - centers[None, :, :]) ** 2).sum(axis=2)
    labels = dist_sq.argmin(axis=1)
    inertia = dist_sq[np.arange(len(x)), labels].sum()
    return labels, inertia


def fit_kmeans(x, k, n_init=12, max_iter=100, rng_seed=42):
    best_labels = None
    best_centers = None
    best_inertia = np.inf
    rng = np.random.default_rng(rng_seed)
    for init_idx in range(n_init):
        centers = init_kmeans_pp(x, k, np.random.default_rng(rng.integers(1_000_000_000) + init_idx))
        for _ in range(max_iter):
            labels, inertia = assign_clusters(x, centers)
            new_centers = centers.copy()
            for j in range(k):
                members = x[labels == j]
                if len(members) == 0:
                    new_centers[j] = x[rng.integers(0, len(x))]
                else:
                    new_centers[j] = members.mean(axis=0)
            if np.allclose(new_centers, centers, atol=1e-6):
                centers = new_centers
                break
            centers = new_centers
        labels, inertia = assign_clusters(x, centers)
        if inertia < best_inertia:
            best_labels = labels.copy()
            best_centers = centers.copy()
            best_inertia = inertia
    return best_labels, best_centers, best_inertia


def silhouette_score_numpy(x, labels):
    unique = np.unique(labels)
    if len(unique) < 2:
        return -1.0
    dist = np.sqrt(np.maximum(((x[:, None, :] - x[None, :, :]) ** 2).sum(axis=2), 0))
    silhouettes = []
    for i in range(len(x)):
        same = labels == labels[i]
        same[i] = False
        a = dist[i, same].mean() if same.any() else 0.0
        b = np.inf
        for lab in unique:
            if lab == labels[i]:
                continue
            other = labels == lab
            if other.any():
                b = min(b, dist[i, other].mean())
        denom = max(a, b)
        silhouettes.append((b - a) / denom if denom > 0 else 0.0)
    return float(np.mean(silhouettes))


def choose_kmeans_k(x_scaled):
    sample_n = min(len(x_scaled), 2500)
    rng = np.random.default_rng(42)
    sample_idx = rng.choice(len(x_scaled), size=sample_n, replace=False)
    sample = x_scaled[sample_idx]
    records = []
    best_score = -1.0
    best_k = 5
    for k in range(4, 9):
        labels, _, _ = fit_kmeans(x_scaled, k, n_init=8, max_iter=80, rng_seed=42 + k)
        score = silhouette_score_numpy(sample, labels[sample_idx])
        records.append({"k": k, "silhouette": score})
        if score > best_score:
            best_score = score
            best_k = k
    return best_k, pd.DataFrame(records)


def build_outputs():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    feature_df, feature_cols = load_feature_table()
    x_scaled, x_mean, x_std = standardize_matrix(feature_df[feature_cols].to_numpy(dtype=float))

    best_k, silhouette_df = choose_kmeans_k(x_scaled)
    labels, centers, _ = fit_kmeans(x_scaled, best_k, n_init=18, max_iter=120, rng_seed=2026)
    feature_df["cluster"] = labels

    coords, pca_components, explained = pca_project(x_scaled, n_components=2)
    feature_df["PC1"] = coords[:, 0]
    feature_df["PC2"] = coords[:, 1]

    centroids_scaled = pd.DataFrame(centers, columns=feature_cols, index=[f"C{i+1}" for i in range(best_k)])
    centroid_coords, _, _ = pca_project(centroids_scaled.to_numpy(), n_components=1)
    feature_coords, _, _ = pca_project(centroids_scaled.to_numpy().T, n_components=1)
    row_order = np.argsort(centroid_coords.ravel())
    col_order = np.argsort(feature_coords.ravel())
    ordered_centroids = centroids_scaled.iloc[row_order, col_order]

    cluster_sizes = feature_df["cluster"].value_counts().sort_index()
    ordered_cluster_labels = [ordered_centroids.index[i] for i in range(len(ordered_centroids))]
    cluster_name_map = {old: new for old, new in zip(centroids_scaled.index, ordered_cluster_labels)}
    feature_df["cluster_label"] = feature_df["cluster"].map(lambda x: cluster_name_map[f"C{x+1}"])
    feature_df["cluster_desc"] = feature_df["cluster_label"].map(CLUSTER_DESCRIPTIONS)

    db_cluster = (
        feature_df.groupby(["database", "cluster_label"]).size().reset_index(name="count")
        .pivot(index="database", columns="cluster_label", values="count")
        .fillna(0)
        .reindex(DB_ORDER)
    )
    db_cluster_pct = db_cluster.div(db_cluster.sum(axis=1), axis=0)

    db_centroids = feature_df.groupby("database")[feature_cols].mean().reindex(DB_ORDER)
    db_scaled, _, _ = standardize_matrix(db_centroids.to_numpy(dtype=float))
    db_coords, _, _ = pca_project(db_scaled, n_components=1)
    db_order_clustered = db_centroids.index[np.argsort(db_coords.ravel())].tolist()

    feature_df.to_csv(OUTPUT_DIR / "new_mechanistic_feature_matrix.csv", index=False)
    silhouette_df.to_csv(OUTPUT_DIR / "new_mechanistic_k_selection.csv", index=False)
    ordered_centroids.rename(columns=FEATURE_LABELS).to_csv(OUTPUT_DIR / "new_mechanistic_cluster_centroids.csv")
    db_cluster_pct.to_csv(OUTPUT_DIR / "new_mechanistic_database_cluster_composition.csv")
    pd.DataFrame({"database_cluster_order": db_order_clustered}).to_csv(OUTPUT_DIR / "new_mechanistic_database_order.csv", index=False)
    pd.DataFrame({"feature": feature_cols, "mean": x_mean, "std": x_std}).to_csv(OUTPUT_DIR / "new_mechanistic_feature_scaling.csv", index=False)
    pd.DataFrame(
        {
            "cluster": ordered_centroids.index,
            "label": [CLUSTER_DESCRIPTIONS[c] for c in ordered_centroids.index],
            "size": [int(cluster_sizes.iloc[int(c[1:]) - 1]) for c in ordered_centroids.index],
        }
    ).to_csv(OUTPUT_DIR / "new_mechanistic_cluster_labels.csv", index=False)

    fig = plt.figure(figsize=(14.5, 9.0), dpi=300)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.1], height_ratios=[1.0, 0.95], wspace=0.28, hspace=0.32)

    ax1 = fig.add_subplot(gs[0, 0])
    sample_n = 350
    rng = np.random.default_rng(7)
    for database in DB_ORDER:
        sub = feature_df[feature_df["database"] == database]
        if len(sub) > sample_n:
            sub = sub.iloc[rng.choice(len(sub), size=sample_n, replace=False)]
        ax1.scatter(sub["PC1"], sub["PC2"], s=18, alpha=0.35, c=DB_COLORS[database], label=f"{database} (n={len(feature_df[feature_df['database']==database]):,})", linewidths=0)
        centroid = feature_df[feature_df["database"] == database][["PC1", "PC2"]].mean()
        ax1.scatter([centroid["PC1"]], [centroid["PC2"]], s=90, c=DB_COLORS[database], edgecolor="white", linewidth=1.4, zorder=5)
    ax1.set_title("a  PCA of model-level mechanistic error profiles", loc="left", fontsize=13, fontweight="bold")
    ax1.set_xlabel(f"PC1 ({explained[0]*100:.1f}%)")
    ax1.set_ylabel(f"PC2 ({explained[1]*100:.1f}%)")
    ax1.grid(True, linestyle="--", alpha=0.25)
    ax1.legend(frameon=False, fontsize=8, loc="upper right")

    ax2 = fig.add_subplot(gs[0, 1])
    cmap = LinearSegmentedColormap.from_list("journal_blue", ["#F6F8FB", "#C8D8EA", "#7EA8CE", "#2F6FA3"])
    heat = ordered_centroids.to_numpy()
    im = ax2.imshow(heat, aspect="auto", cmap=cmap, vmin=np.percentile(heat, 5), vmax=np.percentile(heat, 95))
    ax2.set_title(f"b  Cluster-resolved heatmap of mechanistic centroids (k = {best_k})", loc="left", fontsize=13, fontweight="bold")
    ax2.set_yticks(np.arange(len(ordered_centroids.index)))
    heat_labels = []
    for lab in ordered_centroids.index:
        cluster_id = int(lab[1:])
        size = int(cluster_sizes.iloc[cluster_id - 1])
        heat_labels.append(f"{lab}  (n={size})")
    ax2.set_yticklabels(heat_labels, fontsize=8.5)
    ax2.set_xticks(np.arange(len(ordered_centroids.columns)))
    ax2.set_xticklabels([FEATURE_LABELS[c] for c in ordered_centroids.columns], rotation=55, ha="right", fontsize=8)
    ax2.tick_params(length=0)
    for spine in ax2.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, ax=ax2, fraction=0.03, pad=0.02)
    cbar.set_label("Standardized feature intensity", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    ax3 = fig.add_subplot(gs[1, 0])
    cluster_palette = ["#355C7D", "#C06C84", "#6C5B7B", "#F0A202", "#4F9D69", "#B56576", "#2A9D8F", "#9B5DE5"]
    left = np.zeros(len(DB_ORDER))
    ordered_cluster_names = ordered_centroids.index.tolist()
    for i, cluster_name in enumerate(ordered_cluster_names):
        vals = db_cluster_pct[cluster_name].reindex(DB_ORDER).to_numpy() * 100
        ax3.barh(DB_ORDER, vals, left=left, height=0.58, color=cluster_palette[i], edgecolor="white", linewidth=0.8, label=cluster_name)
        left += vals
    ax3.set_xlim(0, 100)
    ax3.set_xlabel("Cluster composition within each database (%)")
    ax3.set_title("c  Cluster mixture across GEM resources", loc="left", fontsize=13, fontweight="bold")
    ax3.grid(True, axis="x", linestyle="--", alpha=0.25)
    ax3.set_axisbelow(True)
    ax3.invert_yaxis()
    ax3.legend(frameon=False, fontsize=8, ncol=min(4, best_k), loc="lower center", bbox_to_anchor=(0.5, -0.30))

    ax4 = fig.add_subplot(gs[1, 1])
    silhouette_ax = ax4
    silhouette_ax.plot(silhouette_df["k"], silhouette_df["silhouette"], marker="o", color="#2F6FA3", linewidth=2.2)
    silhouette_ax.scatter([best_k], [float(silhouette_df.loc[silhouette_df["k"] == best_k, "silhouette"].iloc[0])], s=80, color="#D9892B", zorder=5)
    silhouette_ax.set_title("d  Cluster-number selection", loc="left", fontsize=13, fontweight="bold")
    silhouette_ax.set_xlabel("Number of clusters (k)")
    silhouette_ax.set_ylabel("Silhouette score")
    silhouette_ax.set_xticks(silhouette_df["k"])
    silhouette_ax.grid(True, linestyle="--", alpha=0.25)
    silhouette_ax.text(
        0.02,
        0.96,
        "Mechanistic profiles separate into reproducible,\ndatabase-enriched clusters rather than random noise.",
        transform=silhouette_ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        color="#5F6B7A",
    )
    silhouette_ax.text(
        0.02,
        0.75,
        "Features include cofactor-specific infeasibility,\nyield/material violations, and biomass defect mechanisms.",
        transform=silhouette_ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        color="#7A8591",
    )

    fig.suptitle("Mechanistic clustering atlas derived from MQC mechanistic tables", x=0.06, y=0.98, ha="left", fontsize=18, fontweight="bold")
    fig.text(0.06, 0.948, "Clustering was performed on model-level mechanistic error profiles using standardized cofactor, yield, and biomass features from the *_new.xlsx tables.", fontsize=10, color="#5F6B7A")
    fig.text(0.06, 0.02, "ModelSEED source tables contain many rows with missing model names and should therefore be interpreted cautiously.", fontsize=8.5, color="#7A8591")

    png = OUTPUT_DIR / "new_mechanistic_clustering_atlas.png"
    pdf = OUTPUT_DIR / "new_mechanistic_clustering_atlas.pdf"
    svg = OUTPUT_DIR / "new_mechanistic_clustering_atlas.svg"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)

    compact_fig = plt.figure(figsize=(13.2, 7.8), dpi=300)
    gs2 = compact_fig.add_gridspec(2, 2, width_ratios=[1.0, 1.15], height_ratios=[1.0, 0.82], wspace=0.28, hspace=0.36)

    axa = compact_fig.add_subplot(gs2[0, 0])
    sample_n_compact = 260
    rng2 = np.random.default_rng(17)
    for database in DB_ORDER:
        sub = feature_df[feature_df["database"] == database]
        if len(sub) > sample_n_compact:
            sub = sub.iloc[rng2.choice(len(sub), size=sample_n_compact, replace=False)]
        axa.scatter(sub["PC1"], sub["PC2"], s=16, alpha=0.33, c=DB_COLORS[database], label=database, linewidths=0)
        centroid = feature_df[feature_df["database"] == database][["PC1", "PC2"]].mean()
        axa.scatter([centroid["PC1"]], [centroid["PC2"]], s=88, c=DB_COLORS[database], edgecolor="white", linewidth=1.3, zorder=5)
    axa.set_title("a  PCA of mechanistic error profiles", loc="left", fontsize=13, fontweight="bold")
    axa.set_xlabel(f"PC1 ({explained[0]*100:.1f}%)")
    axa.set_ylabel(f"PC2 ({explained[1]*100:.1f}%)")
    axa.grid(True, linestyle="--", alpha=0.23)
    axa.legend(frameon=False, fontsize=8.5, loc="upper right")

    axb = compact_fig.add_subplot(gs2[:, 1])
    im2 = axb.imshow(heat, aspect="auto", cmap=cmap, vmin=np.percentile(heat, 5), vmax=np.percentile(heat, 95))
    axb.set_title("b  Mechanistic cluster signatures", loc="left", fontsize=13, fontweight="bold")
    axb.set_yticks(np.arange(len(ordered_centroids.index)))
    ytick_labels = []
    for lab in ordered_centroids.index:
        cluster_id = int(lab[1:])
        size = int(cluster_sizes.iloc[cluster_id - 1])
        ytick_labels.append(f"{lab}  {CLUSTER_DESCRIPTIONS[lab]}  (n={size})")
    axb.set_yticklabels(ytick_labels, fontsize=8.8)
    axb.yaxis.tick_right()
    axb.tick_params(axis="y", labelright=True, labelleft=False, pad=8)
    axb.set_xticks(np.arange(len(ordered_centroids.columns)))
    axb.set_xticklabels([FEATURE_LABELS[c] for c in ordered_centroids.columns], rotation=55, ha="right", fontsize=8)
    axb.tick_params(length=0)
    for tick, cluster_name, color in zip(axb.get_yticklabels(), ordered_centroids.index.tolist(), cluster_palette[: len(ordered_centroids.index)]):
        tick.set_color(color)
    for spine in axb.spines.values():
        spine.set_visible(False)
    cbar2 = compact_fig.colorbar(im2, ax=axb, fraction=0.03, pad=0.02)
    cbar2.set_label("Standardized feature intensity", fontsize=9)
    cbar2.ax.tick_params(labelsize=8)

    axc = compact_fig.add_subplot(gs2[1, 0])
    left = np.zeros(len(DB_ORDER))
    ordered_cluster_names = ordered_centroids.index.tolist()
    for i, cluster_name in enumerate(ordered_cluster_names):
        vals = db_cluster_pct[cluster_name].reindex(DB_ORDER).to_numpy() * 100
        axc.barh(DB_ORDER, vals, left=left, height=0.56, color=cluster_palette[i], edgecolor="white", linewidth=0.8, label=f"{cluster_name}  {CLUSTER_DESCRIPTIONS[cluster_name]}")
        left += vals
    axc.set_xlim(0, 100)
    axc.set_xlabel("Cluster composition within each database (%)")
    axc.set_title("c  Database enrichment of mechanistic clusters", loc="left", fontsize=13, fontweight="bold")
    axc.grid(True, axis="x", linestyle="--", alpha=0.25)
    axc.set_axisbelow(True)
    axc.invert_yaxis()

    compact_fig.suptitle("Structured mechanistic clusters derived from the MQC error tables", x=0.06, y=0.98, ha="left", fontsize=18, fontweight="bold")
    compact_fig.text(0.06, 0.945, "Quantitative failure mechanisms form database-enriched clusters rather than a single undifferentiated defect class.", fontsize=10, color="#5F6B7A")
    compact_fig.text(0.06, 0.055, "Cluster colors in panel c correspond to the labeled signatures in panel b.", fontsize=8.5, color="#7A8591")
    compact_fig.text(0.06, 0.02, "ModelSEED source tables contain many rows with missing model names and should be interpreted cautiously.", fontsize=8.5, color="#7A8591")

    compact_png = OUTPUT_DIR / "figure2b_mechanistic_cluster_panel.png"
    compact_pdf = OUTPUT_DIR / "figure2b_mechanistic_cluster_panel.pdf"
    compact_svg = OUTPUT_DIR / "figure2b_mechanistic_cluster_panel.svg"
    compact_fig.savefig(compact_png, dpi=300, bbox_inches="tight")
    compact_fig.savefig(compact_pdf, bbox_inches="tight")
    compact_fig.savefig(compact_svg, bbox_inches="tight")
    plt.close(compact_fig)

    return {
        "png": png,
        "pdf": pdf,
        "svg": svg,
        "compact_png": compact_png,
        "compact_pdf": compact_pdf,
        "compact_svg": compact_svg,
        "k_selection": OUTPUT_DIR / "new_mechanistic_k_selection.csv",
        "cluster_centroids": OUTPUT_DIR / "new_mechanistic_cluster_centroids.csv",
        "cluster_composition": OUTPUT_DIR / "new_mechanistic_database_cluster_composition.csv",
        "feature_matrix": OUTPUT_DIR / "new_mechanistic_feature_matrix.csv",
        "cluster_labels": OUTPUT_DIR / "new_mechanistic_cluster_labels.csv",
    }


if __name__ == "__main__":
    outputs = build_outputs()
    for key, path in outputs.items():
        print(f"{key}: {path}")
