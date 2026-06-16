from __future__ import annotations
import argparse
import ast
import logging
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import squareform
from skbio import DistanceMatrix
from skbio.stats.distance import permanova
from statsmodels.stats.multitest import multipletests

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

PATTERNS: List[Tuple[str, str, str, str]] = [
    # --- diet as position 1 (12 patterns) ---
    ("diet",        "microbe",      "host",         "diet_microbe_host"),
    ("diet",        "microbe",      "hostdiet",     "diet_microbe_hostdiet"),
    ("diet",        "microbe",      "hostmicrobe",  "diet_microbe_hostmicrobe"),
    ("diet",        "microbe",      "all",          "diet_microbe_all"),
    ("diet",        "microbediet",  "host",         "diet_microbediet_host"),
    ("diet",        "microbediet",  "hostdiet",     "diet_microbediet_hostdiet"),
    ("diet",        "microbediet",  "hostmicrobe",  "diet_microbediet_hostmicrobe"),
    ("diet",        "microbediet",  "all",          "diet_microbediet_all"),
    ("diet",        "all",          "host",         "diet_all_host"),
    ("diet",        "all",          "hostdiet",     "diet_all_hostdiet"),
    ("diet",        "all",          "hostmicrobe",  "diet_all_hostmicrobe"),
    ("diet",        "all",          "all",          "diet_all_all"),
    # --- microbediet as position 1 (12 patterns) ---
    ("microbediet", "microbe",      "host",         "microbediet_microbe_host"),
    ("microbediet", "microbe",      "hostdiet",     "microbediet_microbe_hostdiet"),
    ("microbediet", "microbe",      "hostmicrobe",  "microbediet_microbe_hostmicrobe"),
    ("microbediet", "microbe",      "all",          "microbediet_microbe_all"),
    ("microbediet", "microbediet",  "host",         "microbediet_microbediet_host"),
    ("microbediet", "microbediet",  "hostdiet",     "microbediet_microbediet_hostdiet"),
    ("microbediet", "microbediet",  "hostmicrobe",  "microbediet_microbediet_hostmicrobe"),
    ("microbediet", "microbediet",  "all",          "microbediet_microbediet_all"),
    ("microbediet", "all",          "host",         "microbediet_all_host"),
    ("microbediet", "all",          "hostdiet",     "microbediet_all_hostdiet"),
    ("microbediet", "all",          "hostmicrobe",  "microbediet_all_hostmicrobe"),
    ("microbediet", "all",          "all",          "microbediet_all_all"),
    # --- all as position 1 (12 patterns) ---
    ("all",         "microbe",      "host",         "all_microbe_host"),
    ("all",         "microbe",      "hostdiet",     "all_microbe_hostdiet"),
    ("all",         "microbe",      "hostmicrobe",  "all_microbe_hostmicrobe"),
    ("all",         "microbe",      "all",          "all_microbe_all"),
    ("all",         "microbediet",  "host",         "all_microbediet_host"),
    ("all",         "microbediet",  "hostdiet",     "all_microbediet_hostdiet"),
    ("all",         "microbediet",  "hostmicrobe",  "all_microbediet_hostmicrobe"),
    ("all",         "microbediet",  "all",          "all_microbediet_all"),
    ("all",         "all",          "host",         "all_all_host"),
    ("all",         "all",          "hostdiet",     "all_all_hostdiet"),
    ("all",         "all",          "hostmicrobe",  "all_all_hostmicrobe"),
    ("all",         "all",          "all",          "all_all_all"),
]

# The single focal pattern processed individually per edge
FOCAL_PATTERN_KEY = "diet_microbe_host"

# The 35 remaining patterns whose KOs are pooled (per edge) into an aggregate
AGGREGATE_PATTERN_KEYS = [key for _, _, _, key in PATTERNS if key != FOCAL_PATTERN_KEY]

# Human-readable label for each pattern
PATTERN_LABELS: Dict[str, str] = {
    # diet as position 1
    "diet_microbe_host":                    "diet → microbe → host",
    "diet_microbe_hostdiet":                "diet → microbe → hostdiet",
    "diet_microbe_hostmicrobe":             "diet → microbe → hostmicrobe",
    "diet_microbe_all":                     "diet → microbe → all",
    "diet_microbediet_host":                "diet → microbediet → host",
    "diet_microbediet_hostdiet":            "diet → microbediet → hostdiet",
    "diet_microbediet_hostmicrobe":         "diet → microbediet → hostmicrobe",
    "diet_microbediet_all":                 "diet → microbediet → all",
    "diet_all_host":                        "diet → all → host",
    "diet_all_hostdiet":                    "diet → all → hostdiet",
    "diet_all_hostmicrobe":                 "diet → all → hostmicrobe",
    "diet_all_all":                         "diet → all → all",
    # microbediet as position 1
    "microbediet_microbe_host":             "microbediet → microbe → host",
    "microbediet_microbe_hostdiet":         "microbediet → microbe → hostdiet",
    "microbediet_microbe_hostmicrobe":      "microbediet → microbe → hostmicrobe",
    "microbediet_microbe_all":              "microbediet → microbe → all",
    "microbediet_microbediet_host":         "microbediet → microbediet → host",
    "microbediet_microbediet_hostdiet":     "microbediet → microbediet → hostdiet",
    "microbediet_microbediet_hostmicrobe":  "microbediet → microbediet → hostmicrobe",
    "microbediet_microbediet_all":          "microbediet → microbediet → all",
    "microbediet_all_host":                 "microbediet → all → host",
    "microbediet_all_hostdiet":             "microbediet → all → hostdiet",
    "microbediet_all_hostmicrobe":          "microbediet → all → hostmicrobe",
    "microbediet_all_all":                  "microbediet → all → all",
    # all as position 1
    "all_microbe_host":                     "all → microbe → host",
    "all_microbe_hostdiet":                 "all → microbe → hostdiet",
    "all_microbe_hostmicrobe":              "all → microbe → hostmicrobe",
    "all_microbe_all":                      "all → microbe → all",
    "all_microbediet_host":                 "all → microbediet → host",
    "all_microbediet_hostdiet":             "all → microbediet → hostdiet",
    "all_microbediet_hostmicrobe":          "all → microbediet → hostmicrobe",
    "all_microbediet_all":                  "all → microbediet → all",
    "all_all_host":                         "all → all → host",
    "all_all_hostdiet":                     "all → all → hostdiet",
    "all_all_hostmicrobe":                  "all → all → hostmicrobe",
    "all_all_all":                          "all → all → all",
}

# Edge column names present in the graph CSVs
EDGE_COLS = ("edge1_KOs", "edge2_KOs")
EDGE_LABELS = ("edge1", "edge2")


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def csv_to_inputs(metadata: str, paths_col: str, names_col: str) -> Tuple[List[str], List[str]]:
    """Read metadata CSV and return paths and names."""
    df = pd.read_csv(metadata)
    for col in (paths_col, names_col):
        if col not in df.columns:
            raise ValueError(
                f"Column '{col}' not found in metadata CSV ({metadata}). "
                f"Available: {list(df.columns)}"
            )
    paths_list = df[paths_col].astype(str).tolist()
    names_list = df[names_col].astype(str).tolist()
    if len(paths_list) != len(names_list):
        raise ValueError("Number of paths and names differ in metadata CSV.")
    return paths_list, names_list


def get_graphs(paths: List[str], names: List[str]) -> Dict[str, pd.DataFrame]:
    """Read each CSV path into a dataframe keyed by the corresponding name."""
    if len(paths) != len(names):
        raise ValueError("Number of names not equal to the number of paths provided.")
    graph_dict: Dict[str, pd.DataFrame] = {}
    for idx, p in enumerate(paths):
        name = names[idx]
        if name in graph_dict:
            raise ValueError(f"Duplicate name detected: '{name}'. Names must be unique.")
        pth = Path(p)
        if not pth.exists():
            raise FileNotFoundError(f"Graph CSV not found for sample '{name}': {p}")
        graph_dict[name] = pd.read_csv(pth)
    return graph_dict


def subset_graphs(
    graph_dict: Dict[str, pd.DataFrame]
) -> Dict[str, Dict[str, pd.DataFrame]]:
    """Split each sample's graph dataframe into the 12 pattern subsets."""
    required_cols = {"compound1_origin", "compound2_origin", "compound3_origin"}
    for name, df in graph_dict.items():
        missing = required_cols - set(df.columns)
        if missing:
            raise ValueError(
                f"Missing required columns in graph dataframe for sample '{name}': "
                f"{missing}. Found: {set(df.columns)}"
            )

    pattern_dicts: Dict[str, Dict[str, pd.DataFrame]] = {
        key: {} for _, _, _, key in PATTERNS
    }

    for name, df in graph_dict.items():
        for c1_origin, c2_origin, c3_origin, key in PATTERNS:
            mask = (
                (df["compound1_origin"] == c1_origin) &
                (df["compound2_origin"] == c2_origin) &
                (df["compound3_origin"] == c3_origin)
            )
            pattern_dicts[key][name] = df[mask]

    return pattern_dicts


# ---------------------------------------------------------------------------
# KO extraction — now per-edge
# ---------------------------------------------------------------------------

def _safe_literal_eval(item) -> List[str]:
    """Parse a string representation of a list; return empty list on failure."""
    if item is None or (isinstance(item, float) and np.isnan(item)):
        return []
    if isinstance(item, (list, tuple, set)):
        return list(item)
    s = str(item).strip()
    if s == "":
        return []
    try:
        parsed = ast.literal_eval(s)
        if isinstance(parsed, (list, tuple, set)):
            return list(parsed)
        return [parsed]
    except Exception:
        trimmed = s.strip("[](){} ")
        if trimmed == "":
            return []
        return [x.strip().strip("'\"") for x in trimmed.split(",") if x.strip()]


def get_kos_per_edge(
    graph_dict: Dict[str, pd.DataFrame],
    edge_col: str,
) -> Dict[str, List[str]]:
    """Extract KOs from a single edge column for each sample.

    Parameters
    ----------
    graph_dict:
        Mapping of sample name → (filtered) DataFrame for one pattern.
    edge_col:
        One of 'edge1_KOs' or 'edge2_KOs'.

    Returns
    -------
    Dict mapping sample name → sorted unique KO list.
    """
    kos_dict: Dict[str, List[str]] = {}
    for name, df in graph_dict.items():
        parsed: List[str] = []
        if edge_col in df.columns:
            for entry in df[edge_col].tolist():
                parsed.extend(_safe_literal_eval(entry))
        else:
            logging.warning(
                f"Edge column '{edge_col}' not found in dataframe for '{name}'. "
                f"Using empty list."
            )
        kos_dict[name] = sorted(set(parsed))
    return kos_dict


def aggregate_kos_across_patterns(
    pattern_dicts: Dict[str, Dict[str, pd.DataFrame]],
    pattern_keys: List[str],
    edge_col: str,
) -> Dict[str, List[str]]:
    """Pool KOs from *multiple* patterns into a single per-sample KO set.

    For each sample the union of KOs across all specified patterns (for the
    given edge column) is returned.

    Parameters
    ----------
    pattern_dicts:
        Full dict[pattern_key -> dict[sample_name -> DataFrame]].
    pattern_keys:
        Which patterns to pool together.
    edge_col:
        'edge1_KOs' or 'edge2_KOs'.

    Returns
    -------
    Dict mapping sample name → sorted unique KO list (union across patterns).
    """
    # Collect all sample names present in any of the requested patterns
    all_names: set = set()
    for key in pattern_keys:
        all_names.update(pattern_dicts[key].keys())

    pooled: Dict[str, set] = {name: set() for name in all_names}

    for key in pattern_keys:
        pat_kos = get_kos_per_edge(pattern_dicts[key], edge_col=edge_col)
        for name, ko_list in pat_kos.items():
            pooled[name].update(ko_list)

    return {name: sorted(kos) for name, kos in pooled.items()}


# ---------------------------------------------------------------------------
# Similarity / clustering
# ---------------------------------------------------------------------------

def jaccard(a: set, b: set) -> float:
    """Jaccard similarity. Returns 1.0 when both sets are empty."""
    union = a | b
    if not union:
        return 1.0
    return len(a & b) / len(union)


def calculate_similarity_matrix(
    pattern_kos: Dict[str, List[str]]
) -> Tuple[np.ndarray, List[str]]:
    """Return square Jaccard similarity matrix and ordered labels."""
    labels = list(pattern_kos.keys())
    n = len(labels)
    sets = {name: set(kos) for name, kos in pattern_kos.items()}
    matrix = np.zeros((n, n), dtype=float)
    for i, ni in enumerate(labels):
        for j, nj in enumerate(labels):
            matrix[i, j] = jaccard(sets[ni], sets[nj])
    return matrix, labels


def cluster_matrix(
    matrix: np.ndarray, labels: List[str]
) -> Tuple[np.ndarray, List[str], np.ndarray]:
    """Hierarchical average-linkage clustering."""
    if matrix.size == 0:
        return matrix, labels, np.array([])

    distance_matrix = 1.0 - matrix

    if len(distance_matrix) == 1:
        return matrix, labels, np.array([])

    condensed = squareform(distance_matrix, checks=True)
    Z = linkage(condensed, method="average")
    leaf_order = leaves_list(Z)
    clustered_matrix = matrix[np.ix_(leaf_order, leaf_order)]
    clustered_labels = [labels[i] for i in leaf_order]
    return clustered_matrix, clustered_labels, Z


# ---------------------------------------------------------------------------
# Statistical testing
# ---------------------------------------------------------------------------

def stat_test(
    pattern_kos: Dict[str, List[str]],
    metadata: pd.DataFrame,
    group_col: str,
    permutations: int = 5000,
    seed: int = 5,
) -> dict:
    """PERMANOVA on the Jaccard distance matrix."""
    matrix, labels = calculate_similarity_matrix(pattern_kos)

    if len(labels) < 2:
        raise ValueError("Need at least 2 samples to run PERMANOVA.")

    if group_col not in metadata.columns:
        raise KeyError(
            f"Grouping column '{group_col}' not found in metadata. "
            f"Available columns: {list(metadata.columns)}"
        )

    missing_samples = set(labels) - set(metadata.index)
    if missing_samples:
        raise KeyError(f"Some samples missing in metadata for PERMANOVA: {missing_samples}")

    group_series = metadata.loc[labels, group_col]

    counts = group_series.value_counts()
    small_groups = counts[counts < 2]
    if not small_groups.empty:
        logging.warning(
            f"Dropping groups with <2 samples in '{group_col}': {small_groups.to_dict()}"
        )
        valid_groups = counts[counts >= 2].index
        group_series = group_series[group_series.isin(valid_groups)]
        valid_labels = group_series.index.tolist()
        idx = [labels.index(lab) for lab in valid_labels]
        matrix = matrix[np.ix_(idx, idx)]
        labels = valid_labels

    remaining_counts = group_series.value_counts()
    if len(remaining_counts) < 2:
        raise ValueError(
            f"After filtering small groups, fewer than 2 groups remain in '{group_col}'. "
            f"Remaining: {remaining_counts.to_dict()}"
        )

    dm = DistanceMatrix(1.0 - matrix, labels)
    res = permanova(
        distance_matrix=dm,
        grouping=group_series,
        permutations=permutations,
        seed=seed,
    )

    pseudo_f = res["test statistic"]
    p_value = res.get("p-value", np.nan)
    n = len(group_series)
    k = group_series.nunique()
    r2 = (pseudo_f * (k - 1)) / (pseudo_f * (k - 1) + (n - k))

    return {
        "pseudo_F": pseudo_f,
        "p_value": p_value,
        "R2": r2,
        "groups": group_series.value_counts().to_dict(),
        "n_samples": len(group_series),
    }


def fdr_correction(pvalues: List[float]) -> List[float]:
    """Benjamini-Hochberg FDR correction."""
    if len(pvalues) == 0:
        return []
    _, corrected, _, _ = multipletests(pvalues, alpha=0.05, method="fdr_bh")
    return corrected


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plotting(
    pattern_kos: Dict[str, List[str]],
    title: str,
    file_stem: str,
    output: str,
) -> None:
    """Create and save heatmap, dendrogram, and similarity matrix CSV."""
    output_dir = Path(output)
    output_dir.mkdir(parents=True, exist_ok=True)

    matrix, labels = calculate_similarity_matrix(pattern_kos)

    if len(labels) == 0:
        logging.info(f"No samples for '{title}'; skipping plots.")
        return

    ordered_matrix, ordered_labels, Z = cluster_matrix(matrix=matrix, labels=labels)
    n = len(ordered_labels)

    # Heatmap
    plt.figure(figsize=(max(6, n * 0.5), max(4, n * 0.5)))
    plt.imshow(ordered_matrix, cmap="viridis", vmin=0.0, vmax=1.0, aspect="auto")
    plt.xticks(range(n), ordered_labels, rotation=90)
    plt.yticks(range(n), ordered_labels)
    plt.colorbar(label="Jaccard Similarity")
    plt.title(f"Clustered Jaccard Similarity Heatmap\n{title}")
    plt.tight_layout()
    heatmap_file = output_dir / f"{file_stem}_Heatmap.png"
    plt.savefig(heatmap_file)
    plt.close()
    logging.info(f"Saved heatmap to {heatmap_file}")

    # Dendrogram
    if n > 1 and Z.size:
        plt.figure(figsize=(max(6, n * 0.2), 4))
        dendrogram(Z, labels=ordered_labels, leaf_rotation=90)
        plt.title(f"Hierarchical Clustering Dendrogram\n{title}")
        plt.tight_layout()
        dend_file = output_dir / f"{file_stem}_Dendrogram.png"
        plt.savefig(dend_file)
        plt.close()
        logging.info(f"Saved dendrogram to {dend_file}")
    else:
        logging.info(f"Not enough samples for dendrogram '{title}' (n={n}).")

    # Similarity matrix CSV
    df = pd.DataFrame(matrix, columns=labels, index=labels)
    matrix_file = output_dir / f"SimilarityMatrix_{file_stem}.csv"
    df.to_csv(matrix_file)
    logging.info(f"Saved similarity matrix to {matrix_file}")


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def summary(
    pattern_kos: Dict[str, List[str]],
    title: str,
    file_stem: str,
    stat: bool,
    metadata: pd.DataFrame,
    groups: List[str],
    output: str,
) -> None:
    """Write a summary file with shared/unique KOs and optional PERMANOVA results."""
    output_dir = Path(output)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_sets = [set(v) for v in pattern_kos.values()]
    intersection_set = set.intersection(*all_sets) if all_sets else set()

    unique_dict = {
        name: sorted(set(ko_list) - intersection_set)
        for name, ko_list in pattern_kos.items()
    }

    summary_file = output_dir / f"{file_stem}_Summary.txt"
    with open(summary_file, "w") as fh:
        fh.write(f"##### SUMMARY: {title} #####\n")
        fh.write(f"Number of KOs shared across all samples: {len(intersection_set)}\n")

        if stat:
            permanova_results = []
            permanova_groups = []

            for group in groups:
                try:
                    res = stat_test(
                        pattern_kos=pattern_kos,
                        metadata=metadata,
                        group_col=group,
                    )
                    permanova_results.append(res)
                    permanova_groups.append(group)
                except Exception as e:
                    logging.warning(
                        f"PERMANOVA failed for '{title}', group '{group}': {e}"
                    )

            fh.write("\n##### PERMANOVA RESULTS #####\n")
            raw_pvalues = []
            for group, res in zip(permanova_groups, permanova_results):
                fh.write(f"\n--- Group: {group} ---\n")
                fh.write(f"Samples used: {res['n_samples']}\n")
                fh.write(f"Group sizes: {res['groups']}\n")
                fh.write(f"Pseudo-F: {res['pseudo_F']:.4f}\n")
                fh.write(f"R²: {res['R2']:.4f}\n")
                fh.write(f"P-value (raw): {res['p_value']:.6f}\n")
                raw_pvalues.append(res["p_value"])

            if raw_pvalues:
                corrected = fdr_correction(raw_pvalues)
                fh.write("\n##### FDR CORRECTED P-VALUES #####\n")
                for group, raw_p, adj_p in zip(permanova_groups, raw_pvalues, corrected):
                    fh.write(
                        f"{group}: raw p = {raw_p:.6f}, FDR-adjusted p = {adj_p:.6f}\n"
                    )

        fh.write("\n##### UNIQUE KOs PER SAMPLE #####\n")
        for name, unique_list in unique_dict.items():
            fh.write(f"Unique KOs to {name}: {len(unique_list)}\n")

        fh.write("\n\n\n###### LISTS OF KOs ######\n\n\n")
        fh.write(f"Intersection KOs: {sorted(intersection_set)}\n\n")
        for name, unique_list in unique_dict.items():
            fh.write(f"Unique KOs for {name}: {unique_list}\n\n")

    logging.info(f"Saved summary to {summary_file}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare graph results across samples using per-edge Jaccard similarity.\n"
            "\n"
            "Produces two sets of outputs, each split by edge (edge1 / edge2):\n"
            "  1. Focal pattern  : diet → microbe → host\n"
            "  2. Aggregate      : all 35 remaining patterns pooled per sample\n"
        )
    )
    parser.add_argument("-m", "--metadata", required=True,
                        help="Metadata CSV containing file paths and sample names")
    parser.add_argument("-p", "--paths", required=True,
                        help="Column name containing file paths")
    parser.add_argument("-n", "--names", required=True,
                        help="Column name containing sample names")
    parser.add_argument("-s", "--stat_test", action="store_true",
                        help="Run PERMANOVA statistical test for group comparison")
    parser.add_argument("-g", "--groups", default="",
                        help="Comma-separated metadata column names for PERMANOVA "
                             "(e.g. cohort,diet,location)")
    parser.add_argument("-o", "--output", required=True,
                        help="Output directory for plots and summary files")
    args = parser.parse_args()

    # ---- Load metadata ----
    md = pd.read_csv(args.metadata).set_index(args.names)

    # ---- Load and subset graphs ----
    paths, names = csv_to_inputs(
        metadata=args.metadata, paths_col=args.paths, names_col=args.names
    )
    graphs_dict = get_graphs(paths=paths, names=names)
    pattern_dicts = subset_graphs(graph_dict=graphs_dict)

    groups = [g.strip() for g in args.groups.split(",") if g.strip()]

    # ================================================================
    # BLOCK 1 — Focal pattern: diet → microbe → host, per edge
    # ================================================================
    focal_df_dict = pattern_dicts[FOCAL_PATTERN_KEY]   # dict[sample -> DataFrame]
    focal_label   = PATTERN_LABELS[FOCAL_PATTERN_KEY]  # "diet → microbe → host"

    for edge_col, edge_label in zip(EDGE_COLS, EDGE_LABELS):
        focal_kos = get_kos_per_edge(focal_df_dict, edge_col=edge_col)

        title     = f"{focal_label} | {edge_label}"
        file_stem = f"{FOCAL_PATTERN_KEY}_{edge_label}"

        logging.info(f"Processing: {title}")
        plotting(focal_kos, title=title, file_stem=file_stem, output=args.output)
        summary(
            pattern_kos=focal_kos,
            title=title,
            file_stem=file_stem,
            stat=args.stat_test,
            metadata=md,
            groups=groups,
            output=args.output,
        )

    # ================================================================
    # BLOCK 2 — Aggregate: 11 remaining patterns pooled, per edge
    # ================================================================
    agg_pattern_names = " | ".join(PATTERN_LABELS[k] for k in AGGREGATE_PATTERN_KEYS)
    logging.info(
        f"Aggregating {len(AGGREGATE_PATTERN_KEYS)} patterns: {AGGREGATE_PATTERN_KEYS}"
    )

    for edge_col, edge_label in zip(EDGE_COLS, EDGE_LABELS):
        agg_kos = aggregate_kos_across_patterns(
            pattern_dicts=pattern_dicts,
            pattern_keys=AGGREGATE_PATTERN_KEYS,
            edge_col=edge_col,
        )

        title     = f"Aggregated (35 patterns) | {edge_label}"
        file_stem = f"aggregated_35patterns_{edge_label}"

        logging.info(f"Processing: {title}")
        plotting(agg_kos, title=title, file_stem=file_stem, output=args.output)
        summary(
            pattern_kos=agg_kos,
            title=title,
            file_stem=file_stem,
            stat=args.stat_test,
            metadata=md,
            groups=groups,
            output=args.output,
        )


if __name__ == "__main__":
    main()