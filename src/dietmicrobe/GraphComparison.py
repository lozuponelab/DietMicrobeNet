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


def csv_to_inputs(metadata: str, paths_col: str, names_col: str) -> Tuple[List[str], List[str]]:
    """Read metadata CSV and return paths and names"""
    df = pd.read_csv(metadata)
    for col in (paths_col, names_col):
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in metadata CSV ({metadata}). Available: {list(df.columns)}")

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
        df = pd.read_csv(pth)
        graph_dict[name] = df

    return graph_dict


def subset_graphs(graph_dict: Dict[str, pd.DataFrame]) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
    """Split graph results into main pattern (diet->microbe) and the rest 

    Validates required columns exist; raises informative errors if not.
    """
    required_cols = {"compound1_origin", "compound2_origin"}
    for name, df in graph_dict.items():
        if not required_cols.issubset(df.columns):
            raise ValueError(f"Missing required columns in graph dataframe for sample '{name}'. Required: {required_cols}. Found: {set(df.columns)}")

    food_microbe = {}
    least_restrictive = {}

    for name, df in graph_dict.items():
        food_microbe[name] = df[(df["compound1_origin"] == "food") & (df["compound2_origin"] == "microbe")]
        least_restrictive[name] = df[(df["compound1_origin"] != "food") & (df["compound2_origin"] != "microbe")]
        
    return food_microbe, least_restrictive


def _safe_literal_eval(item):
    """Try to parse a string representation of a list; return empty list on failure."""
    if item is None or (isinstance(item, float) and np.isnan(item)):
        return []
    # If it's already a list/tuple, return as list
    if isinstance(item, (list, tuple, set)):
        return list(item)
    # Strings: try literal_eval, fallback to simple parsing
    s = str(item).strip()
    if s == "":
        return []
    try:
        parsed = ast.literal_eval(s)
        if isinstance(parsed, (list, tuple, set)):
            return list(parsed)
        # if a single value (str) returned, wrap in list
        return [parsed]
    except Exception:
        # fallback: remove outer brackets/quotes and split on commas
        trimmed = s.strip("[](){} ")
        if trimmed == "":
            return []
        return [x.strip().strip("'\"") for x in trimmed.split(",") if x.strip()]


def get_kos(graph_dict: Dict[str, pd.DataFrame], ko_column_name: str = "KOs") -> Dict[str, List[str]]:
    """Extract KOs for each sample from the graph dataframes.

    Ensures unique, sorted KOs for determinism.
    """
    kos_dict: Dict[str, List[str]] = {}
    for name, df in graph_dict.items():
        if ko_column_name not in df.columns:
            logging.warning(f"'{ko_column_name}' column not found in dataframe for '{name}'. Using empty list.")
            kos_dict[name] = []
            continue

        ko_series = df[ko_column_name].tolist()
        parsed_lists = []
        for entry in ko_series:
            parsed = _safe_literal_eval(entry)
            parsed_lists.extend(parsed)

        # unique and sorted
        ko_set_sorted = sorted(set(parsed_lists))
        kos_dict[name] = ko_set_sorted

    return kos_dict


def jaccard(a: set, b: set) -> float:
    """Jaccard similarity for sets. Returns 1.0 when both empty (by convention)."""
    union = a | b
    if not union:
        return 1.0
    return len(a & b) / len(union)


def calculate_similarity_matrix(pattern_kos: Dict[str, List[str]]) -> Tuple[np.ndarray, List[str]]:
    """Return Jaccard similarity matrix (square) and labels.

    Internally converts lists to sets once for speed.
    """
    labels = list(pattern_kos.keys())
    n = len(labels)
    sets = {name: set(kos) for name, kos in pattern_kos.items()}
    matrix = np.zeros((n, n), dtype=float)

    for i, name_i in enumerate(labels):
        si = sets[name_i]
        for j, name_j in enumerate(labels):
            sj = sets[name_j]
            matrix[i, j] = jaccard(si, sj)

    return matrix, labels


def cluster_matrix(matrix: np.ndarray, labels: List[str]) -> Tuple[np.ndarray, List[str], np.ndarray]:
    """Cluster using hierarchical average linkage. Returns ordered square matrix, labels, and linkage Z."""
    if matrix.size == 0:
        return matrix, labels, np.array([])

    # distance matrix (0 on diagonal)
    distance_matrix = 1.0 - matrix

    # For clustering, SciPy expects condensed form:
    if len(distance_matrix) == 1:
        # Single item: nothing to cluster
        Z = np.array([])
        leaf_order = np.array([0], dtype=int)
    else:
        condensed = squareform(distance_matrix, checks=True)
        Z = linkage(condensed, method="average")
        leaf_order = leaves_list(Z)

    clustered_matrix = distance_matrix.copy()  # reuse shape
    if clustered_matrix.size:
        clustered_matrix = matrix[np.ix_(leaf_order, leaf_order)]
    clustered_labels = [labels[i] for i in leaf_order]

    return clustered_matrix, clustered_labels, Z


def stat_test(pattern_dict: Dict[str, List[str]],
              metadata: pd.DataFrame,
              group_col: str,
              permutations: int = 5000,
              seed: int = 5):

    """
    Run PERMANOVA on the Jaccard distance matrix for the given pattern.

    Automatically:
    - Removes groups with < 2 samples
    - Skips test if fewer than 2 valid groups remain
    """

    matrix, labels = calculate_similarity_matrix(pattern_dict)

    if len(labels) < 2:
        raise ValueError("Need at least 2 samples to run PERMANOVA.")

    # ---- Validate grouping column ----
    if group_col not in metadata.columns:
        raise KeyError(
            f"Grouping column '{group_col}' not found in metadata. "
            f"Available columns: {list(metadata.columns)}"
        )

    # ---- Validate samples ----
    missing_samples = set(labels) - set(metadata.index)
    if missing_samples:
        raise KeyError(
            f"Some samples missing in metadata for PERMANOVA: {missing_samples}"
        )

    group_series = metadata.loc[labels, group_col]

    # ---- Automatic group size detection ----
    counts = group_series.value_counts()
    small_groups = counts[counts < 2]

    if not small_groups.empty:
        logging.warning(
            f"Dropping groups with <2 samples in '{group_col}': "
            f"{small_groups.to_dict()}"
        )

        valid_groups = counts[counts >= 2].index
        group_series = group_series[group_series.isin(valid_groups)]

        # subset matrix to remaining samples
        valid_labels = group_series.index.tolist()
        idx = [labels.index(lab) for lab in valid_labels]
        matrix = matrix[np.ix_(idx, idx)]
        labels = valid_labels

    # ---- Check if enough groups remain ----
    remaining_counts = group_series.value_counts()

    if len(remaining_counts) < 2:
        raise ValueError(
            f"After filtering small groups, fewer than 2 groups remain "
            f"in '{group_col}'. Cannot run PERMANOVA.\n"
            f"Remaining groups: {remaining_counts.to_dict()}"
        )

    dm = DistanceMatrix(1.0 - matrix, labels)

    res = permanova(
        distance_matrix=dm,
        grouping=group_series,
        permutations=permutations,
        seed=seed
    )

    # Extract values 
    pseudo_f = res.get("test statistic", np.nan)
    p_value = res.get("p-value", np.nan)

    # Compute R² manually
    # F = pseudo-F
    # k = number of groups
    # N = number of samples
    pseudo_f = res["test statistic"]
    n = len(group_series)
    k = group_series.nunique()
    r2 = (pseudo_f * (k - 1)) / (pseudo_f * (k - 1) + (n - k))

    return {
        "pseudo_F": pseudo_f,
        "p_value": p_value,
        "R2": r2,
        "groups": group_series.value_counts().to_dict(),
        "n_samples": len(group_series)
    }

def fdr_correction(pvalues: List[float]):
    """Benjamini-Hochberg FDR correction."""
    if len(pvalues) == 0:
        return []

    _, corrected, _, _ = multipletests(
        pvalues,
        alpha=0.05,
        method="fdr_bh"
    )

    return corrected

def plotting(pattern_dict: Dict[str, List[str]], pattern_name: str, output: str):
    """Create and save heatmap and dendrogram for the given pattern."""
    output_dir = Path(output)
    output_dir.mkdir(parents=True, exist_ok=True)

    matrix, labels = calculate_similarity_matrix(pattern_dict)

    if len(labels) == 0:
        logging.info(f"No samples for pattern '{pattern_name}'; skipping plots.")
        return

    ordered_matrix, ordered_labels, Z = cluster_matrix(matrix=matrix, labels=labels)

    n = len(ordered_labels)
    # Heatmap
    plt.figure(figsize=(max(6, n * 0.5), max(4, n * 0.5)))
    plt.imshow(ordered_matrix, cmap="viridis", vmin=0.0, vmax=1.0, aspect="auto")
    plt.xticks(range(n), ordered_labels, rotation=90)
    plt.yticks(range(n), ordered_labels)
    plt.colorbar(label="Jaccard Similarity")
    plt.title("Clustered Jaccard Similarity Heatmap: " + pattern_name)
    plt.tight_layout()
    heatmap_file = output_dir / f"{pattern_name.replace(' ', '')}_GraphComparisons_Heatmap.png"
    plt.savefig(heatmap_file)
    plt.close()
    logging.info(f"Saved heatmap to {heatmap_file}")

    # Dendrogram (only if more than 1 sample)
    if n > 1 and Z.size:
        plt.figure(figsize=(max(6, n * 0.2), 4))
        dendrogram(Z, labels=ordered_labels, leaf_rotation=90)
        plt.title("Hierarchical Clustering Dendrogram: " + pattern_name)
        plt.tight_layout()
        dend_file = output_dir / f"{pattern_name.replace(' ', '')}_GraphComparisons_Dendrogram.png"
        plt.savefig(dend_file)
        plt.close()
        logging.info(f"Saved dendrogram to {dend_file}")
    else:
        logging.info(f"Not enough samples to produce dendrogram for pattern '{pattern_name}' (n={n}).")

    # save similarity matrix 
    logging.info(f"Saved similarity matrix to {output_dir}/SimilarityMatrix_{pattern_name}.csv ")
    df = pd.DataFrame(matrix, columns=labels)
    df.to_csv(f"{output_dir}/SimilarityMatrix_{pattern_name}.csv", index=False)


def summary(pattern_dict: Dict[str, List[str]], pattern_name: str, stat: bool, metadata: pd.DataFrame, groups:list, output: str):
    """Write a summary file including intersection KOs, unique KOs per sample, and optional PERMANOVA results."""
    output_dir = Path(output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Intersection across all sets
    all_sets = [set(v) for v in pattern_dict.values()]
    if all_sets:
        intersection_set = set.intersection(*all_sets)
    else:
        intersection_set = set()

    unique_dict = {
        name: sorted(set(ko_list) - intersection_set)
        for name, ko_list in pattern_dict.items()
    }

    summary_file = output_dir / f"{pattern_name.replace(' ', '')}_GraphComparisons_Summary.txt"
    with open(summary_file, "w") as fh:
        fh.write(f"##### SUMMARY FOR PATTERN: {pattern_name} #####\n")
        fh.write(f"Number of KOs shared: {len(intersection_set)}\n")

        if stat:
            permanova_results = []
            permanova_groups = []

            for group in groups:
                try:
                    stat_results = stat_test(
                        pattern_dict=pattern_dict,
                        metadata=metadata,
                        group_col=group
                    )

                    permanova_results.append(stat_results)
                    permanova_groups.append(group)

                except Exception as e:
                    logging.warning(f"PERMANOVA failed for pattern '{pattern_name}', group '{group}': {e}")

            # ---- Write Raw Results ----
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

            # ---- FDR Correction ----
            if len(raw_pvalues) > 0:
                corrected_pvalues = fdr_correction(raw_pvalues)

                fh.write("\n##### FDR CORRECTED P-VALUES #####\n")

                for group, raw_p, adj_p in zip(permanova_groups, raw_pvalues, corrected_pvalues):
                    fh.write(
                        f"{group}: raw p = {raw_p:.6f}, "
                        f"FDR-adjusted p = {adj_p:.6f}\n"
                    )
        
        fh.write("\n##### UNIQUE KOs #####\n")
        for name, unique_list in unique_dict.items():
            fh.write(f"Unique KOs to {name}: {len(unique_list)}\n")

        fh.write("\n\n\n ###### LISTS OF KOs ######\n\n\n")
        fh.write(f"Intersection KOs: {sorted(intersection_set)}\n\n")

        for name, unique_list in unique_dict.items():
            fh.write(f"Unique KOs for {name}: {unique_list}\n\n")

    logging.info(f"Saved summary to {summary_file}")


def main():
    parser = argparse.ArgumentParser(description="Compare graph results across samples using KOs and Jaccard similarity.")
    parser.add_argument("-m", "--metadata", required=True, help="Metadata CSV containing file paths and names")
    parser.add_argument("-p", "--paths", required=True, help="Name of column containing file paths")
    parser.add_argument("-n", "--names", required=True, help="Name of column containing names of graphs (e.g., sampleID)")
    parser.add_argument("-s", "--stat_test", action="store_true", help="If statistical test for group comparison wanted include this parameter")
    parser.add_argument("-g", "--groups", help="Names of columns for use in PERMANOVA, if multiple separate by a comma e.g., cohort,diet,location", default="")
    parser.add_argument("-o", "--output", required=True, help="Output directory for plots and summary files")
    parser.add_argument("--ko_column", help="Name of KOs column in graph CSVs (default: 'KOs')", default="KOs")
    args = parser.parse_args()

    md = pd.read_csv(args.metadata)
    md = md.set_index(args.names) # must index by names 
    paths, names = csv_to_inputs(metadata=args.metadata, paths_col=args.paths, names_col=args.names)

    graphs_dict = get_graphs(paths=paths, names=names)
    food_microbe_dict, least_restrictive_dict = subset_graphs(graph_dict=graphs_dict)

    food_microbe_kos = get_kos(food_microbe_dict, ko_column_name=args.ko_column)
    rest_patterns_kos = get_kos(least_restrictive_dict, ko_column_name=args.ko_column)

    patterns = [food_microbe_kos, rest_patterns_kos]
    pattern_names = ["Food to Microbe", "Least Restrictive Patterns"]

    groups = [g.strip() for g in args.groups.split(',') if g.strip()]
    for pat_dict, pat_name in zip(patterns, pattern_names):
        plotting(pat_dict, pat_name, args.output)
        summary(pattern_dict=pat_dict, pattern_name=pat_name, stat=args.stat_test, metadata=md, groups=groups, output=args.output)


if __name__ == "__main__":
    main()