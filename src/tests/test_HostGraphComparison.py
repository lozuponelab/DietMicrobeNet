"""
Unit tests for graph_comparisons_edge.py

Run with:
    python -m unittest test_graph_comparisons_edge_unittest.py -v
  or:
    pytest test_graph_comparisons_edge_unittest.py -v

The script under test must be importable as 'graph_comparisons_edge'
(i.e. graph_comparisons_edge.py must be in the same directory or on PYTHONPATH).
"""
from __future__ import annotations

import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from dietmicrobehost import host_GraphComparison as gc


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _make_metadata_csv(tmp: Path, graph_paths: list, names: list,
                       groups: list) -> Path:
    meta_path = tmp / "metadata.csv"
    pd.DataFrame({
        "paths":  graph_paths,
        "names":  names,
        "groups": groups,
    }).to_csv(meta_path, index=False)
    return meta_path


# ===========================================================================
# Main test class
# ===========================================================================

class TestGraphComparisonsEdge(unittest.TestCase):

    # -----------------------------------------------------------------------
    # setUp / tearDown
    # -----------------------------------------------------------------------

    def setUp(self):
        """Create a temporary directory with 4 sample graph CSVs and a metadata CSV."""
        self.tmpdir = tempfile.mkdtemp()
        tmp = Path(self.tmpdir)

        # Four samples with slightly different KO sets so similarity < 1.0
        graph_data = {
            "graph1.csv": pd.DataFrame([
                {"compound1_origin": "diet",        "compound2_origin": "microbe",
                 "compound3_origin": "host",         "edge1_KOs": "['K00001', 'K00002']",
                 "edge2_KOs": "['K00010']"},
                {"compound1_origin": "diet",        "compound2_origin": "microbediet",
                 "compound3_origin": "hostdiet",     "edge1_KOs": "['K00003']",
                 "edge2_KOs": "['K00011']"},
                {"compound1_origin": "microbediet", "compound2_origin": "microbediet",
                 "compound3_origin": "all",          "edge1_KOs": "['K00004']",
                 "edge2_KOs": "['K00012']"},
                {"compound1_origin": "all",         "compound2_origin": "all",
                 "compound3_origin": "all",          "edge1_KOs": "['K00005']",
                 "edge2_KOs": "['K00013']"},
            ]),
            "graph2.csv": pd.DataFrame([
                {"compound1_origin": "diet",        "compound2_origin": "microbe",
                 "compound3_origin": "host",         "edge1_KOs": "['K00001', 'K00006']",
                 "edge2_KOs": "['K00010', 'K00014']"},
                {"compound1_origin": "diet",        "compound2_origin": "microbediet",
                 "compound3_origin": "hostdiet",     "edge1_KOs": "['K00007']",
                 "edge2_KOs": "['K00015']"},
                {"compound1_origin": "microbediet", "compound2_origin": "microbediet",
                 "compound3_origin": "all",          "edge1_KOs": "['K00004', 'K00008']",
                 "edge2_KOs": "['K00016']"},
                {"compound1_origin": "all",         "compound2_origin": "all",
                 "compound3_origin": "all",          "edge1_KOs": "['K00009']",
                 "edge2_KOs": "['K00017']"},
            ]),
            "graph3.csv": pd.DataFrame([
                {"compound1_origin": "diet",        "compound2_origin": "microbe",
                 "compound3_origin": "host",         "edge1_KOs": "['K00001']",
                 "edge2_KOs": "['K00010', 'K00018']"},
                {"compound1_origin": "diet",        "compound2_origin": "microbediet",
                 "compound3_origin": "hostdiet",     "edge1_KOs": "['K00003', 'K00019']",
                 "edge2_KOs": "['K00020']"},
                {"compound1_origin": "microbediet", "compound2_origin": "microbediet",
                 "compound3_origin": "all",          "edge1_KOs": "['K00021']",
                 "edge2_KOs": "['K00022']"},
                {"compound1_origin": "all",         "compound2_origin": "all",
                 "compound3_origin": "all",          "edge1_KOs": "['K00023']",
                 "edge2_KOs": "['K00024']"},
            ]),
            "graph4.csv": pd.DataFrame([
                {"compound1_origin": "diet",        "compound2_origin": "microbe",
                 "compound3_origin": "host",         "edge1_KOs": "['K00001', 'K00025']",
                 "edge2_KOs": "['K00026']"},
                {"compound1_origin": "diet",        "compound2_origin": "microbediet",
                 "compound3_origin": "hostdiet",     "edge1_KOs": "['K00027']",
                 "edge2_KOs": "['K00028']"},
                {"compound1_origin": "microbediet", "compound2_origin": "microbediet",
                 "compound3_origin": "all",          "edge1_KOs": "['K00029']",
                 "edge2_KOs": "['K00030']"},
                {"compound1_origin": "all",         "compound2_origin": "all",
                 "compound3_origin": "all",          "edge1_KOs": "['K00031']",
                 "edge2_KOs": "['K00032']"},
            ]),
        }

        self.graph_paths = []
        for filename, df in graph_data.items():
            path = tmp / filename
            df.to_csv(path, index=False)
            self.graph_paths.append(str(path))

        self.names  = ["sample1", "sample2", "sample3", "sample4"]
        self.groups = ["A", "A", "B", "B"]

        self.metadata_csv = _make_metadata_csv(
            tmp, self.graph_paths, self.names, self.groups
        )

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    # -----------------------------------------------------------------------
    # Convenience helpers
    # -----------------------------------------------------------------------

    def _load_pattern_dicts(self):
        paths, names = gc.csv_to_inputs(
            str(self.metadata_csv), paths_col="paths", names_col="names"
        )
        graph_dict    = gc.get_graphs(paths, names)
        pattern_dicts = gc.subset_graphs(graph_dict)
        return pattern_dicts

    def _metadata_df(self):
        return pd.read_csv(self.metadata_csv).set_index("names")

    # -----------------------------------------------------------------------
    # 1. Module-level constants  (all accessed via gc.)
    # -----------------------------------------------------------------------

    def test_patterns_total_count(self):
        self.assertEqual(len(gc.PATTERNS), 36)

    def test_patterns_keys_are_unique(self):
        keys = [key for _, _, _, key in gc.PATTERNS]
        self.assertEqual(len(keys), len(set(keys)))

    def test_focal_pattern_key_value(self):
        self.assertEqual(gc.FOCAL_PATTERN_KEY, "diet_microbe_host")

    def test_focal_pattern_in_patterns(self):
        all_keys = {key for _, _, _, key in gc.PATTERNS}
        self.assertIn(gc.FOCAL_PATTERN_KEY, all_keys)

    def test_aggregate_keys_count(self):
        self.assertEqual(len(gc.AGGREGATE_PATTERN_KEYS), 35)

    def test_aggregate_keys_exclude_focal(self):
        self.assertNotIn(gc.FOCAL_PATTERN_KEY, gc.AGGREGATE_PATTERN_KEYS)

    def test_all_patterns_have_labels(self):
        for _, _, _, key in gc.PATTERNS:
            self.assertIn(key, gc.PATTERN_LABELS,
                          msg=f"'{key}' missing from PATTERN_LABELS")

    def test_edge_cols_constant(self):
        self.assertEqual(gc.EDGE_COLS, ("edge1_KOs", "edge2_KOs"))

    # -----------------------------------------------------------------------
    # 2. csv_to_inputs
    # -----------------------------------------------------------------------

    def test_csv_to_inputs_returns_correct_lengths(self):
        paths, names = gc.csv_to_inputs(
            str(self.metadata_csv), paths_col="paths", names_col="names"
        )
        self.assertEqual(len(paths), 4)
        self.assertEqual(len(names), 4)

    def test_csv_to_inputs_first_name(self):
        _, names = gc.csv_to_inputs(
            str(self.metadata_csv), paths_col="paths", names_col="names"
        )
        self.assertEqual(names[0], "sample1")

    def test_csv_to_inputs_missing_column_raises(self):
        with self.assertRaises(ValueError):
            gc.csv_to_inputs(
                str(self.metadata_csv), paths_col="NONEXISTENT", names_col="names"
            )

    # -----------------------------------------------------------------------
    # 3. get_graphs
    # -----------------------------------------------------------------------

    def test_get_graphs_returns_all_samples(self):
        paths, names = gc.csv_to_inputs(
            str(self.metadata_csv), paths_col="paths", names_col="names"
        )
        graph_dict = gc.get_graphs(paths, names)
        for name in self.names:
            self.assertIn(name, graph_dict)

    def test_get_graphs_values_are_dataframes(self):
        paths, names = gc.csv_to_inputs(
            str(self.metadata_csv), paths_col="paths", names_col="names"
        )
        graph_dict = gc.get_graphs(paths, names)
        for name, df in graph_dict.items():
            self.assertIsInstance(df, pd.DataFrame,
                                  msg=f"'{name}' is not a DataFrame")

    def test_get_graphs_duplicate_name_raises(self):
        p = self.graph_paths[0]
        with self.assertRaises(ValueError):
            gc.get_graphs([p, p], ["dup", "dup"])

    def test_get_graphs_missing_file_raises(self):
        with self.assertRaises(FileNotFoundError):
            gc.get_graphs(["/no/such/file.csv"], ["s1"])

    # -----------------------------------------------------------------------
    # 4. subset_graphs
    # -----------------------------------------------------------------------

    def test_subset_graphs_returns_all_36_pattern_keys(self):
        paths, names = gc.csv_to_inputs(
            str(self.metadata_csv), paths_col="paths", names_col="names"
        )
        graph_dict    = gc.get_graphs(paths, names)
        pattern_dicts = gc.subset_graphs(graph_dict)
        expected_keys = {key for _, _, _, key in gc.PATTERNS}
        self.assertEqual(set(pattern_dicts.keys()), expected_keys)

    def test_subset_graphs_focal_pattern_has_correct_rows(self):
        paths, names = gc.csv_to_inputs(
            str(self.metadata_csv), paths_col="paths", names_col="names"
        )
        graph_dict    = gc.get_graphs(paths, names)
        pattern_dicts = gc.subset_graphs(graph_dict)
        focal_df = pattern_dicts[gc.FOCAL_PATTERN_KEY]["sample1"]
        self.assertTrue(
            (focal_df["compound1_origin"] == "diet").all() and
            (focal_df["compound2_origin"] == "microbe").all() and
            (focal_df["compound3_origin"] == "host").all()
        )

    def test_subset_graphs_all_samples_present_in_focal(self):
        pattern_dicts = self._load_pattern_dicts()
        self.assertEqual(
            set(pattern_dicts[gc.FOCAL_PATTERN_KEY].keys()),
            set(self.names)
        )

    def test_subset_graphs_missing_origin_column_raises(self):
        bad_df = pd.DataFrame([{"edge1_KOs": "['K1']", "edge2_KOs": "['K2']"}])
        with self.assertRaises(ValueError):
            gc.subset_graphs({"s1": bad_df})

    # -----------------------------------------------------------------------
    # 5. get_kos_per_edge
    # -----------------------------------------------------------------------

    def test_get_kos_per_edge_returns_list_per_sample(self):
        pattern_dicts = self._load_pattern_dicts()
        focal_kos = gc.get_kos_per_edge(
            pattern_dicts[gc.FOCAL_PATTERN_KEY], "edge1_KOs"
        )
        for name in self.names:
            self.assertIsInstance(focal_kos[name], list)

    def test_get_kos_per_edge_edge1_non_empty(self):
        pattern_dicts = self._load_pattern_dicts()
        focal_kos = gc.get_kos_per_edge(
            pattern_dicts[gc.FOCAL_PATTERN_KEY], "edge1_KOs"
        )
        self.assertTrue(len(focal_kos["sample1"]) > 0)

    def test_get_kos_per_edge_edge1_and_edge2_differ(self):
        pattern_dicts = self._load_pattern_dicts()
        kos_e1 = gc.get_kos_per_edge(
            pattern_dicts[gc.FOCAL_PATTERN_KEY], "edge1_KOs"
        )
        kos_e2 = gc.get_kos_per_edge(
            pattern_dicts[gc.FOCAL_PATTERN_KEY], "edge2_KOs"
        )
        self.assertNotEqual(kos_e1["sample1"], kos_e2["sample1"])

    def test_get_kos_per_edge_result_is_sorted(self):
        pattern_dicts = self._load_pattern_dicts()
        focal_kos = gc.get_kos_per_edge(
            pattern_dicts[gc.FOCAL_PATTERN_KEY], "edge1_KOs"
        )
        for name, ko_list in focal_kos.items():
            self.assertEqual(ko_list, sorted(ko_list),
                             msg=f"KOs for '{name}' are not sorted")

    def test_get_kos_per_edge_no_duplicates(self):
        df = pd.DataFrame([
            {"edge1_KOs": "['K00001', 'K00001']", "edge2_KOs": "[]"},
            {"edge1_KOs": "['K00001']",            "edge2_KOs": "[]"},
        ])
        result = gc.get_kos_per_edge({"s1": df}, "edge1_KOs")
        self.assertEqual(result["s1"].count("K00001"), 1)

    def test_get_kos_per_edge_missing_column_returns_empty(self):
        df = pd.DataFrame([{"edge1_KOs": "['K00001']"}])  # no edge2_KOs
        result = gc.get_kos_per_edge({"s1": df}, "edge2_KOs")
        self.assertEqual(result["s1"], [])

    # -----------------------------------------------------------------------
    # 6. aggregate_kos_across_patterns
    # -----------------------------------------------------------------------

    def test_aggregate_kos_covers_all_samples(self):
        pattern_dicts = self._load_pattern_dicts()
        agg = gc.aggregate_kos_across_patterns(
            pattern_dicts, gc.AGGREGATE_PATTERN_KEYS, "edge1_KOs"
        )
        self.assertEqual(set(agg.keys()), set(self.names))

    def test_aggregate_kos_non_empty(self):
        pattern_dicts = self._load_pattern_dicts()
        agg = gc.aggregate_kos_across_patterns(
            pattern_dicts, gc.AGGREGATE_PATTERN_KEYS, "edge1_KOs"
        )
        for name in self.names:
            self.assertTrue(len(agg[name]) > 0,
                            msg=f"Aggregate KOs for '{name}' are empty")

    def test_aggregate_kos_is_superset_of_single_pattern(self):
        pattern_dicts = self._load_pattern_dicts()
        agg = gc.aggregate_kos_across_patterns(
            pattern_dicts, gc.AGGREGATE_PATTERN_KEYS, "edge1_KOs"
        )
        single = gc.get_kos_per_edge(
            pattern_dicts[gc.AGGREGATE_PATTERN_KEYS[0]], "edge1_KOs"
        )
        for name in self.names:
            self.assertGreaterEqual(
                len(agg[name]), len(single[name]),
                msg=f"Aggregate KOs for '{name}' smaller than single-pattern KOs"
            )

    def test_aggregate_kos_result_sorted(self):
        pattern_dicts = self._load_pattern_dicts()
        agg = gc.aggregate_kos_across_patterns(
            pattern_dicts, gc.AGGREGATE_PATTERN_KEYS, "edge1_KOs"
        )
        for name, ko_list in agg.items():
            self.assertEqual(ko_list, sorted(ko_list),
                             msg=f"Aggregate KOs for '{name}' are not sorted")

    # -----------------------------------------------------------------------
    # 7. jaccard
    # -----------------------------------------------------------------------

    def test_jaccard_both_empty(self):
        self.assertEqual(gc.jaccard(set(), set()), 1.0)

    def test_jaccard_identical(self):
        self.assertEqual(gc.jaccard({"KO1", "KO2"}, {"KO1", "KO2"}), 1.0)

    def test_jaccard_disjoint(self):
        self.assertEqual(gc.jaccard({"KO1"}, {"KO2"}), 0.0)

    def test_jaccard_partial_overlap(self):
        self.assertAlmostEqual(gc.jaccard({"KO1", "KO2"}, {"KO2", "KO3"}), 1 / 3)

    def test_jaccard_symmetry(self):
        a, b = {"KO1", "KO2"}, {"KO2", "KO3", "KO4"}
        self.assertAlmostEqual(gc.jaccard(a, b), gc.jaccard(b, a))

    # -----------------------------------------------------------------------
    # 8. calculate_similarity_matrix
    # -----------------------------------------------------------------------

    def test_similarity_matrix_shape(self):
        kos = {"s1": ["K00001", "K00002"], "s2": ["K00002", "K00003"]}
        matrix, labels = gc.calculate_similarity_matrix(kos)
        self.assertEqual(matrix.shape, (2, 2))
        self.assertEqual(len(labels), 2)

    def test_similarity_matrix_diagonal_is_one(self):
        kos = {"s1": ["K00001"], "s2": ["K00002", "K00003"]}
        matrix, labels = gc.calculate_similarity_matrix(kos)
        for i in range(len(labels)):
            self.assertAlmostEqual(matrix[i, i], 1.0)

    def test_similarity_matrix_is_symmetric(self):
        kos = {"s1": ["K00001", "K00002"], "s2": ["K00002", "K00003"]}
        matrix, _ = gc.calculate_similarity_matrix(kos)
        np.testing.assert_array_almost_equal(matrix, matrix.T)

    def test_similarity_matrix_known_value(self):
        # intersection={K00002}, union={K00001,K00002,K00003} → 1/3
        kos = {"s1": ["K00001", "K00002"], "s2": ["K00002", "K00003"]}
        matrix, labels = gc.calculate_similarity_matrix(kos)
        i, j = labels.index("s1"), labels.index("s2")
        self.assertAlmostEqual(matrix[i, j], 1 / 3)

    # -----------------------------------------------------------------------
    # 9. cluster_matrix
    # -----------------------------------------------------------------------

    def test_cluster_matrix_shape_preserved(self):
        kos = {"s1": ["K00001", "K00002"], "s2": ["K00002", "K00003"]}
        matrix, labels = gc.calculate_similarity_matrix(kos)
        ordered_matrix, ordered_labels, Z = gc.cluster_matrix(matrix, labels)
        self.assertEqual(ordered_matrix.shape, (2, 2))
        self.assertEqual(len(ordered_labels), 2)

    def test_cluster_matrix_labels_are_permutation(self):
        kos = {"s1": ["K1"], "s2": ["K2"], "s3": ["K1", "K2"]}
        matrix, labels = gc.calculate_similarity_matrix(kos)
        _, ordered_labels, _ = gc.cluster_matrix(matrix, labels)
        self.assertEqual(sorted(ordered_labels), sorted(labels))

    def test_cluster_matrix_single_sample(self):
        matrix = np.array([[1.0]])
        ordered_matrix, ordered_labels, Z = gc.cluster_matrix(matrix, ["s1"])
        self.assertEqual(ordered_labels, ["s1"])
        self.assertEqual(Z.size, 0)

    # -----------------------------------------------------------------------
    # 10. _safe_literal_eval
    # -----------------------------------------------------------------------

    def test_safe_literal_eval_list_string(self):
        result = gc._safe_literal_eval("['K00001', 'K00002']")
        self.assertEqual(result, ["K00001", "K00002"])

    def test_safe_literal_eval_empty_list(self):
        self.assertEqual(gc._safe_literal_eval("[]"), [])

    def test_safe_literal_eval_none(self):
        self.assertEqual(gc._safe_literal_eval(None), [])

    def test_safe_literal_eval_nan(self):
        self.assertEqual(gc._safe_literal_eval(float("nan")), [])

    def test_safe_literal_eval_empty_string(self):
        self.assertEqual(gc._safe_literal_eval(""), [])

    def test_safe_literal_eval_native_list(self):
        self.assertEqual(gc._safe_literal_eval(["K00001"]), ["K00001"])

    # -----------------------------------------------------------------------
    # 11. fdr_correction
    # -----------------------------------------------------------------------

    def test_fdr_correction_empty(self):
        self.assertEqual(gc.fdr_correction([]), [])

    def test_fdr_correction_length_preserved(self):
        result = gc.fdr_correction([0.01, 0.05, 0.2])
        self.assertEqual(len(result), 3)

    def test_fdr_correction_values_in_range(self):
        for val in gc.fdr_correction([0.001, 0.01, 0.05, 0.5]):
            self.assertGreaterEqual(val, 0.0)
            self.assertLessEqual(val, 1.0)

    def test_fdr_correction_single_value(self):
        result = gc.fdr_correction([0.04])
        self.assertAlmostEqual(result[0], 0.04)

    # -----------------------------------------------------------------------
    # 12. stat_test
    # -----------------------------------------------------------------------

    def test_stat_test_returns_required_keys(self):
        pattern_kos = {
            "sample1": ["K00001", "K00002"],
            "sample2": ["K00001", "K00003"],
            "sample3": ["K00004", "K00005"],
            "sample4": ["K00004", "K00006"],
        }
        metadata = pd.DataFrame(
            {"groups": ["A", "A", "B", "B"]},
            index=["sample1", "sample2", "sample3", "sample4"],
        )
        result = gc.stat_test(
            pattern_kos, metadata=metadata, group_col="groups",
            permutations=99, seed=0
        )
        for key in ("pseudo_F", "p_value", "R2", "groups", "n_samples"):
            self.assertIn(key, result)

    def test_stat_test_p_value_in_range(self):
        pattern_kos = {
            "sample1": ["K00001", "K00002"],
            "sample2": ["K00001", "K00003"],
            "sample3": ["K00004", "K00005"],
            "sample4": ["K00004", "K00006"],
        }
        metadata = pd.DataFrame(
            {"groups": ["A", "A", "B", "B"]},
            index=["sample1", "sample2", "sample3", "sample4"],
        )
        result = gc.stat_test(
            pattern_kos, metadata=metadata, group_col="groups",
            permutations=99, seed=0
        )
        self.assertGreaterEqual(result["p_value"], 0.0)
        self.assertLessEqual(result["p_value"], 1.0)

    def test_stat_test_r2_in_range(self):
        pattern_kos = {
            "sample1": ["K00001", "K00002"],
            "sample2": ["K00001", "K00003"],
            "sample3": ["K00004", "K00005"],
            "sample4": ["K00004", "K00006"],
        }
        metadata = pd.DataFrame(
            {"groups": ["A", "A", "B", "B"]},
            index=["sample1", "sample2", "sample3", "sample4"],
        )
        result = gc.stat_test(
            pattern_kos, metadata=metadata, group_col="groups",
            permutations=99, seed=0
        )
        self.assertGreaterEqual(result["R2"], 0.0)
        self.assertLessEqual(result["R2"], 1.0)

    def test_stat_test_missing_group_column_raises(self):
        pattern_kos = {"s1": ["K1"], "s2": ["K2"]}
        metadata = pd.DataFrame({"groups": ["A", "B"]}, index=["s1", "s2"])
        with self.assertRaises(KeyError):
            gc.stat_test(pattern_kos, metadata=metadata, group_col="WRONG",
                         permutations=99, seed=0)

    def test_stat_test_single_sample_raises(self):
        pattern_kos = {"s1": ["K1"]}
        metadata = pd.DataFrame({"groups": ["A"]}, index=["s1"])
        with self.assertRaises(ValueError):
            gc.stat_test(pattern_kos, metadata=metadata, group_col="groups",
                         permutations=99, seed=0)

    def test_stat_test_reproducible_with_seed(self):
        pattern_kos = {
            "sample1": ["K00001", "K00002"],
            "sample2": ["K00001", "K00003"],
            "sample3": ["K00004", "K00005"],
            "sample4": ["K00004", "K00006"],
        }
        metadata = pd.DataFrame(
            {"groups": ["A", "A", "B", "B"]},
            index=["sample1", "sample2", "sample3", "sample4"],
        )
        r1 = gc.stat_test(pattern_kos, metadata=metadata, group_col="groups",
                          permutations=99, seed=42)
        r2 = gc.stat_test(pattern_kos, metadata=metadata, group_col="groups",
                          permutations=99, seed=42)
        self.assertAlmostEqual(r1["p_value"], r2["p_value"])

    # -----------------------------------------------------------------------
    # 13. plotting
    # -----------------------------------------------------------------------

    def test_plotting_creates_heatmap(self):
        pattern_kos = {"sample1": ["K00001"], "sample2": ["K00001", "K00002"]}
        outdir = Path(self.tmpdir) / "plots"
        gc.plotting(pattern_kos, title="diet → microbe → host | edge1",
                    file_stem="focal_edge1", output=str(outdir))
        self.assertTrue((outdir / "focal_edge1_Heatmap.png").exists())

    def test_plotting_creates_dendrogram(self):
        pattern_kos = {"sample1": ["K00001"], "sample2": ["K00001", "K00002"]}
        outdir = Path(self.tmpdir) / "plots"
        gc.plotting(pattern_kos, title="t", file_stem="focal_edge1",
                    output=str(outdir))
        self.assertTrue((outdir / "focal_edge1_Dendrogram.png").exists())

    def test_plotting_creates_similarity_csv(self):
        pattern_kos = {"sample1": ["K00001"], "sample2": ["K00001", "K00002"]}
        outdir = Path(self.tmpdir) / "plots"
        gc.plotting(pattern_kos, title="t", file_stem="focal_edge1",
                    output=str(outdir))
        self.assertTrue((outdir / "SimilarityMatrix_focal_edge1.csv").exists())

    def test_plotting_at_least_one_file_created(self):
        pattern_kos = {"sample1": ["K00001"], "sample2": ["K00002"]}
        outdir = Path(self.tmpdir) / "plots2"
        gc.plotting(pattern_kos, title="t", file_stem="stem", output=str(outdir))
        self.assertTrue(len(list(outdir.glob("*"))) >= 1)

    def test_plotting_empty_input_no_crash(self):
        outdir = Path(self.tmpdir) / "empty_plots"
        try:
            gc.plotting({}, title="empty", file_stem="empty", output=str(outdir))
        except Exception as exc:
            self.fail(f"plotting() raised {exc} on empty input")

    def test_plotting_no_dendrogram_for_single_sample(self):
        outdir = Path(self.tmpdir) / "single_plots"
        gc.plotting({"s1": ["K1"]}, title="t", file_stem="single",
                    output=str(outdir))
        self.assertFalse((outdir / "single_Dendrogram.png").exists())

    # -----------------------------------------------------------------------
    # 14. summary
    # -----------------------------------------------------------------------

    def test_summary_creates_file(self):
        pattern_kos = {"sample1": ["K00001"], "sample2": ["K00001", "K00002"]}
        outdir = Path(self.tmpdir) / "summary"
        gc.summary(pattern_kos, title="diet → microbe → host | edge1",
                   file_stem="focal_edge1", stat=False,
                   metadata=pd.DataFrame(), groups=[], output=str(outdir))
        self.assertEqual(len(list(outdir.glob("*Summary.txt"))), 1)

    def test_summary_contains_title(self):
        pattern_kos = {"sample1": ["K00001"], "sample2": ["K00001", "K00002"]}
        outdir = Path(self.tmpdir) / "summary2"
        gc.summary(pattern_kos, title="diet → microbe → host | edge1",
                   file_stem="focal_edge1", stat=False,
                   metadata=pd.DataFrame(), groups=[], output=str(outdir))
        text = (outdir / "focal_edge1_Summary.txt").read_text()
        self.assertIn("diet → microbe → host | edge1", text)

    def test_summary_intersection_count(self):
        pattern_kos = {"sample1": ["K00001"], "sample2": ["K00001", "K00002"]}
        outdir = Path(self.tmpdir) / "summary3"
        gc.summary(pattern_kos, title="t", file_stem="stem", stat=False,
                   metadata=pd.DataFrame(), groups=[], output=str(outdir))
        text = (outdir / "stem_Summary.txt").read_text()
        self.assertIn("Number of KOs shared across all samples: 1", text)

    def test_summary_no_permanova_when_stat_false(self):
        pattern_kos = {"sample1": ["K00001"], "sample2": ["K00002"]}
        outdir = Path(self.tmpdir) / "summary4"
        gc.summary(pattern_kos, title="t", file_stem="stem", stat=False,
                   metadata=pd.DataFrame(), groups=[], output=str(outdir))
        text = (outdir / "stem_Summary.txt").read_text()
        self.assertNotIn("PERMANOVA", text)

    def test_summary_permanova_present_when_stat_true(self):
        pattern_kos = {
            "sample1": ["K00001", "K00002"],
            "sample2": ["K00001", "K00003"],
            "sample3": ["K00004", "K00005"],
            "sample4": ["K00004", "K00006"],
        }
        metadata = pd.DataFrame(
            {"groups": ["A", "A", "B", "B"]},
            index=["sample1", "sample2", "sample3", "sample4"],
        )
        outdir = Path(self.tmpdir) / "summary5"
        gc.summary(pattern_kos, title="t", file_stem="stem", stat=True,
                   metadata=metadata, groups=["groups"], output=str(outdir))
        text = (outdir / "stem_Summary.txt").read_text()
        self.assertIn("PERMANOVA", text)

    def test_summary_fdr_section_present_when_stat_true(self):
        pattern_kos = {
            "sample1": ["K00001", "K00002"],
            "sample2": ["K00001", "K00003"],
            "sample3": ["K00004", "K00005"],
            "sample4": ["K00004", "K00006"],
        }
        metadata = pd.DataFrame(
            {"groups": ["A", "A", "B", "B"]},
            index=["sample1", "sample2", "sample3", "sample4"],
        )
        outdir = Path(self.tmpdir) / "summary6"
        gc.summary(pattern_kos, title="t", file_stem="stem", stat=True,
                   metadata=metadata, groups=["groups"], output=str(outdir))
        text = (outdir / "stem_Summary.txt").read_text()
        self.assertIn("FDR", text)

    # -----------------------------------------------------------------------
    # 15. Integration: full pipeline for focal and aggregate
    # -----------------------------------------------------------------------

    def test_integration_focal_edge1_pipeline(self):
        pattern_dicts = self._load_pattern_dicts()
        focal_kos = gc.get_kos_per_edge(
            pattern_dicts[gc.FOCAL_PATTERN_KEY], "edge1_KOs"
        )
        outdir = Path(self.tmpdir) / "integ_focal_e1"
        gc.plotting(focal_kos, title="focal edge1", file_stem="focal_e1",
                    output=str(outdir))
        gc.summary(focal_kos, title="focal edge1", file_stem="focal_e1",
                   stat=False, metadata=pd.DataFrame(), groups=[],
                   output=str(outdir))
        self.assertTrue((outdir / "focal_e1_Heatmap.png").exists())
        self.assertTrue((outdir / "focal_e1_Summary.txt").exists())

    def test_integration_focal_edge2_pipeline(self):
        pattern_dicts = self._load_pattern_dicts()
        focal_kos = gc.get_kos_per_edge(
            pattern_dicts[gc.FOCAL_PATTERN_KEY], "edge2_KOs"
        )
        outdir = Path(self.tmpdir) / "integ_focal_e2"
        gc.plotting(focal_kos, title="focal edge2", file_stem="focal_e2",
                    output=str(outdir))
        gc.summary(focal_kos, title="focal edge2", file_stem="focal_e2",
                   stat=False, metadata=pd.DataFrame(), groups=[],
                   output=str(outdir))
        self.assertTrue((outdir / "focal_e2_Heatmap.png").exists())
        self.assertTrue((outdir / "focal_e2_Summary.txt").exists())

    def test_integration_aggregate_edge1_pipeline(self):
        pattern_dicts = self._load_pattern_dicts()
        agg_kos = gc.aggregate_kos_across_patterns(
            pattern_dicts, gc.AGGREGATE_PATTERN_KEYS, "edge1_KOs"
        )
        outdir = Path(self.tmpdir) / "integ_agg_e1"
        gc.plotting(agg_kos, title="agg edge1", file_stem="agg_e1",
                    output=str(outdir))
        gc.summary(agg_kos, title="agg edge1", file_stem="agg_e1",
                   stat=False, metadata=pd.DataFrame(), groups=[],
                   output=str(outdir))
        self.assertTrue((outdir / "agg_e1_Heatmap.png").exists())
        self.assertTrue((outdir / "agg_e1_Summary.txt").exists())

    def test_integration_aggregate_edge2_pipeline(self):
        pattern_dicts = self._load_pattern_dicts()
        agg_kos = gc.aggregate_kos_across_patterns(
            pattern_dicts, gc.AGGREGATE_PATTERN_KEYS, "edge2_KOs"
        )
        outdir = Path(self.tmpdir) / "integ_agg_e2"
        gc.plotting(agg_kos, title="agg edge2", file_stem="agg_e2",
                    output=str(outdir))
        gc.summary(agg_kos, title="agg edge2", file_stem="agg_e2",
                   stat=False, metadata=pd.DataFrame(), groups=[],
                   output=str(outdir))
        self.assertTrue((outdir / "agg_e2_Heatmap.png").exists())
        self.assertTrue((outdir / "agg_e2_Summary.txt").exists())

    def test_integration_focal_and_aggregate_kos_differ(self):
        pattern_dicts = self._load_pattern_dicts()
        focal_kos = gc.get_kos_per_edge(
            pattern_dicts[gc.FOCAL_PATTERN_KEY], "edge1_KOs"
        )
        agg_kos = gc.aggregate_kos_across_patterns(
            pattern_dicts, gc.AGGREGATE_PATTERN_KEYS, "edge1_KOs"
        )
        for name in self.names:
            self.assertGreaterEqual(len(agg_kos[name]), len(focal_kos[name]))

    def test_integration_similarity_matrix_from_real_data(self):
        pattern_dicts = self._load_pattern_dicts()
        focal_kos = gc.get_kos_per_edge(
            pattern_dicts[gc.FOCAL_PATTERN_KEY], "edge1_KOs"
        )
        matrix, labels = gc.calculate_similarity_matrix(focal_kos)
        self.assertEqual(matrix.shape, (len(self.names), len(self.names)))
        for i in range(len(labels)):
            self.assertAlmostEqual(matrix[i, i], 1.0)


if __name__ == "__main__":
    unittest.main()