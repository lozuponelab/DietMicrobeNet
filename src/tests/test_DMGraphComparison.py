import unittest
import tempfile
import shutil
import os
from pathlib import Path
import pandas as pd
import numpy as np

# Import  module
from dietmicrobe import GraphComparison as gc


class TestGraphCompare(unittest.TestCase):

    def setUp(self):
        """Create a temporary directory and 4 example graph CSVs."""
        self.tmpdir = tempfile.mkdtemp()
        tmp = Path(self.tmpdir)

        # Define 4 example graphs
        graph_data = {
            "graph1.csv": pd.DataFrame({
                "compound1_origin": ["food", "food", "both", "food"],
                "compound2_origin": ["microbe", "both", "both", "microbe"],
                "KOs": ["['KO1', 'KO2']",
                        "['KO2']",
                        "['KO3']",
                        "['KO4']"]
            }),

            "graph2.csv": pd.DataFrame({
                "compound1_origin": ["food", "both", "both"],
                "compound2_origin": ["microbe", "both", "both"],
                "KOs": ["['KO1']",
                        "['KO3', 'KO5']",
                        "['KO6']"]
            }),

            "graph3.csv": pd.DataFrame({
                "compound1_origin": ["both", "both", "food"],
                "compound2_origin": ["both", "microbe", "both"],
                "KOs": ["['KO2', 'KO7']",
                        "['KO7']",
                        "['KO8']"]
            }),

            "graph4.csv": pd.DataFrame({
                "compound1_origin": ["food", "both", "both", "food"],
                "compound2_origin": ["microbe", "both", "both", "both"],
                "KOs": ["['KO1', 'KO9']",
                        "['KO10']",
                        "['KO3']",
                        "['KO4', 'KO11']"]
            }),
        }

        # Write 4 graph files
        graph_paths = []
        for filename, df in graph_data.items():
            path = tmp / filename
            df.to_csv(path, index=False)
            graph_paths.append(str(path))

        # Create metadata CSV for the 4 graphs
        self.metadata_csv = tmp / "metadata.csv"
        metadata_df = pd.DataFrame({
            "paths": graph_paths,           # matches csv_to_inputs signature
            "names": ["sample1", "sample2", "sample3", "sample4"],
            "groups": ["A", "A", "B", "B"]
        })
        metadata_df.to_csv(self.metadata_csv, index=False)


    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    # -------------------------------------------------------------------------
    # TESTS
    # -------------------------------------------------------------------------

    def test_csv_to_inputs(self):
        paths, names = gc.csv_to_inputs(
            metadata=str(self.metadata_csv),
            paths_col="paths",
            names_col="names"
        )
        self.assertEqual(len(paths), 4)
        self.assertEqual(names[0], "sample1")

    def test_get_graphs(self):
        paths, names = gc.csv_to_inputs(
            metadata=str(self.metadata_csv), paths_col="paths",
            names_col="names"
        )
        graph_dict = gc.get_graphs(paths, names)
        self.assertIn("sample1", graph_dict)
        self.assertIsInstance(graph_dict["sample1"], pd.DataFrame)

    def test_subset_graphs(self):
        paths, names = gc.csv_to_inputs(
            metadata=str(self.metadata_csv), paths_col="paths",
            names_col="names"
        )
        graph_dict = gc.get_graphs(paths, names)
        fm, rest = gc.subset_graphs(graph_dict)
        # all 4 graphs have at least one FM/rest entry
        self.assertTrue(len(fm) > 0)
        self.assertTrue(len(rest) > 0)

    def test_get_kos(self):
        paths, names = gc.csv_to_inputs(
            metadata=str(self.metadata_csv), paths_col="paths",
            names_col="names"
        )
        graph_dict = gc.get_graphs(paths, names)
        fm, _ = gc.subset_graphs(graph_dict)
        kos_dict = gc.get_kos(fm)

        first_sample = list(kos_dict.keys())[0]
        self.assertIsInstance(kos_dict[first_sample], list)
        self.assertTrue(len(kos_dict[first_sample]) > 0)

    def test_jaccard(self):
        self.assertEqual(gc.jaccard(set(), set()), 1.0)
        self.assertEqual(gc.jaccard({"a"}, {"a"}), 1.0)
        self.assertEqual(gc.jaccard({"a"}, {"b"}), 0.0)
        self.assertAlmostEqual(gc.jaccard({"a", "b"}, {"b", "c"}), 1/3)

    def test_similarity_matrix(self):
        pattern = {"s1": ["KO1", "KO2"], "s2": ["KO2"]}
        matrix, labels = gc.calculate_similarity_matrix(pattern)
        self.assertEqual(matrix.shape, (2, 2))
        self.assertEqual(labels, ["s1", "s2"])

    def test_cluster_matrix(self):
        pattern = {"s1": ["KO1", "KO2"], "s2": ["KO2"]}
        matrix, labels = gc.calculate_similarity_matrix(pattern)
        ordered_matrix, ordered_labels, Z = gc.cluster_matrix(matrix, labels)
        self.assertEqual(ordered_matrix.shape, (2, 2))
        self.assertEqual(len(ordered_labels), 2)

    def test_stat_test(self):
        import pandas as pd

        # Pattern dictionary
        pattern = {
            "sample1": ["KO1", "KO2"],
            "sample2": ["KO2", "KO3"],
            "sample3": ["KO1", "KO2"],  # ensure ≥2 samples per group
            "sample4": ["KO2", "KO3"]
        }

        # Groups dictionary
        groups = {
            "sample1": "A",
            "sample2": "A",
            "sample3": "B",
            "sample4": "B"
        }

        # Convert to metadata DataFrame
        metadata = pd.DataFrame({
            "group": groups
        })
        metadata.index.name = "sample"
        metadata.index = list(groups.keys())

        # Run the PERMANOVA test
        result = gc.stat_test(pattern, metadata=metadata, group_col="group")

        # Check output
        self.assertIn("pseudo_F", result)
        self.assertIn("p_value", result)
        self.assertIn("R2", result)


    def test_plotting(self):
        pattern = {
            "sample1": ["KO1"],
            "sample2": ["KO1", "KO2"]
        }
        outdir = Path(self.tmpdir) / "plots"
        gc.plotting(pattern, "Food to Microbe", str(outdir))

        files = list(outdir.glob("*"))
        self.assertTrue(len(files) >= 1)

    def test_summary(self):
        import pandas as pd
        from pathlib import Path

        pattern = {
            "sample1": ["KO1"],
            "sample2": ["KO1", "KO2"]
        }
        groups = {"sample1": "A", "sample2": "A"}  # no PERMANOVA expected

        metadata = pd.DataFrame({"group": groups})
        metadata.index.name = "sample"
        metadata.index = list(groups.keys())

        outdir = Path(self.tmpdir) / "summary"
        gc.summary(pattern, "Food to Microbe", True, metadata, ["group"], str(outdir))

        files = list(outdir.glob("*Summary.txt"))
        self.assertEqual(len(files), 1)
        text = files[0].read_text()
        self.assertIn("SUMMARY FOR PATTERN: Food to Microbe", text)



if __name__ == "__main__":
    unittest.main()
