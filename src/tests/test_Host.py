import unittest
import pandas as pd
import networkx as nx
from dietmicrobehost import host_nodes_edges as hne
from dietmicrobehost import host_run_graph as hrg

# create dummy data 
food_meta_df = pd.DataFrame({
    'kegg_id': ['C1', 'C2', 'C8'], 
    'name': ['apple', 'apple', 'cattle'],
    'name_scientific': ['Malus pumila', 'Malus pumila', 'Bos taurus'],
    'food_frequency': [60, 60, 40]
})

microbe_meta = pd.DataFrame({
    'KO': ['K00001', 'K00001', 'K00004'],
    'taxonomy': ['org1', 'org2', 'org3'],
    'Abundance_RPKs': [5, 10, 50]
})

host_meta = pd.DataFrame({
    'kos': ['K00002', 'K00004'],
    'abundance': [1, 2]
})

host_meta_duplicates = pd.DataFrame({
    'kos': ['K00002', 'K00002', 'K00004'],
    'abundance': [1.0, 3.0, 2.0]
})

rn_dict = { 
    'rn1': {'ORTHOLOGY': [['K00001']], 'EQUATION': [['C1', 'C2'], ['C3']]},
    'rn5': {'ORTHOLOGY': [['K00004']], 'EQUATION': [['C8'], ['C9']]}
}

def build_test_graph():
    """Build a small directed graph covering all 12 patterns and some non-matching edges."""
    G = nx.DiGraph()

    # nodes
    G.add_node("C1", origin="diet",        assoc_food=["apple"], freq=60)
    G.add_node("C2", origin="microbediet", assoc_food=["banana"], freq=40)
    G.add_node("C3", origin="microbe",     assoc_food=None, freq=None)
    G.add_node("C4", origin="microbediet", assoc_food=None, freq=None)
    G.add_node("C5", origin="host",        assoc_food=None, freq=None)
    G.add_node("C6", origin="hostdiet",    assoc_food=None, freq=None)
    G.add_node("C7", origin="hostmicrobe", assoc_food=None, freq=None)
    G.add_node("C8", origin="all",         assoc_food=None, freq=None)
    G.add_node("C9", origin="none",        assoc_food=None, freq=None)  # should never match

    edge_attrs = {
        "reaction": "rn1",
        "KOs": ["K00001"],
        "organisms": "org1",
        "m_abundance": 10.0,
        "h_abundance": 5.0
    }
    edge_attrs2 = {
        "reaction": "rn2",
        "KOs": ["K00002"],
        "organisms": "org2",
        "m_abundance": 20.0,
        "h_abundance": 8.0
    }

    # diet → microbe → {host, hostdiet, hostmicrobe, all}  (patterns 1-4)
    G.add_edge("C1", "C3", **edge_attrs)
    G.add_edge("C3", "C5", **edge_attrs2)  # pattern 1
    G.add_edge("C3", "C6", **edge_attrs2)  # pattern 2
    G.add_edge("C3", "C7", **edge_attrs2)  # pattern 3
    G.add_edge("C3", "C8", **edge_attrs2)  # pattern 4

    # diet → microbediet → {host, hostdiet, hostmicrobe, all}  (patterns 5-8)
    G.add_edge("C1", "C4", **edge_attrs)
    G.add_edge("C4", "C5", **edge_attrs2)  # pattern 5
    G.add_edge("C4", "C6", **edge_attrs2)  # pattern 6
    G.add_edge("C4", "C7", **edge_attrs2)  # pattern 7
    G.add_edge("C4", "C8", **edge_attrs2)  # pattern 8

    # microbediet → microbediet → {host, hostdiet, hostmicrobe, all}  (patterns 9-12)
    G.add_edge("C2", "C4", **edge_attrs)   # reuses C4 as middle node
    # C4 already connects to C5-C8 above, so patterns 9-12 are covered

    # non-matching edges (should never appear in results)
    G.add_edge("C1", "C9", **edge_attrs)   # c2 origin is none
    G.add_edge("C9", "C5", **edge_attrs)   # c1 origin is none
    G.add_edge("C3", "C9", **edge_attrs2)  # c3 origin is none

    return G

class MyTestCase(unittest.TestCase): 
    def test_org_abundance(self): 
        abundance, orgs = hne.make_organisms_abundance_dict(microbe_meta_clean=microbe_meta, 
                                            microbe_abundance_col_name='Abundance_RPKs', e_weights=True, orgs=True)
        expected_abundance = {'K00001': 15, 'K00004': 50}
        expected_orgs = {'K00001': 'org1, org2', 'K00004': 'org3'}  
        
        self.assertEqual(abundance, expected_abundance)
        self.assertEqual(orgs, expected_orgs)

    def test_org_abundance_no_weights_no_orgs(self):
        abundance, orgs = hne.make_organisms_abundance_dict(microbe_meta_clean=microbe_meta,
                                            microbe_abundance_col_name='Abundance_RPKs',
                                            e_weights=False, orgs=False)
        self.assertIsNone(abundance)
        self.assertIsNone(orgs)

    def test_host_abundance(self):
        host_ko_abundance = hne.make_host_abundance_dict(host_meta)
        expected_abundance = {'K00002':1, 'K00004':2}
        self.assertEqual(host_ko_abundance, expected_abundance)

    def test_host_abundance_sums_duplicates(self):
        """duplicate KOs should be summed and trigger a warning."""
        with self.assertWarns(UserWarning):
            host_ko_abundance = hne.make_host_abundance_dict(host_meta_duplicates)
        expected_abundance = {'K00002': 4.0, 'K00004': 2.0}
        self.assertEqual(host_ko_abundance, expected_abundance)

    def test_get_organisms_and_abundance(self): 
        kos = ['K00003', 'K00004']
        abundance, orgs = hne.make_organisms_abundance_dict(microbe_meta_clean=microbe_meta, 
                                          microbe_abundance_col_name='Abundance_RPKs', e_weights=True, orgs=True)

        o = hne.get_organisms(kos, orgs)
        a = hne.get_abundance(kos, abundance)
        
        expected_o = ['org3']
        expected_a = 50 

        self.assertCountEqual(o , expected_o)
        self.assertEqual(a, expected_a)

    def test_edge_creation(self): 
        m_abundance, orgs = hne.make_organisms_abundance_dict(microbe_meta_clean=microbe_meta, 
                                        microbe_abundance_col_name='Abundance_RPKs', 
                                        e_weights=True, orgs=True)
        
        h_abundance = hne.make_host_abundance_dict(h_meta_clean=host_meta)

        edge_df, microbe_comps, host_comps, all_rxn_comps = hne.build_edges_df(
                                    rxns=rn_dict, 
                                    orgs=True, 
                                    e_weights=True, 
                                    m_ko_abundance=m_abundance,
                                    m_ko_organisms=orgs,
                                    h_ko_abundance=h_abundance)

        # test how many edges would be made, expect 3
        self.assertEqual(len(edge_df), 3) 
        
        # test dataframe output 
        edge_columns = ['compound1', 'compound2', 'reaction', 'KOs', 'organisms', 'm_abundance', 'h_abundance']
        self.assertEqual(list(edge_df.columns), edge_columns)

        # test that edge weights become NAs 
        edge_wo_weights, _, _, _ = hne.build_edges_df(rxns=rn_dict, 
                                    orgs=True, 
                                    e_weights=False,
                                    m_ko_abundance=m_abundance,
                                    m_ko_organisms=orgs,
                                    h_ko_abundance=h_abundance)
        self.assertTrue(edge_wo_weights['m_abundance'].isnull().all())

        # test that organisms become NAs 
        edge_wo_orgs, microbe_comps, host_comps, all_rxn_comps = hne.build_edges_df(rxns=rn_dict, 
                                    orgs=False, 
                                    e_weights=True,
                                    m_ko_abundance=m_abundance,
                                    m_ko_organisms=orgs,
                                    h_ko_abundance=h_abundance)
        self.assertTrue(edge_wo_orgs['organisms'].isnull().all())

        # test number of organism compounds 
        self.assertEqual(len(microbe_comps), 2)

        # test number of all compounds 
        self.assertEqual(len(all_rxn_comps), 5)

    def test_node_origins(self):
        """test all 8 origin combinations using existing dummy data."""
        m_abundance, orgs = hne.make_organisms_abundance_dict(microbe_meta_clean=microbe_meta,
                                            microbe_abundance_col_name='Abundance_RPKs',
                                            e_weights=True, orgs=True)
        h_abundance = hne.make_host_abundance_dict(host_meta)

        _, microbe_comps, host_comps, all_rxn_comps = hne.build_edges_df(
                                    rxns=rn_dict,
                                    orgs=True,
                                    e_weights=True,
                                    m_ko_abundance=m_abundance,
                                    m_ko_organisms=orgs,
                                    h_ko_abundance=h_abundance)

        node_df = hne.build_nodes_df(
                    food_meta=food_meta_df,
                    microbe_comps=microbe_comps,
                    host_comps=host_comps,
                    all_rxn_comps=all_rxn_comps,
                    frequency=True)

        # helper to get origin for a compound
        def get_origin(compound):
            return node_df[node_df['compound'] == compound]['origin'].iloc[0]

        # C1, C2 are in food and are reactants only -> diet
        self.assertEqual(get_origin('C1'), 'diet')
        self.assertEqual(get_origin('C2'), 'diet')

        # C8 is in food and is a reactant only -> diet
        self.assertEqual(get_origin('C8'), 'diet')

        # C3 is produced by K00001 (microbe only) -> microbe
        self.assertEqual(get_origin('C3'), 'microbe')

        # C9 is produced by K00004 (both microbe and host) -> hostmicrobe
        self.assertEqual(get_origin('C9'), 'hostmicrobe')

def test_node_food_associations(self):
    """compounds with diet origin should have food names and frequency populated."""
    _, microbe_comps, host_comps, all_rxn_comps = hne.build_edges_df(
                                rxns=rn_dict,
                                orgs=True,
                                e_weights=True,
                                m_ko_abundance=hne.make_host_abundance_dict(host_meta),
                                m_ko_organisms={},
                                h_ko_abundance=hne.make_host_abundance_dict(host_meta))

    node_df = hne.build_nodes_df(
                food_meta=food_meta_df,
                microbe_comps=microbe_comps,
                host_comps=host_comps,
                all_rxn_comps=all_rxn_comps,
                frequency=True)

    c1_row = node_df[node_df['compound'] == 'C1'].iloc[0]

    # food name should be populated
    self.assertEqual(c1_row['assoc_food'], ['apple'])

    # frequency should be summed (food_meta_df has two apple rows each with freq 60, 
    # but drop_duplicates on name means only one counted -> 60)
    self.assertEqual(c1_row['freq'], 60)

def test_node_non_diet_has_no_food_info(self):
    """compounds with no diet origin should have NA for food and frequency."""
    _, microbe_comps, host_comps, all_rxn_comps = hne.build_edges_df(
                                rxns=rn_dict,
                                orgs=True,
                                e_weights=True,
                                m_ko_abundance=hne.make_host_abundance_dict(host_meta),
                                m_ko_organisms={},
                                h_ko_abundance=hne.make_host_abundance_dict(host_meta))

    node_df = hne.build_nodes_df(
                food_meta=food_meta_df,
                microbe_comps=microbe_comps,
                host_comps=host_comps,
                all_rxn_comps=all_rxn_comps,
                frequency=True)

    c3_row = node_df[node_df['compound'] == 'C3'].iloc[0]
    self.assertTrue(pd.isna(c3_row['assoc_food']))
    self.assertTrue(pd.isna(c3_row['freq']))

def test_node_count(self):
    """total nodes should equal unique compounds across reactions and food."""
    _, microbe_comps, host_comps, all_rxn_comps = hne.build_edges_df(
                                rxns=rn_dict,
                                orgs=True,
                                e_weights=True,
                                m_ko_abundance=hne.make_host_abundance_dict(host_meta),
                                m_ko_organisms={},
                                h_ko_abundance=hne.make_host_abundance_dict(host_meta))

    node_df = hne.build_nodes_df(
                food_meta=food_meta_df,
                microbe_comps=microbe_comps,
                host_comps=host_comps,
                all_rxn_comps=all_rxn_comps,
                frequency=False)

    # rn_dict has C1,C2,C3,C8,C9 from reactions + food adds no new compounds 
    # (C1, C2, C8 already in reactions)
    self.assertEqual(len(node_df), 5)

def _run_pattern_match(self, G):
        """Run the pattern matching loop and return results list."""
        results = []
        valid_n1 = {"diet", "microbediet"}
        valid_n2 = {"microbe", "microbediet"}
        valid_n3 = {"host", "hostdiet", "hostmicrobe", "all"}

        for c1, c2, edge1 in G.edges(data=True):
            n1 = G.nodes[c1]
            n2 = G.nodes[c2]
            if n1["origin"] not in valid_n1 or n2["origin"] not in valid_n2:
                continue
            for c3 in G.successors(c2):
                n3 = G.nodes[c3]
                if n3["origin"] in valid_n3:
                    edge2 = G[c2][c3]
                    hrg.append_result(results, c1, n1, c2, n2, edge1, c3, n3, edge2)

        return results

def test_no_false_positives(self):
    """nodes with origin none should never appear in results."""
    G = build_test_graph()
    results = self._run_pattern_match(G)
    for row in results:
        self.assertNotEqual(row["compound1_origin"], "none")
        self.assertNotEqual(row["compound2_origin"], "none")
        self.assertNotEqual(row["compound3_origin"], "none")

def test_all_12_patterns_matched(self):
    """all 12 origin combination patterns should appear in results."""
    G = build_test_graph()
    results = self._run_pattern_match(G)

    found_patterns = {
        (r["compound1_origin"], r["compound2_origin"], r["compound3_origin"])
        for r in results
    }

    expected_patterns = {
        ("diet",        "microbe",     "host"),
        ("diet",        "microbe",     "hostdiet"),
        ("diet",        "microbe",     "hostmicrobe"),
        ("diet",        "microbe",     "all"),
        ("diet",        "microbediet", "host"),
        ("diet",        "microbediet", "hostdiet"),
        ("diet",        "microbediet", "hostmicrobe"),
        ("diet",        "microbediet", "all"),
        ("microbediet", "microbediet", "host"),
        ("microbediet", "microbediet", "hostdiet"),
        ("microbediet", "microbediet", "hostmicrobe"),
        ("microbediet", "microbediet", "all"),
    }

    self.assertEqual(found_patterns, expected_patterns)

def test_empty_graph_returns_no_results(self):
    """an empty graph should produce no results."""
    G = nx.DiGraph()
    results = self._run_pattern_match(G)
    self.assertEqual(len(results), 0)

def test_no_matching_origins_returns_no_results(self):
    """a graph with only none-origin nodes should produce no results."""
    G = nx.DiGraph()
    G.add_node("X", origin="none", assoc_food=None, freq=None)
    G.add_node("Y", origin="none", assoc_food=None, freq=None)
    G.add_edge("X", "Y", reaction="rn1", KOs=[], organisms=None,
                m_abundance=0, h_abundance=0)
    results = self._run_pattern_match(G)
    self.assertEqual(len(results), 0)

if __name__ == "__main__":
    unittest.main()