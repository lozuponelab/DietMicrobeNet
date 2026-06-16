import unittest
import pandas as pd
from foodb_proc import foodb_nodes_edges as ne

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

rn_dict = { 
    'rn1': {'ORTHOLOGY': [['K00001']], 'EQUATION': [['C1', 'C2'], ['C3']]},
    'rn5': {'ORTHOLOGY': [['K00004']], 'EQUATION': [['C8'], ['C9']]}
}

class MyTestCase(unittest.TestCase): 
    def test_org_abundance(self): 
        abundance, orgs = ne.make_organisms_abundance_dict(microbe_meta_clean=microbe_meta, 
                                          abundance_column='Abundance_RPKs', e_weights=True, orgs=True)
        expected_abundance = {'K00001': 15, 'K00004': 50}
        expected_orgs = {'K00001': 'org1, org2', 'K00004': 'org3'}  
        
        self.assertEqual(abundance, expected_abundance)
        self.assertEqual(orgs, expected_orgs)

    def test_get_organisms_and_abundance(self): 
        kos = ['K00003', 'K00004']
        abundance, orgs = ne.make_organisms_abundance_dict(microbe_meta_clean=microbe_meta, 
                                          abundance_column='Abundance_RPKs', e_weights=True, orgs=True)

        o = ne.get_organisms(kos, orgs)
        a = ne.get_abundance(kos, abundance)
        
        expected_o = ['org3']
        expected_a = 50 

        self.assertCountEqual(o , expected_o)
        self.assertEqual(a, expected_a)
    
    def test_edge_creation(self): 
        abundance, orgs = ne.make_organisms_abundance_dict(microbe_meta_clean=microbe_meta, 
                                          abundance_column='Abundance_RPKs', e_weights=True, orgs=True)

        edge_df, org_comps, all_rxn_comps = ne.build_edges_df(rxns=rn_dict, 
                                    orgs=True, 
                                    e_weights=True,
                                    ko_abundance=abundance, 
                                    ko_organisms=orgs)

        # test how many edges would be made, expect 3
        self.assertEqual(len(edge_df), 3) 
        
        # test dataframe output 
        edge_columns = ['compound1','compound2','reaction','KOs','organisms','abundance']
        self.assertEqual(list(edge_df.columns), edge_columns)

        # test that edge weights become NAs 
        edge_wo_weights, org_comps, all_rxn_comps = ne.build_edges_df(rxns=rn_dict, 
                                    orgs=True, 
                                    e_weights=False,
                                    ko_abundance=abundance, 
                                    ko_organisms=orgs)
        self.assertTrue(edge_wo_weights['abundance'].isnull().all())

        # test that organisms become NAs 
        edge_wo_orgs, org_comps, all_rxn_comps = ne.build_edges_df(rxns=rn_dict, 
                                    orgs=False, 
                                    e_weights=True,
                                    ko_abundance=abundance, 
                                    ko_organisms=orgs)
        self.assertTrue(edge_wo_orgs['organisms'].isnull().all())

        # test number of organism compounds 
        self.assertEqual(len(org_comps), 2)

        # test number of all compounds 
        self.assertEqual(len(all_rxn_comps), 5)

    def test_node_creation(self):
        abundance, orgs = ne.make_organisms_abundance_dict(microbe_meta_clean=microbe_meta, 
                                          abundance_column='Abundance_RPKs', e_weights=True, orgs=True)
        edge_df, org_comps, all_rxn_comps = ne.build_edges_df(rxns=rn_dict, 
                                    orgs=True, 
                                    e_weights=True,
                                    ko_abundance=abundance, 
                                    ko_organisms=orgs) 

        node_df = ne.build_nodes_df(food_meta=food_meta_df, org_comps=org_comps, all_rxn_comps=all_rxn_comps,
                                    frequency=True)

        # test how many nodes would be made, expect 9
        self.assertEqual(len(node_df), 5) 

        # test that all frequencies are either 60, 40, 100 or NA
        self.assertTrue(node_df['freq'].isin([40, 60, 100, pd.NA]).all())

        # make sure that all microbial compounds have no associated food 
        associations = node_df[node_df['origin']=='microbe']
        check = associations['assoc_food'].isna().all() 
        self.assertTrue(check)

        # make sure all microbial compound have no associated freq 
        associations = node_df[node_df['origin']=='microbe']
        check = associations['freq'].isna().all() 
        self.assertTrue(check)

        # make sure all food associated compounds have food assoc 
        associations = node_df[node_df['origin']!='microbe']
        check = associations['assoc_food'].isna().all() 
        self.assertFalse(check)

        # make sure all food compound have no associated freq 
        associations = node_df[node_df['origin']!='microbe']
        check = associations['freq'].isna().all() 
        self.assertFalse(check)

        # make sure weights are na when frequency = False
        node_df = ne.build_nodes_df(food_meta=food_meta_df, org_comps=org_comps, all_rxn_comps=all_rxn_comps,
                                    frequency=False)
        check = node_df['freq'].isna().all() 
        self.assertTrue(check)

if __name__ == "__main__":
    unittest.main()