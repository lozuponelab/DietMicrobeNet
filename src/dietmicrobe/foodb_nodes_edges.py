import json
import pandas as pd 

def data_read_in(reaction_path:str, food_path:str, m_meta:str, 
                 e_weights:bool, orgs:bool, abundance_column:str): 
    """read in all data needed for node and edge dataframe creation 

    Args:
        reaction_path (str): path to AMON rn_dict.json file 
        food_path (str): path to CSV containing food information 
        m_meta (str): path to CSV of microbe KO metadata (taxonomy and abundance)
        e_weights (bool): True or False whether edge weights (abundance) should be included
        orgs (bool): True or False whether organism information will be added 
        abundance_column (str): name of the column where abundance information is located 

    Returns:
        pandas dataframes and a dictionary 
    """

    # get reaction dictionary from AMON output 
    with open(reaction_path, 'r') as json_rn: 
        rxns = json.load(json_rn)

    # get food associated compounds 
    f_comp_df = pd.read_csv(food_path)

    # get microbe metadata 
    microbe_meta = pd.read_csv(m_meta)
    if e_weights and orgs: # both abundance and organisms 
        m_meta_clean = microbe_meta.dropna(subset=['KO', 'taxonomy', abundance_column]).copy() # remove Nan
    elif orgs: # just organisms 
        m_meta_clean = microbe_meta.dropna(subset=['KO', 'taxonomy']).copy() # remove Nan
    elif e_weights: # just abundance 
        m_meta_clean = microbe_meta.dropna(subset=['KO', abundance_column]).copy() # remove Nan
    else: # neither abundance nor organisms 
        m_meta_clean = microbe_meta.dropna(subset=['KO']).copy() # remove Nan
    m_meta_clean['taxonomy'] = m_meta_clean['taxonomy'].astype(str) # convert taxonomy to string 
    print('NAs have been removed from microbe metadata')

    return rxns, f_comp_df, m_meta_clean

def make_organisms_abundance_dict(microbe_meta_clean, abundance_column:str, e_weights:bool, orgs:bool):
    """associated microbe taxonomy and abundance values with KOs 

    Args:
        microbe_meta_clean (pandas df): dataframe containing KOs, taxonomy and abundance 
        abundance_column (str): name of the column containing abundance information 
        e_weights (bool): True or false whether edge weights (abundance will be included)
        orgs (bool): True or false wether organisms will be included in edge information 
    
    Returns:
        two dictionaries: {KO:abundance} and {KO: organisms}
    """
    if e_weights:
        # create dict of KOs with abundance 
        ko_abundance_sum = microbe_meta_clean.groupby('KO', as_index=False).sum(abundance_column)
        ko_abundance_dict = ko_abundance_sum.set_index('KO')[abundance_column].to_dict()
    else:
        ko_abundance_dict = None

    if orgs:
        # create dict of KOs with list of associated organisms -> removed duplicates and sorted alphabetically
        ko_orgs_df = microbe_meta_clean.groupby('KO')['taxonomy'].agg(lambda x: ', '.join(sorted(set(x)))).reset_index()
        ko_orgs_dict = ko_orgs_df.set_index('KO')['taxonomy'].to_dict()
    else:
        ko_orgs_dict = None

    return ko_abundance_dict, ko_orgs_dict

def get_organisms(kos:list, ko_organisms:dict): 
    """from a list of KOs get all associated organisms 

    Args:
        kos (list): list of KOs 
        ko_organisms (dict): dictionary where keys are KOs and values are organisms 

    Returns:
        list: organisms associated with KOs
    """
    # initialize list of organisms 
    organisms = []
    # for each ko get organisms
    organisms.extend(ko_organisms[ko] for ko in kos if ko in ko_organisms)
    
    # return set for no duplications 
    return list(set(organisms))

def get_abundance(kos:list, ko_abundance:dict): 
    """from a list of KOs get total summed abundance 

    Args:
        kos (list): list of KOs
        ko_abundance (dict): dictionary where keys are KOs and values are a number 

    Returns:
        float: total abundance of the list of KOs
    """
    # returning sum of abundance for all KOs
    return sum(ko_abundance[ko] for ko in kos if ko in ko_abundance)

def build_edges_df(rxns:dict, 
                   orgs:bool, 
                   e_weights:bool, 
                   ko_organisms:dict, 
                   ko_abundance:dict):
    """build a dataframe containing edge information

    Args:
        rxns (dict): dictionary of reactions from AMON outputs 
        orgs (bool): whether organism information will be included 
        e_weights (bool): whether abundance information will be included 
        ko_organisms (dict): keys are kos and values a list of associated organisms 
        ko_abundance (dict): keys are kos and values are numbers representing their abundance 

    Returns:
        pandas df: each row represents an edge and associated information 
    """
    edges_list = []
    org_comps = []
    all_rxn_comps = []

    for rxn, info in rxns.items():

        # for edge creation need reaction equation and associated KOs
        equation = info['EQUATION']

        orthology = info['ORTHOLOGY']
        kos = {orthology[index][0] for index in range(len(orthology))}

        # Optional fields: organisms and abundance 
        organisms = get_organisms(kos, ko_organisms) if orgs else pd.NA
        abundance = get_abundance(kos, ko_abundance) if e_weights else pd.NA

        reactants = equation[0]
        products = equation[1]

        org_comps += products
        all_rxn_comps += reactants + products

        for r in reactants:
            for p in products:
                edges_list.append({
                    'compound1': r,
                    'compound2': p,
                    'reaction': rxn,
                    'KOs': list(kos),
                    'organisms': organisms,
                    'abundance': abundance
                })

    return pd.DataFrame(edges_list), list(set(org_comps)), list(set(all_rxn_comps))

def build_nodes_df(food_meta, org_comps:list, all_rxn_comps:list, frequency:bool):
    """build a dataframe containing information about all nodes 

    Args:
        food_meta (pandas df): dataframe containing food KOs and their associated food items
        org_comps (list): list compounds with organism origin 
        all_rxn_comps (list): list of compounds associated with reactions 
        frequency (bool): whether food frequency will be taken into account 

    Returns:
        pandas df: each row represents a node and associated information 
    """
    nodes_list = []

    food_comps = set(food_meta['kegg_id'].to_list())
    all_comps = set(all_rxn_comps + list(food_comps))

    for compound in all_comps:

        # Default to NA
        foods = pd.NA
        freq = pd.NA

        # set origin of compound 
        if compound in food_comps and compound in org_comps:
            origin = 'both'
        elif compound in food_comps:
            origin = 'food'
        elif compound in org_comps:
            origin = 'microbe'
        else:
            origin = 'none'

        # Only assign food associations if origin is food or both
        if origin in ('food', 'both'):
            subset = food_meta[food_meta['kegg_id'] == compound]
            foods = list(subset['name'].unique())
            
            if frequency: # if nodes are weights by frequency of consumption 
                unique_foods = subset.drop_duplicates(subset=['name'])
                freq = unique_foods['food_frequency'].sum()

        nodes_list.append({
            'compound': compound,
            'origin': origin,
            'assoc_food': foods,
            'freq': freq
        })

    return pd.DataFrame(nodes_list)

def summarize_res(rxns:dict, nodes_df, edges_df, directory:str): 
    """create a summary of information for the network 

    Args:
        rxns (dict): keys are reactions and values ae KOs 
        nodes_df (pandas df): dataframe where each row is a node 
        edges_df (pandas df): dataframe where each row is an edge 
        directory (str): file path 
    """
    with open(directory, 'w') as file_object: 
        
        # how many reactions were gathered 
        file_object.write(f'{len(rxns.keys())} number of reactions were found\n')
    
        # how many nodes were made 
        file_object.write(f'{len(nodes_df)} nodes were created\n')

        # how many edges 
        file_object.write(f'{len(edges_df)} edges were created\n\n')

        # how many nodes with none, both, microbial, and food origins 
        file_object.write('Number of nodes per origin category\n')
        vc_str = "\n".join(f"{idx}: {val}" for idx, val in nodes_df['origin'].value_counts().items())
        file_object.write(vc_str)
