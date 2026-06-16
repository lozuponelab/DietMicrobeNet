import json 
import pandas as pd 
import warnings

def data_read_in(
        reaction_path:str, 
        food_path:str, 
        m_meta:str, 
        h_meta:str,
        e_weights:bool, 
        orgs:bool, 
        microbe_abundance_col_name:str
): 
    """read in all data needed for node and edge dataframe creation 

    Args:
        reaction_path (str): path to AMON rn_dict.json file 
        food_path (str): path to CSV containing food information 
        m_meta (str): path to CSV of microbe KO metadata (taxonomy and abundance)
        h_meta (str): path to CSV of host KO metadata (abundance)
        e_weights (bool): True or False whether edge weights (abundance) should be included
        orgs (bool): True or False whether organism information will be added 
        microbe_abundance_col_name (str): name of the column where abundance information is located 

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
        m_meta_clean = microbe_meta.dropna(subset=['KO', 'taxonomy', microbe_abundance_col_name]).copy() # remove Nan
    elif orgs: # just organisms 
        m_meta_clean = microbe_meta.dropna(subset=['KO', 'taxonomy']).copy() # remove Nan
    elif e_weights: # just abundance 
        m_meta_clean = microbe_meta.dropna(subset=['KO', microbe_abundance_col_name]).copy() # remove Nan
    else: # neither abundance nor organisms 
        m_meta_clean = microbe_meta.dropna(subset=['KO']).copy() # remove Nan
    m_meta_clean['taxonomy'] = m_meta_clean['taxonomy'].astype(str) # convert taxonomy to string 
    print('NAs have been removed from microbe metadata')

    # get host metadata 
    h_meta_clean = pd.read_csv(h_meta)

    # Check column count
    if h_meta_clean.shape[1] != 2:
        raise ValueError(f"Expected 2 columns, got {h_meta_clean.shape[1]}")

    # Validate column 1 (position 0): non-null values must be strings
    h_meta_clean.iloc[:, 0] = h_meta_clean.iloc[:, 0].where(
        h_meta_clean.iloc[:, 0].isnull(), h_meta_clean.iloc[:, 0].astype(str)
    )

    # Validate column 2 (position 1): non-null values must be numeric (raises on bad values)
    non_null_mask = h_meta_clean.iloc[:, 1].notnull()
    try:
        pd.to_numeric(h_meta_clean.iloc[:, 1][non_null_mask], errors="raise")
    except ValueError as e:
        raise ValueError(f"Column 2 contains non-numeric values: {e}")

    h_meta_clean.iloc[:, 1] = pd.to_numeric(h_meta_clean.iloc[:, 1])

    # warn and remove NAs only after confirming column 2 values are valid
    count_no_abundance = h_meta_clean.iloc[:, 1].isna().sum()
    if count_no_abundance > 0:
        warnings.warn(f"There are {count_no_abundance} host KOs with no associated abundance. These will be removed.")
        h_meta_clean = h_meta_clean.dropna(subset=[h_meta_clean.columns[1]]).reset_index(drop=True)

    return rxns, f_comp_df, m_meta_clean, h_meta_clean
def make_organisms_abundance_dict(microbe_meta_clean, 
                                  microbe_abundance_col_name:str, 
                                  e_weights:bool, 
                                  orgs:bool):
    """associated microbe taxonomy and abundance values with KOs 

    Args:
        microbe_meta_clean (pandas df): dataframe containing KOs, taxonomy and abundance 
        microbe_abundance_col_name (str): name of the column containing abundance information 
        e_weights (bool): True or false whether edge weights (abundance will be included)
        orgs (bool): True or false wether organisms will be included in edge information 
    
    Returns:
        two dictionaries: {KO:abundance} and {KO: organisms}
    """
    if e_weights:
        # create dict of KOs with abundance 
        ko_abundance_sum = microbe_meta_clean.groupby('KO', as_index=False).sum(microbe_abundance_col_name)
        ko_abundance_dict = ko_abundance_sum.set_index('KO')[microbe_abundance_col_name].to_dict()
    else:
        ko_abundance_dict = None

    if orgs:
        # create dict of KOs with list of associated organisms -> removed duplicates and sorted alphabetically
        ko_orgs_df = microbe_meta_clean.groupby('KO')['taxonomy'].agg(lambda x: ', '.join(sorted(set(x)))).reset_index()
        ko_orgs_dict = ko_orgs_df.set_index('KO')['taxonomy'].to_dict()
    else:
        ko_orgs_dict = None

    return ko_abundance_dict, ko_orgs_dict

def make_host_abundance_dict(h_meta_clean): 
    """convert host metadata to dictionary, sums abundance if duplicated KOs, will warn if there are duplicates 

    Args:
        h_meta_clean (pandas df): created during data_read_in() or pandas dataframe of host metadata 

    Returns:
        dict: host kos and abundance {ko:abundance}
    """
    # warn if duplicate KOs
    num_duplicates = h_meta_clean.iloc[:, 0].duplicated().sum()
    if num_duplicates > 0:
         warnings.warn(f"There are {num_duplicates} KOs that are duplicated.")

    # create dict of KOs with summed abundance
    first_col = h_meta_clean.columns[0]
    second_col = h_meta_clean.columns[1]
    ko_abundance_sum = h_meta_clean.groupby(first_col)[second_col].sum()
    return ko_abundance_sum.to_dict()

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
                   m_ko_organisms:dict, 
                   m_ko_abundance:dict, 
                   h_ko_abundance:dict):

    edges_list = []
    microbe_comps = []
    host_comps = []
    all_rxn_comps = []

    for rxn, info in rxns.items():

        equation = info['EQUATION']
        orthology = info['ORTHOLOGY']
        kos = {orthology[index][0] for index in range(len(orthology))}

        # Optional fields: organisms and microbe abundance 
        organisms = get_organisms(kos, m_ko_organisms) if orgs else pd.NA
        m_abundance = get_abundance(kos, m_ko_abundance) if e_weights else pd.NA
        h_abundance = get_abundance(kos, h_ko_abundance)

        reactants = equation[0]
        products = equation[1]

        all_rxn_comps += reactants + products

        # track products by source separately
        if any(ko in m_ko_abundance for ko in kos):
            microbe_comps += products
        if any(ko in h_ko_abundance for ko in kos):
            host_comps += products

        for r in reactants:
            for p in products:
                edges_list.append({
                    'compound1': r,
                    'compound2': p,
                    'reaction': rxn,
                    'KOs': list(kos),
                    'organisms': organisms,
                    'm_abundance': m_abundance, 
                    'h_abundance': h_abundance
                })

    return pd.DataFrame(edges_list), list(set(microbe_comps)), list(set(host_comps)), list(set(all_rxn_comps))

# TODO: build nodes df 
def build_nodes_df(food_meta, 
                   microbe_comps: list,
                   host_comps: list,
                   all_rxn_comps: list, 
                   frequency: bool):

    nodes_list = []

    food_comps = set(food_meta['kegg_id'].to_list())
    microbe_comps = set(microbe_comps)
    host_comps = set(host_comps)
    all_comps = set(all_rxn_comps + list(food_comps))

    origin_map = {
        frozenset():                            'none',
        frozenset({'diet'}):                    'diet',
        frozenset({'microbe'}):                 'microbe',
        frozenset({'host'}):                    'host',
        frozenset({'host', 'diet'}):            'hostdiet',
        frozenset({'host', 'microbe'}):         'hostmicrobe',
        frozenset({'microbe', 'diet'}):         'microbediet',
        frozenset({'host', 'microbe', 'diet'}): 'all',
    }

    for compound in all_comps:

        sources = set()
        if compound in food_comps:
            sources.add('diet')
        if compound in microbe_comps:
            sources.add('microbe')
        if compound in host_comps:
            sources.add('host')

        origin = origin_map[frozenset(sources)]

        foods = pd.NA
        freq = pd.NA

        if 'diet' in sources:
            subset = food_meta[food_meta['kegg_id'] == compound]
            foods = list(subset['name'].unique())
            if frequency:
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