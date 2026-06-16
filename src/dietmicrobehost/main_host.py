import host_nodes_edges as hne
import argparse as arg 
import matplotlib.pyplot as plt 
import seaborn as sns 
import warnings

def main(): 
    parser = arg.ArgumentParser(description='Create two dataframes of node and edge information from food metabolomes and AMON outputs of both host and microbes')
   
    parser.add_argument('--f', type=str, required=True, 
        help='file path to food compound information taken from comp_FoodDB.R')
    parser.add_argument('--r', type=str, required=True, 
        help='file path to rn_dict.json AMON output')
    parser.add_argument('--h_meta', type=str, required=True, 
        help='filepath to CSV containing host KOs and abundance information')
    parser.add_argument('--m_meta', type=str, required=True, 
        help='filepath to CSV containing KOs, taxonomy, and abundance information')
    parser.add_argument('--n_weights', action='store_true',  
        help='True or False whether node weights (food frequencies) will be used')
    parser.add_argument('--e_weights', action='store_true', 
        help='True or False whether edge weights (abundance measures) will be used')
    parser.add_argument('--org', action='store_true', 
        help='True or False whether organism information is to be included')
    parser.add_argument('--microbe_a', type=str, required=False, default='',
        help='column in microbe metadata corresponding to the abundance information')
    parser.add_argument('--o', type=str, required=True, 
        help='directory where outputs will go')

    args = parser.parse_args()

    # read in the data 
    rxns, f_comp_df, m_meta_clean, h_meta_clean = hne.data_read_in(
                reaction_path=args.r, 
                food_path=args.f,
                m_meta=args.m_meta,
                h_meta=args.h_meta,
                e_weights=args.e_weights,
                orgs=args.org,
                microbe_abundance_col_name=args.microbe_a)

    # get microbe ko:org and ko:abundance dicts
    microbe_ko_abundance_dict, microbe_ko_orgs_dict = hne.make_organisms_abundance_dict(
                                                                microbe_meta_clean=m_meta_clean,
                                                                microbe_abundance_col_name=args.microbe_a, 
                                                                e_weights=args.e_weights,
                                                                orgs=args.org)

    # get host:abundance dict 
    host_ko_abundance_dict = hne.make_host_abundance_dict(h_meta_clean)

    # create edge dataframe 
    edge_df, microbe_comps, host_comps, all_rxn_comps = hne.build_edges_df(
                                        rxns=rxns, 
                                        orgs=args.org, 
                                        e_weights=args.e_weights, 
                                        m_ko_abundance=microbe_ko_abundance_dict,
                                        m_ko_organisms=microbe_ko_orgs_dict,
                                        h_ko_abundance=host_ko_abundance_dict)

    if len(edge_df) == 0:
        warnings.warn('There are no edges found, please check files to make sure they are correct ones')

    # create node dataframe 
    node_df = hne.build_nodes_df(
                            food_meta=f_comp_df, 
                            microbe_comps=microbe_comps, 
                            host_comps=host_comps, 
                            all_rxn_comps=all_rxn_comps, 
                            frequency=args.n_weights)

    if len(node_df) == 0:
        warnings.warn('There are no nodes found, please check files to make sure they are correct ones')

    # summarize results 
    hne.summarize_res(rxns=rxns, nodes_df=node_df, edges_df=edge_df, directory=args.o + 'network_summary.txt')

    # create distribution plot of food frequency
    if args.n_weights:
        plt.figure(figsize=(8, 6))
        sns.histplot(node_df['freq'].dropna(), kde=True)
        plt.title('Distribution of Food Frequencies Associated with Nodes')
        plt.xlabel('Food Frequency Value')
        plt.ylabel('Frequency/Density')
        plt.savefig(args.o + 'FoodFrequencyDistribution.png')
        plt.close()

    # create distribution plots of microbe and host abundance
    if args.e_weights:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        sns.histplot(edge_df['m_abundance'].dropna(), kde=True, ax=axes[0])
        axes[0].set_title('Distribution of Microbe KO Abundance')
        axes[0].set_xlabel('Microbe Abundance')
        axes[0].set_ylabel('Frequency/Density')

        sns.histplot(edge_df['h_abundance'].dropna(), kde=True, ax=axes[1])
        axes[1].set_title('Distribution of Host KO Abundance')
        axes[1].set_xlabel('Host Abundance')
        axes[1].set_ylabel('Frequency/Density')

        plt.tight_layout()
        plt.savefig(args.o + 'AbundanceDistribution.png')
        plt.close()

    # create CSVs 
    node_df.to_csv(args.o + 'nodes_df.csv', index=False)
    edge_df.to_csv(args.o + 'edges_df.csv', index=False)

    exit()

if __name__ == "__main__":
    main()