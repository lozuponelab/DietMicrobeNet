import foodb_nodes_edges as ne 
import argparse as arg
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import warnings

def main():
    parser = arg.ArgumentParser(description='Create two dataframes of node and edge information from food metabolomes and AMON outputs')
   
    parser.add_argument('--f', type=str,required=True, 
        help='file path to food compound information taken from comp_FoodDB.R')
    parser.add_argument('--r', type=str, required=True, 
        help='file path to rn_dict.json AMON output')
    parser.add_argument('--m_meta', type=str, required=True, 
        help='filepath to CSV containing KOs, taxonomy, and abundance information')
    parser.add_argument('--n_weights', action='store_true',  
        help='True or False wether node weights (food frequencies) will be used')
    parser.add_argument('--e_weights', action='store_true', 
        help='True or False wether edge weights (abundance measures) will be used')
    parser.add_argument('--org', action='store_true', 
        help='True or False whether organism information is to be included')
    parser.add_argument('--a', type=str, required=False, default='',
        help='column in microbe metadata corresponding to the abundance information')
    parser.add_argument('--o', type=str, required=True, 
        help='directory where outputs will go')

    args = parser.parse_args()

    # read in data 
    rxns, f_meta, m_meta = ne.data_read_in(reaction_path=args.r, food_path=args.f, m_meta=args.m_meta, 
                                           e_weights=args.e_weights, orgs=args.org, abundance_column=args.a)
    
    # get taxonomy and abundance for each KO
    ko_abundance, ko_taxonomy = ne.make_organisms_abundance_dict(microbe_meta_clean=m_meta, abundance_column=args.a, 
                                                                 e_weights=args.e_weights, orgs=args.org)
    
    # create edges and get lists of compounds 
    edges_df, org_comps, all_rxn_comps = ne.build_edges_df(rxns=rxns, orgs=args.org, e_weights=args.e_weights, 
                                                           ko_abundance=ko_abundance, ko_organisms=ko_taxonomy)
    if len(edges_df) == 0:
        warnings.warn('There are no edges found, please check files to make sure they are correct ones')


    # create nodes 
    nodes_df = ne.build_nodes_df(food_meta=f_meta, org_comps=org_comps, all_rxn_comps=all_rxn_comps, frequency=args.n_weights)
    if len(nodes_df) == 0:
        warnings.warn('There are no nodes found, please check files to make sure they are correct ones')
    
    # summarize results 
    ne.summarize_res(rxns=rxns, nodes_df=nodes_df, edges_df=edges_df, directory=args.o + 'network_summary.txt')
    
    # create distribution plots of food frequency and abundance 
    if args.n_weights:
        plt.figure(figsize=(8,6))
        sns.histplot(nodes_df['freq'], kde=True)
        plt.title('Distribution of Food Frequencies Associated with Nodes')
        plt.xlabel('Food Frequency Value')
        plt.ylabel('Frequency/Density')
        plt.savefig(args.o + 'M_FoodFrequencyDistribution.png')
        plt.close()

    if args.e_weights:
        plt.figure(figsize=(8,6))
        sns.histplot(edges_df['abundance'], kde=True)
        plt.title('Distribution of KO Abundance Associated with Edges')
        plt.xlabel('Abundance')
        plt.ylabel('Frequency/Density')
        plt.savefig(args.o + 'M_AbundanceDistribution.png')
        plt.close()
    
    # create CSVs 
    nodes_df.to_csv(args.o + 'M_nodes_df.csv', index=False)
    edges_df.to_csv(args.o + 'M_edges_df.csv', index=False)

    exit()

if __name__ == "__main__":
    main()