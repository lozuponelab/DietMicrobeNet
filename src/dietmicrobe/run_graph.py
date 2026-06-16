# Original script built using Neo4j which is much more efficient w >1M nodes
# most graphs created here are small so it was replaced with networkx

import pandas as pd
import networkx as nx
from tqdm import tqdm
import argparse as arg


# ------------------------------------------------------
# Pattern query helper
# ------------------------------------------------------

def append_result(results, c1, n1, c2, n2, edge):
    """Add matched edge and node attributes to result list."""
    results.append({
        "compound1_id": c1,
        "compound1_origin": n1["origin"],
        "compound1_assoc_food": n1["assoc_food"],
        "compound1_freq": n1["freq"],

        "reaction": edge["reaction"],
        "KOs": edge["KOs"],
        "organisms": edge["organisms"],
        "abundance": edge["abundance"],

        "compound2_id": c2,
        "compound2_origin": n2["origin"],
        "compound2_assoc_food": n2["assoc_food"],
        "compound2_freq": n2["freq"]
    })


# ------------------------------------------------------
# Main
# ------------------------------------------------------

def main():
    parser = arg.ArgumentParser(description="In-memory graph builder + pattern queries")
    parser.add_argument('--n', required=True, help='Node CSV file')
    parser.add_argument('--e', required=True, help='Edge CSV file')
    parser.add_argument('--o', required=True, help='Output CSV file')
    args = parser.parse_args()

    # Load CSVs
    print("📄 Loading CSV files...")
    nodes_df = pd.read_csv(args.n)
    edges_df = pd.read_csv(args.e)
    print(f" → Loaded {len(nodes_df)} nodes and {len(edges_df)} edges")

    # Create a directed graph
    print("\n🔧 Building graph in memory...")
    G = nx.DiGraph()

    # Add nodes
    for _, row in tqdm(nodes_df.iterrows(), total=len(nodes_df), desc="Adding nodes", ncols=80):
        G.add_node(
            row["compound"],
            origin=row["origin"],
            assoc_food=row["assoc_food"],
            freq=row["freq"]
        )

    # Add edges
    for _, row in tqdm(edges_df.iterrows(), total=len(edges_df), desc="Adding edges", ncols=80):
        G.add_edge(
            row["compound1"],
            row["compound2"],
            reaction=row["reaction"],
            KOs=row["KOs"],
            organisms=row["organisms"],
            abundance=row["abundance"]
        )

    # ------------------------------------------------------
    # Pattern matching (Cypher → Python)
    # ------------------------------------------------------

    print("\n🔍 Running pattern queries...")

    results = []

    for c1, c2, edge in tqdm(G.edges(data=True), total=G.number_of_edges(), ncols=80):

        n1 = G.nodes[c1]
        n2 = G.nodes[c2]

        # 1. food → microbe
        if n1["origin"] == "food" and n2["origin"] == "microbe":
            append_result(results, c1, n1, c2, n2, edge)

        # 2. food → both
        if n1["origin"] == "food" and n2["origin"] == "both":
            append_result(results, c1, n1, c2, n2, edge)

        # 3. both → both
        if n1["origin"] == "both" and n2["origin"] == "both":
            append_result(results, c1, n1, c2, n2, edge)

    # Convert to DataFrame
    df = pd.DataFrame(results)
    print(f"\n📊 Found {len(df)} matching relationships.")

    # Save output
    df.to_csv(args.o, index=False)
    print(f"\n💾 Saved results to: {args.o}")


if __name__ == "__main__":
    main()
