# All 36 possible patterns (node1 → node2 → node3):
#
# --- diet as position 1 (12 patterns) ---
# 1.  diet → microbe     → host
# 2.  diet → microbe     → hostdiet
# 3.  diet → microbe     → hostmicrobe
# 4.  diet → microbe     → all
# 5.  diet → microbediet → host
# 6.  diet → microbediet → hostdiet
# 7.  diet → microbediet → hostmicrobe
# 8.  diet → microbediet → all
# 9.  diet → all         → host
# 10. diet → all         → hostdiet
# 11. diet → all         → hostmicrobe
# 12. diet → all         → all
#
# --- microbediet as position 1 (12 patterns) ---
# 13. microbediet → microbe     → host
# 14. microbediet → microbe     → hostdiet
# 15. microbediet → microbe     → hostmicrobe
# 16. microbediet → microbe     → all
# 17. microbediet → microbediet → host
# 18. microbediet → microbediet → hostdiet
# 19. microbediet → microbediet → hostmicrobe
# 20. microbediet → microbediet → all
# 21. microbediet → all         → host
# 22. microbediet → all         → hostdiet
# 23. microbediet → all         → hostmicrobe
# 24. microbediet → all         → all
#
# --- all as position 1 (12 patterns) ---
# 25. all → microbe     → host
# 26. all → microbe     → hostdiet
# 27. all → microbe     → hostmicrobe
# 28. all → microbe     → all
# 29. all → microbediet → host
# 30. all → microbediet → hostdiet
# 31. all → microbediet → hostmicrobe
# 32. all → microbediet → all
# 33. all → all         → host
# 34. all → all         → hostdiet
# 35. all → all         → hostmicrobe
# 36. all → all         → all


import pandas as pd
import networkx as nx
from tqdm import tqdm
import argparse as arg


# ------------------------------------------------------
# Pattern query helper
# ------------------------------------------------------

def append_result(results, c1, n1, c2, n2, edge1, c3, n3, edge2):
    results.append({
        # compound 1 (diet or microbediet)
        "compound1_id": c1,
        "compound1_origin": n1["origin"],
        "compound1_assoc_food": n1["assoc_food"],
        "compound1_freq": n1["freq"],

        # edge 1 (c1 → c2)
        "edge1_reaction": edge1["reaction"],
        "edge1_KOs": edge1["KOs"],
        "edge1_organisms": edge1["organisms"],
        "edge1_m_abundance": edge1["m_abundance"],
        "edge1_h_abundance": edge1["h_abundance"],

        # compound 2 (microbe or microbediet)
        "compound2_id": c2,
        "compound2_origin": n2["origin"],
        "compound2_assoc_food": n2["assoc_food"],
        "compound2_freq": n2["freq"],

        # edge 2 (c2 → c3)
        "edge2_reaction": edge2["reaction"],
        "edge2_KOs": edge2["KOs"],
        "edge2_organisms": edge2["organisms"],
        "edge2_m_abundance": edge2["m_abundance"],
        "edge2_h_abundance": edge2["h_abundance"],

        # compound 3 (host, hostdiet, hostmicrobe, or all)
        "compound3_id": c3,
        "compound3_origin": n3["origin"],
        "compound3_assoc_food": n3["assoc_food"],
        "compound3_freq": n3["freq"],
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
            m_abundance=row["m_abundance"], 
            h_abundance=row["h_abundance"]
        )

    # ------------------------------------------------------
    # Pattern matching 
    # ------------------------------------------------------

    print("\n🔍 Running pattern queries...")

    results = []

    valid_n1 = {"diet", "microbediet", "all"}
    valid_n2 = {"microbe", "microbediet", "all"}
    valid_n3 = {"host", "hostdiet", "hostmicrobe", "all"}

    for c1, c2, edge1 in tqdm(G.edges(data=True), total=G.number_of_edges(), ncols=80):

        n1 = G.nodes[c1]
        n2 = G.nodes[c2]

        if n1["origin"] not in valid_n1 or n2["origin"] not in valid_n2:
            continue

        for c3 in G.successors(c2):
            n3 = G.nodes[c3]
            if n3["origin"] in valid_n3:
                edge2 = G[c2][c3]  # fetch edge attributes between c2 and c3
                append_result(results, c1, n1, c2, n2, edge1, c3, n3, edge2)

    # Convert to DataFrame
    df = pd.DataFrame(results)
    print(f"\n📊 Found {len(df)} matching relationships.")

    # Save output
    df.to_csv(args.o, index=False)
    print(f"\n💾 Saved results to: {args.o}")


if __name__ == "__main__":
    main()
