#!/usr/bin/env python3
import argparse
import networkx as nx

def read_gfa(gfa_file):
    """
    Read the GFA file and return a tuple of:
      - s_lines: list of S lines (segments)
      - l_lines: list of L lines (links)
    """
    s_lines = []
    l_lines = []
    with open(gfa_file, 'r') as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("S"):
                s_lines.append(line)
            elif line.startswith("L"):
                l_lines.append(line)
    return s_lines, l_lines

def build_graph(s_lines, l_lines):
    """
    Build an undirected graph using S lines for nodes and L lines for edges.
    Only distinct edges are counted.
    """
    G = nx.Graph()
    for line in s_lines:
        parts = line.split("\t")
        if len(parts) >= 2:
            node = parts[1]
            G.add_node(node)
    for line in l_lines:
        parts = line.split("\t")
        if len(parts) >= 6:
            seg1 = parts[1]
            seg2 = parts[3]
            G.add_edge(seg1, seg2)
    return G

def compute_average_degree(G):
    """
    Compute the average degree of graph G.
    For an undirected graph, average_degree = 2 * (#edges) / (#nodes).
    """
    n = G.number_of_nodes()
    if n == 0:
        return 0
    return 2 * G.number_of_edges() / n

def main():
    parser = argparse.ArgumentParser(
        description="Score a GFA graph based on how untangled it is and how many segments remain (lower is better)."
    )
    parser.add_argument("gfa_file", help="Path to the GFA file to be scored.")
    parser.add_argument("--original-nodes", type=int, required=True,
                        help="The original number of segments (nodes) before any removals.")
    parser.add_argument("--verbose", action="store_true", 
                        help="Print detailed metrics.")
    args = parser.parse_args()

    # Read GFA file.
    s_lines, l_lines = read_gfa(args.gfa_file)
    N_current = len(s_lines)
    if args.verbose:
        print(f"Number of segments in current graph: {N_current}")
        print(f"Number of links in current graph: {len(l_lines)}")
    
    # Build graph and compute average degree.
    G = build_graph(s_lines, l_lines)
    avg_degree = compute_average_degree(G)
    if args.verbose:
        print(f"Average degree: {avg_degree:.3f}")
    
    # Compute fraction of segments remaining.
    if args.original_nodes == 0:
        print("Error: The original number of segments must be greater than zero.")
        return
    frac_remaining = N_current / args.original_nodes
    if args.verbose:
        print(f"Fraction of segments remaining: {frac_remaining:.3f}")
    
    # Compute score: lower score means better (untangled and less deletion).
    # Adding 1 to avg_degree avoids division by zero.
    if frac_remaining == 0:
        score = float("inf")
    else:
        score = (avg_degree + 1) / frac_remaining**0.5
    
    print("\nGraph Score (lower = better, meaning more untangled and more segments retained):")
    print(f"{score:.5f}")
    
    if args.verbose:
        print("\nSummary:")
        print(f"  Original segments:   {args.original_nodes}")
        print(f"  Remaining segments:  {N_current}")
        print(f"  Average degree:      {avg_degree:.3f}")
        print(f"  Score:               {score:.5f}")

if __name__ == "__main__":
    main()

