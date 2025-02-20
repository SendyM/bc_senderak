#!/usr/bin/env python3
import argparse
import networkx as nx
import csv
import sys

def read_gfa(gfa_file):
    """
    Read the GFA file and separate S lines (segments) and L lines (links).
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
    Build an undirected NetworkX graph using nodes from S lines and
    edges from L lines.
    
    The S lines are assumed to be in the format:
      S <segment> <sequence> [other tags...]
      
    The L lines are assumed to be in the format:
      L <seg1> <ori1> <seg2> <ori2> <overlap>
      
    Duplicate edges are automatically ignored by the Graph.
    """
    G = nx.Graph()
    # Add nodes from S lines.
    for line in s_lines:
        parts = line.split("\t")
        if len(parts) >= 2:
            node = parts[1]
            G.add_node(node)
    # Add edges from L lines.
    for line in l_lines:
        parts = line.split("\t")
        if len(parts) >= 6:
            seg1 = parts[1]
            seg2 = parts[3]
            G.add_edge(seg1, seg2)
    return G

def compute_metrics(G):
    """
    Compute the degree (distinct edges) and betweenness centrality for each node.
    Returns two dictionaries:
      - degree_dict: {node: degree}
      - bc_dict: {node: betweenness_centrality}
    """
    degree_dict = dict(G.degree())
    bc_dict = nx.betweenness_centrality(G)
    return degree_dict, bc_dict

def main():
    parser = argparse.ArgumentParser(
        description="Output a table with every vertex (segment) from a GFA file along with its distinct degree and betweenness centrality. "
                    "Optionally, sort the table by degree or centrality."
    )
    parser.add_argument("gfa_file", help="Path to the input GFA file.")
    parser.add_argument("--output-file", help="Path to the output CSV file. If not provided, output is printed to stdout.")
    parser.add_argument("--sort-by", choices=["degree", "centrality"], default=None,
                        help="Sort the output table by 'degree' or 'centrality'.")
    parser.add_argument("--descending", action="store_true",
                        help="Sort in descending order (default is ascending).")
    args = parser.parse_args()

    # Read the GFA file.
    s_lines, l_lines = read_gfa(args.gfa_file)
    # Build the graph.
    G = build_graph(s_lines, l_lines)
    # Compute metrics.
    degree_dict, bc_dict = compute_metrics(G)

    # Create table data.
    table = []
    for node in G.nodes():
        table.append({
            "segment": node,
            "degree": degree_dict.get(node, 0),
            "centrality": bc_dict.get(node, 0.0)
        })

    # Sort the table if requested.
    if args.sort_by:
        if args.sort_by == "degree":
            sort_key = lambda x: x["degree"]
        else:  # sort_by centrality
            sort_key = lambda x: x["centrality"]
        table = sorted(table, key=sort_key, reverse=args.descending)

    # Output the table in tab-delimited CSV format.
    fieldnames = ["segment", "degree", "centrality"]
    if args.output_file:
        with open(args.output_file, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for row in table:
                writer.writerow(row)
        print(f"Output written to {args.output_file}")
    else:
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in table:
            writer.writerow(row)

if __name__ == "__main__":
    main()

