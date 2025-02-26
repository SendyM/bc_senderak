#!/usr/bin/env python3
import argparse
import networkx as nx
import csv

def read_gfa(gfa_file):
    s_lines = []
    l_lines = []
    other_lines = []
    with open(gfa_file, 'r') as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("S"):
                s_lines.append(line)
            elif line.startswith("L"):
                l_lines.append(line)
            else:
                other_lines.append(line)
    return s_lines, l_lines, other_lines

def filter_s_lines(s_lines, removed_set):
    new_s_lines = []
    for line in s_lines:
        fields = line.split("\t")
        if len(fields) >= 2:
            seg = fields[1]
            if seg not in removed_set:
                new_s_lines.append(line)
    return new_s_lines

def build_graph(s_lines, l_lines):
    G = nx.Graph()
    for line in s_lines:
        fields = line.split("\t")
        if len(fields) >= 2:
            node = fields[1]
            G.add_node(node)
    for line in l_lines:
        fields = line.split("\t")
        if len(fields) >= 6:
            seg1 = fields[1]
            seg2 = fields[3]
            G.add_edge(seg1, seg2)
    return G

### Reconnect segments using the .geese file ordering ###
def build_links_from_geese(geese_file, removed_set):
    # Dictionary: sample -> list of tuples (atom_nr, class, strand)
    sample_order = {}
    with open(geese_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        # Expected field names: name (or #name), atom_nr, class, strand, start, end
        for row in reader:
            sample = row.get("#name") or row.get("name")
            if not sample:
                continue
            try:
                atom_nr = int(float(row["atom_nr"]))
            except Exception:
                continue
            seg_class = row["class"]
            strand = row["strand"]
            if sample not in sample_order:
                sample_order[sample] = []
            sample_order[sample].append((atom_nr, seg_class, strand))

    new_l_lines_set = set()

    # For each sample, sort entries by atom_nr and then connect consecutive segments
    for sample, entries in sample_order.items():
        sorted_entries = sorted(entries, key=lambda x: x[0])
        filtered_entries = [
            (atom_nr, seg_class, strand)
            for (atom_nr, seg_class, strand) in sorted_entries
            if seg_class not in removed_set
        ]
        for i in range(len(filtered_entries) - 1):
            seg1 = filtered_entries[i][1]
            ori1 = filtered_entries[i][2]
            seg2 = filtered_entries[i+1][1]
            ori2 = filtered_entries[i+1][2]
            new_line = f"L\t{seg1}\t{ori1}\t{seg2}\t{ori2}\t0M"
            new_l_lines_set.add(new_line)
    return list(new_l_lines_set)

### Main iterative removal procedure using geese ordering for bridging ###
def main():
    parser = argparse.ArgumentParser(
        description="Iteratively flag and remove candidate TE nodes from a GFA graph using centrality metrics. "
                    "When reconnecting neighbors after removal, the script uses ordering information from a .geese file "
                    "to determine the proper sequence of segments in each sample."
    )
    parser.add_argument("gfa_file", help="Path to the input GFA file.")
    parser.add_argument("output_gfa", help="Path for the output GFA file (with removed segments and bridged links).")
    parser.add_argument("--geese-file", required=True,
                        help="Path to the .geese file with sample ordering (header: #name, atom_nr, class, strand, start, end).")
    parser.add_argument("--iterations", type=int, default=1,
                        help="Number of iterations to perform (default: 1).")
    parser.add_argument("--degree-threshold", type=int, default=10,
                        help="Minimum degree threshold for flagging a node (default: 10).")
    parser.add_argument("--bc-threshold", type=float, default=0.05,
                        help="Minimum betweenness centrality threshold for flagging a node (default: 0.05).")
    parser.add_argument("--flagged-output-file", type=str, default=None,
                        help="Optional: Output file for detailed information on flagged (removed) nodes.")
    parser.add_argument("--plain-output-file", type=str, default=None,
                        help="Optional: Output file for a plain list of removed node names (one per line).")
    args = parser.parse_args()

    # Read the original GFA file.
    s_lines, l_lines, other_lines = read_gfa(args.gfa_file)
    print(f"Initial graph: {len(s_lines)} segments, {len(l_lines)} links.")

    # This list will record all removed nodes along with the iteration in which they were removed.
    flagged_all = []
    # Cumulative set of removed nodes.
    all_removed = set()

    # Counters for logging removals.
    total_degree_removed = 0
    total_bc_removed = 0
    total_forced_removed = 0

    # Iterative removal loop.
    for iteration in range(args.iterations):
        print(f"\nIteration {iteration+1}:")
        G = build_graph(s_lines, l_lines)
        print(f"  Graph has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")

        degrees = dict(G.degree())
        bc = nx.betweenness_centrality(G)

        candidates = []
        iteration_degree_removed = 0
        iteration_bc_removed = 0

        for node in G.nodes():
            d = degrees.get(node, 0)
            b = bc.get(node, 0.0)
            degree_flag = d >= args.degree_threshold
            bc_flag = b >= args.bc_threshold
            if degree_flag or bc_flag:
                candidates.append((node, d, b))
                if degree_flag:
                    iteration_degree_removed += 1
                if bc_flag:
                    iteration_bc_removed += 1

        if not candidates:
            highest = max(bc, key=bc.get)
            d = degrees.get(highest, 0)
            b = bc.get(highest, 0.0)
            candidates = [(highest, d, b)]
            iteration_forced_removed = 1
            print(f"  No node met the thresholds; forcing removal of '{highest}' (degree={d}, bc={b:.5f}).")
        else:
            iteration_forced_removed = 0
            print(f"  Flagging {len(candidates)} nodes for removal based on thresholds.")

        for cand in candidates:
            flagged_all.append((cand[0], cand[1], cand[2], iteration+1))

        newly_removed = {cand[0] for cand in candidates}
        all_removed |= newly_removed

        total_degree_removed += iteration_degree_removed
        total_bc_removed += iteration_bc_removed
        total_forced_removed += iteration_forced_removed

        print(f"  Removing {len(newly_removed)} new nodes. Total removed so far: {len(all_removed)}")
        s_lines = filter_s_lines(s_lines, all_removed)
        l_lines = build_links_from_geese(args.geese_file, all_removed)
        print(f"  After iteration {iteration+1}: {len(s_lines)} segments remain.")
        if len(s_lines) == 0:
            print("  All segments removed; stopping iterations.")
            break

    # Write final GFA file.
    with open(args.output_gfa, 'w') as out:
        for line in other_lines:
            out.write(line + "\n")
        for line in s_lines:
            out.write(line + "\n")
        for line in l_lines:
            out.write(line + "\n")
    print(f"\nFinal simplified GFA file written to '{args.output_gfa}'.")

    # Log cumulative removal stats.
    print("\nSummary of removals:")
    print(f"  Total nodes removed based on degree threshold: {total_degree_removed}")
    print(f"  Total nodes removed based on betweenness centrality threshold: {total_bc_removed}")
    if total_forced_removed > 0:
        print(f"  Total forced removals: {total_forced_removed}")

    # Optionally, output flagged node information.
    if args.flagged_output_file:
        with open(args.flagged_output_file, 'w') as detail_out:
            detail_out.write("Iteration\tNode\tDegree\tBetweennessCentrality\n")
            for node, d, b, it in flagged_all:
                detail_out.write(f"{it}\t{node}\t{d}\t{b:.5f}\n")
        print(f"Detailed flagged node information written to '{args.flagged_output_file}'.")
    if args.plain_output_file:
        with open(args.plain_output_file, 'w') as plain_out:
            for node, d, b, it in flagged_all:
                plain_out.write(f"{node}\n")
        print(f"Plain list of flagged node names written to '{args.plain_output_file}'.")

if __name__ == "__main__":
    main()

