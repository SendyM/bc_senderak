
#!/usr/bin/env python3
"""Score two GFA graphs (original vs. filtered) with a composite Hairball Index.

Usage:
    python score_graphs.py original.gfa filtered.gfa
"""

import networkx as nx
import sys

# ---------- Helpers ----------------------------------------------------------
def load_gfa(path):
    """
    Parse a GFA file and return a NetworkX undirected graph plus total segment length.
    Nodes carry attribute 'length' (bp).
    """
    G = nx.Graph()
    total_len = 0

    with open(path, "r") as fh:
        for line in fh:
            if line.startswith("S"):                         # Segment
                parts = line.rstrip().split("\t")
                name   = parts[1]
                seq    = parts[2] if len(parts) > 2 else ""
                seglen = len(seq) if seq not in ("", "*") else 0
                G.add_node(name, length=seglen)
                total_len += seglen
            elif line.startswith("L"):                       # Link
                parts = line.rstrip().split("\t")
                u = parts[1]
                v = parts[3]
                G.add_edge(u, v)
    return G, total_len


def graph_metrics(G):
    """
    Compute basic topological metrics needed for scoring.
    Returns dict with:
        edges, nodes, edge_per_node, avg_shortest_path, modularity
    """
    E = G.number_of_edges()
    V = G.number_of_nodes()
    edge_per_node = E / V if V else 0.0

    # average shortest path on largest component (avoid infinities)
    if V > 1 and not nx.is_empty(G):
        largest_cc = max(nx.connected_components(G), key=len)
        sub = G.subgraph(largest_cc).copy()
        try:
            avg_sp = nx.average_shortest_path_length(sub) if sub.number_of_nodes() > 1 else 0.0
        except Exception:
            avg_sp = 0.0
    else:
        avg_sp = 0.0

    # greedy modularity (fast, deterministic)
    if V >= 2:
        communities = list(nx.community.greedy_modularity_communities(G))
        Q = nx.community.modularity(G, communities)
    else:
        Q = 0.0

    return {
        "edges": E,
        "nodes": V,
        "edge_per_node": edge_per_node,
        "avg_shortest_path": avg_sp,
        "modularity": Q,
    }


def compute_score(orig, filt, totals, weights):
    """
    Compute composite Hairball‑like score (lower is better).

    weights: dict with keys
        w_e2n, w_mod, w_sp,
        w_nodes_removed, w_edges_removed, w_len_removed
    """
    # terms computed on filtered graph
    e2n_term  = filt["edge_per_node"]
    mod_term  = 1.0 - filt["modularity"]
    sp_ratio  = (
        filt["avg_shortest_path"] / orig["avg_shortest_path"]
        if orig["avg_shortest_path"] else 1.0
    )

    # removal ratios
    nodes_removed = (orig["nodes"] - filt["nodes"]) / orig["nodes"] if orig["nodes"] else 0.0
    edges_removed = (orig["edges"] - filt["edges"]) / orig["edges"] if orig["edges"] else 0.0
    len_removed   = (
        (totals["orig_len"] - totals["filt_len"]) / totals["orig_len"]
        if totals["orig_len"] else 0.0
    )

    HI = (
        weights["w_e2n"]            * e2n_term +
        weights["w_mod"]            * mod_term +
        weights["w_sp"]             * sp_ratio +
        weights["w_nodes_removed"]  * nodes_removed +
        weights["w_edges_removed"]  * edges_removed +
        weights["w_len_removed"]    * len_removed
    )
    return HI


# ---------- Main -------------------------------------------------------------
def main(orig_gfa, filt_gfa):
    # default weights (can be tweaked)
    weights = {
        "w_e2n":            0.5,
        "w_mod":            0.25,
        "w_sp":             0.012,
        "w_nodes_removed":  0.1,
        "w_edges_removed":  0.1,
        "w_len_removed":    0.05,
    }

    G_orig, len_orig = load_gfa(orig_gfa)
    G_filt, len_filt = load_gfa(filt_gfa)

    m_orig = graph_metrics(G_orig)
    m_filt = graph_metrics(G_filt)

    nodes_removed = (m_orig["nodes"] - m_filt["nodes"]) / m_orig["nodes"] if m_orig["nodes"] else 0.0
    edges_removed = (m_orig["edges"] - m_filt["edges"]) / m_orig["edges"] if m_orig["edges"] else 0.0
    len_removed   = (len_orig - len_filt) / len_orig if len_orig else 0.0

    HI = compute_score(
        m_orig,
        m_filt,
        {"orig_len": len_orig, "filt_len": len_filt},
        weights,
    )

    # -------- Output summary -------------
    print("=== Original graph ===")
    for k, v in m_orig.items():
        print(f"{k:22s}: {v}")
    print(f"{'total_length':22s}: {len_orig}\n")

    print("=== Filtered graph ===")
    for k, v in m_filt.items():
        print(f"{k:22s}: {v}")
    print(f"{'total_length':22s}: {len_filt}\n")

    print("=== Differences ===")
    print(f"nodes_removed        : {m_orig['nodes'] - m_filt['nodes']} "
          f"({nodes_removed*100:.2f}%)")
    print(f"edges_removed        : {m_orig['edges'] - m_filt['edges']} "
          f"({edges_removed*100:.2f}%)")
    print(f"length_removed(bp)   : {len_orig - len_filt} "
          f"({len_removed*100:.2f}%)\n")

    print(f"Composite Hairball score (lower is better): {HI:.4f}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python score_graphs.py original.gfa filtered.gfa")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])

