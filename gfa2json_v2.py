#!/usr/bin/env python3

import argparse
import json
from collections import defaultdict

def extract_color(tags):
    """
    From a list of GFA tags (e.g. ["SR:Z:...", "CL:z:#ff00ff", ...]),
    return the first CL:z: color code (e.g. "#ff00ff"), or None.
    """
    for tag in tags:
        if tag.startswith("CL:z:"):
            return tag.split("CL:z:")[1]
    return None

def read_geese_for_usage(geese_file):
    """
    Parse .geese and build a mapping atom_id -> set(of genomes that mention it).
    """
    node_to_genomes = defaultdict(set)
    with open(geese_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip("\n").split('\t')
            if len(parts) < 6:
                continue
            genome, _, atom, _, _, _ = parts[:6]
            node_to_genomes[atom].add(genome)
    return node_to_genomes

def read_gfa(gfa_file):
    """
    Read a GFA, returning:
      - node_set:    set of all segment IDs
      - edges:       list of {"source": id1, "target": id2}
      - color_map:   dict id -> "#rrggbb" if found, else absent
    """
    node_set = set()
    edges = []
    color_map = {}
    with open(gfa_file, 'r') as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split('\t')
            if line.startswith('S'):
                # S <id> <sequence> [tags...]
                seg = parts[1]
                node_set.add(seg)
                tags = parts[3:]
                color = extract_color(tags)
                if color:
                    color_map[seg] = color
            elif line.startswith('L'):
                # L <src> <ori> <dst> <ori> <overlap> [tags...]
                src, dst = parts[1], parts[3]
                node_set.update([src, dst])
                edges.append({"source": src, "target": dst})
    return node_set, edges, color_map

def build_json(node_set, edges, node_to_genomes, color_map, out_json):
    """
    Combine everything into a single JSON:
      {
        "nodes": [
          {"id": "...", "usageValue": 5, "color": "#ff00ff"},
          ...
        ],
        "links": [ ... ]
      }
    """
    # find overall usageValue range
    usage_vals = [len(node_to_genomes[n]) for n in node_set]
    # but we don't need min/max here in Python, it's for JS
    
    nodes = []
    for n in sorted(node_set):
        usage = len(node_to_genomes.get(n, []))
        node_obj = {
            "id": n,
            "usageValue": usage
        }
        if n in color_map:
            node_obj["color"] = color_map[n]
        nodes.append(node_obj)

    graph = {"nodes": nodes, "links": edges}
    with open(out_json, 'w') as out:
        json.dump(graph, out, indent=2)

    print(f"Written {len(nodes)} nodes and {len(edges)} links to {out_json}")

def main():
    p = argparse.ArgumentParser(
        description="Convert GFA + .geese into a single JSON with color & usageValue"
    )
    p.add_argument("gfa", help="Input GFA file")
    p.add_argument("geese", help="Input .geese file")
    p.add_argument("outjson", help="Output JSON for D3")
    args = p.parse_args()

    node_to_genomes = read_geese_for_usage(args.geese)
    node_set, edges, color_map = read_gfa(args.gfa)
    build_json(node_set, edges, node_to_genomes, color_map, args.outjson)

if __name__ == "__main__":
    main()
