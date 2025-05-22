#!/usr/bin/env python3
"""
filter_geese.py

1) Parse .geese into:
     • occ: atom → list of (genome, start, end, strand)
     • genome_order: genome → list of (atom, strand)
     • raw_lines for later rewriting
2) Load FASTA, build per‐atom metadata (sequence, depth, length, duplicated)
3) Compute each atom’s “unique context” count and globally filter on depth/length/dup
4) Rebuild filtered genome_order
5) Score all (a,b) pairs, **collecting for each** (a,b) a map: context_tuple → set(genomes)
6) **New**: Context‐based filtering: for each pair, for each atom appearing in any context,
   keep it only in the single context in which it appears in the most genomes;
   remove it from all other contexts/genomes.
7) Rewrite a single final .geese applying both the global filter and the context‐based removals.
"""
import argparse
import logging
from collections import defaultdict
from dataclasses import dataclass

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

@dataclass
class GfaSegment:
    name: str
    sequence: str
    depth: int
    length: int
    duplicated: bool
    unique_context: int

# -- FASTA reading --------------------------------------------------------
def read_fasta(fasta_file):
    seqs = {}
    with open(fasta_file) as f:
        cur, buf = None, []
        for line in f:
            line = line.rstrip()
            if not line: continue
            if line.startswith('>'):
                if cur: seqs[cur] = ''.join(buf)
                cur = line[1:].split()[0]; buf = []
            else:
                buf.append(line)
        if cur: seqs[cur] = ''.join(buf)
    return seqs

# -- GEESE parsing --------------------------------------------------------
def parse_geese(geese_file):
    """
    Returns:
      occ: atom -> list of (genome, start, end, strand)
      genome_order: genome -> list of (atom, strand)
      raw_lines: all lines for final rewriting
    """
    occ = defaultdict(list)
    genome_order = defaultdict(list)
    raw_lines = []
    with open(geese_file) as f:
        for L in f:
            raw_lines.append(L.rstrip('\n'))
            if L.startswith('#'): continue
            p = L.rstrip('\n').split('\t')
            if len(p) < 6: continue
            genome, _, atom, strand = p[0], p[1], p[2], p[3]
            try:
                start, end = int(p[4]), int(p[5])
            except ValueError:
                continue
            occ[atom].append((genome, start, end, strand))
            genome_order[genome].append((atom, strand))
    return occ, genome_order, raw_lines

# -- Unique contexts ------------------------------------------------------
def compute_unique_contexts(genome_order):
    """
    For each atom, collect all (prev, next) pairs and count
    how many 'unique' ones after ignoring repeats on either side.
    """
    raw = defaultdict(set)
    for order in genome_order.values():
        for i,(atom,_) in enumerate(order):
            prev_a = order[i-1][0] if i>0 else None
            next_a = order[i+1][0] if i<len(order)-1 else None
            raw[atom].add((prev_a, next_a))
    uniq = {}
    for atom, pairs in raw.items():
        seen_p, seen_n = set(), set()
        cnt = 0
        for p,n in pairs:
            if p not in seen_p and n not in seen_n:
                cnt += 1
                seen_p.add(p); seen_n.add(n)
        uniq[atom] = cnt
    return uniq

# -- Build segment metadata ----------------------------------------------
def revcomp(s):
    return s.translate(str.maketrans('ACGTacgt','TGCAtgca'))[::-1]

def build_segments(occ, seqs):
    """
    Build a GfaSegment per atom:
      • depth = total #occurrences
      • duplicated = does any genome contain >1?
      • seq/length from first occurrence
    """
    segs = []
    for atom, hits in occ.items():
        depth = len(hits)
        by_genome = defaultdict(int)
        for g,_,_,_ in hits:
            by_genome[g] += 1
        dup = any(c>1 for c in by_genome.values())
        g0,s0,e0,st0 = hits[0]
        length = e0 - s0
        seq = ""
        if g0 in seqs:
            frag = seqs[g0][s0:e0]
            seq = revcomp(frag) if st0=='-' else frag
        segs.append(GfaSegment(atom, seq, depth, length, dup, 0))
    return segs

# -- Global filtering -----------------------------------------------------
def filter_atoms(segs, min_depth, max_length, remove_dup, max_unique):
    """
    Which atoms to remove *globally*?
    """
    to_rm = set()
    for s in segs:
        if((remove_dup and s.duplicated)):
            to_rm.add(s.name)
        elif s.length <= max_length:
            if (((s.depth < min_depth) and (s.unique_context > max_unique))):
                to_rm.add(s.name)
    return to_rm

# -- Pair scoring w/ context tracking ------------------------------------
def compute_in_out(genome_order):
    IN, OUT = defaultdict(set), defaultdict(set)
    for order in genome_order.values():
        for i,(a,st) in enumerate(order):
            # respect strand
            if st == '-':
                if i>0:  OUT[a].add(order[i-1][0])
                if i<len(order)-1: IN[a].add(order[i+1][0])
            else:
                if i>0:  IN[a].add(order[i-1][0])
                if i<len(order)-1: OUT[a].add(order[i+1][0])
    return IN, OUT

def find_high_diverse_pairs_with_contexts(genome_order, atom_lengths,
                                          min_in, min_out, max_span):
    """
    For each genome, for each i<j:
      a=order[i], b=order[j]
      require |OUT[a]|>=min_in, |IN[b]|>=min_out, span<=max_span
      record the in-between tuple and which genome it came from.

    Returns:
      contexts_map: (a,b) -> dict { between_tuple : set(of genomes) }
    """
    IN, OUT = compute_in_out(genome_order)
    contexts_map = defaultdict(lambda: defaultdict(set))
    for genome, order in genome_order.items():
        n = len(order)
        # prefix‐sum of lengths
        ps = [0]*(n+1)
        for i,(atom,_) in enumerate(order):
            ps[i+1] = ps[i] + atom_lengths.get(atom,0)
        for i,(a,st) in enumerate(order):
            if len(OUT[a]) < min_in: continue
            for j in range(i+1, n):
                b,_ = order[j]
                if len(IN[b]) < min_out: continue
                span = ps[j+1] - ps[i]
                if span > max_span: break
                between = tuple(order[k][0] for k in range(i+1, j))
                contexts_map[(a,b)][between].add(genome)
    return contexts_map

# -- New: context‐based filtering ----------------------------------------
def context_filter(contexts_map, genome_order):
    """
    For each (a,b) in contexts_map:
      • union_atoms = all atoms appearing in any 'between' tuple
      For each such atom:
        1) Recompute atom_ctx2genomes from genome_order: for every occurrence
           of `atom`, record its (prev,next) context (respecting strand) → set(genomes).
        2) Choose best_ctx = the context with the largest genome set.
        3) In *all* other contexts, for each genome in that genome‐set, mark atom for removal.
    Finally apply all per-genome removals to genome_order in-place.
    Returns per_genome_rm: genome → set(atoms_to_remove).
    """
    per_genome_rm = defaultdict(set)

    for (a, b), between_ctxs in contexts_map.items():
        # union of all the 'between' atoms for this (a,b) pair
        union_atoms = set(x for ctx in between_ctxs for x in ctx)

        for atom in union_atoms:
            # build that atom's true contexts across *all* genome_order
            atom_ctx2gen = defaultdict(set)
            for genome, order in genome_order.items():
                for idx, (at, strand) in enumerate(order):
                    if at != atom:
                        continue
                    prev_atom = order[idx-1][0] if idx > 0 else None
                    next_atom = order[idx+1][0] if idx < len(order)-1 else None
                    # respect strand
                    if strand == '+':
                        ctx = (prev_atom, next_atom)
                    else:
                        ctx = (next_atom, prev_atom)
                    atom_ctx2gen[ctx].add(genome)

            if not atom_ctx2gen:
                continue

            # pick the single best context (max number of genomes)
            best_ctx, best_genomes = max(
                atom_ctx2gen.items(),
                key=lambda item: len(item[1])
            )

            # remove `atom` from all *other* contexts' genomes
            for ctx, genomes in atom_ctx2gen.items():
                if ctx == best_ctx:
                    continue
                for g in genomes:
                    per_genome_rm[g].add(atom)

    # apply all per-genome removals
    for genome, lst in genome_order.items():
        genome_order[genome] = [
            (at, st) for (at, st) in lst
            if at not in per_genome_rm.get(genome, ())
        ]

    return per_genome_rm

# -- Final .geese rewrite -------------------------------------------------
def rewrite_final_geese(raw_lines, global_rm, per_genome_rm, out_file):
    with open(out_file, 'w') as f:
        for L in raw_lines:
            if L.startswith('#'):
                f.write(L + "\n")
                continue
            p = L.split('\t')
            if len(p) < 3:
                f.write(L + "\n")
                continue
            genome, atom = p[0], p[2]
            if atom in global_rm or atom in per_genome_rm.get(genome, set()):
                continue
            f.write(L + "\n")

# -- main -----------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description="filter_geese.py with context‐based atom removal")
    p.add_argument("geese_in")
    p.add_argument("fasta")
    p.add_argument("-o", "--output", required=True,
                   help="Final filtered .geese")
    p.add_argument("--min-depth",  type=int, default=10)
    p.add_argument("--max-length", type=int, default=10000)
    p.add_argument("--remove-dup",  action="store_true")
    p.add_argument("--max-unique", type=int, default=1)
    p.add_argument("--pair-min-in",   type=int, default=2)
    p.add_argument("--pair-min-out",  type=int, default=2)
    p.add_argument("--pair-max-span", type=int, default=70000)
    p.add_argument("--iterations",   type=int, default=1, help="Repeat global filter this many times, recomputing depths/contexts")
    args = p.parse_args()

 # 1) initial load
    seqs = read_fasta(args.fasta)
    occ, genome_order, raw = parse_geese(args.geese_in)

    global_rm = set()
    per_genome_rm = {}

    for it in range(args.iterations):
        logging.info(f"=== ITERATION {it+1} ===")

        # 2) rebuild segment metadata & unique_context
        #    (we only need occ → segs → unique_context; sequence coords not needed
        #     for repeated passes, so we reuse occ with dummy coords)
        segs = build_segments(occ, seqs)
        uniq = compute_unique_contexts(genome_order)
        for s in segs:
            s.unique_context = uniq.get(s.name, 0)

        # 3) global filter
        new_global_rm = filter_atoms(
            segs,
            args.min_depth,
            args.max_length,
            args.remove_dup,
            args.max_unique
        )
        if not new_global_rm:
            logging.info("No more atoms to remove globally; stopping iterations.")
            break

        logging.info(f"Globally removing {len(new_global_rm)} atoms")
        global_rm |= new_global_rm

        # apply global removal to genome_order
        for g, lst in genome_order.items():
            genome_order[g] = [(a,st) for a,st in lst if a not in new_global_rm]

        # rebuild a minimal occ for next round (coords unused)
        occ = defaultdict(list)
        for g, lst in genome_order.items():
            for a,st in lst:
                occ[a].append((g,0,0,st))

    # 5) pair scoring & context-based removals (unchanged)
    atom_lengths = {s.name: s.length for s in build_segments(parse_geese(args.geese_in)[0], seqs)
                    if s.name not in global_rm}
    
    contexts_map = find_high_diverse_pairs_with_contexts(
        genome_order, atom_lengths,
        args.pair_min_in, args.pair_min_out, args.pair_max_span
    )
    logging.info(f"Found {len(contexts_map)} candidate pairs")
    per_genome_rm = context_filter(contexts_map, genome_order)
    logging.info(f"Context-based removal: atoms removed in {sum(len(v) for v in per_genome_rm.values())} genomes")

    # 7) final write
    rewrite_final_geese(raw, global_rm, per_genome_rm, args.output)
    logging.info(f"Wrote filtered .geese → {args.output}")

if __name__=="__main__":
    main()
