#!/usr/bin/env python3
import sys

def parse_fasta(fasta_file):
    """
    Reads the FASTA file and returns a dictionary mapping sequence IDs
    to their type ('chromosome', 'plasmid', or None if not marked).
    """
    seq_types = {}
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                line = line.strip()
                # Remove the leading '>' and any surrounding punctuation.
                header = line[1:]
                tokens = header.split()
                seq_id = tokens[0]
                # Check for explicit markers in the header.
                if 'chromosome=true' in header:
                    seq_types[seq_id] = 'chromosome'
                elif 'plasmid=true' in header:
                    seq_types[seq_id] = 'plasmid'
                else:
                    seq_types[seq_id] = None
    return seq_types

def parse_geese(geese_file, seq_types):
    """
    Reads the .geese file and returns a dictionary mapping atom IDs (from the 'class' field)
    to a set of sequence types (chromosome, plasmid) in which the atom appears.
    """
    atom_dict = {}
    with open(geese_file, 'r') as f:
        header = f.readline()  # skip the header line
        for line in f:
            if not line.strip():
                continue
            fields = line.strip().split()
            # The expected columns are:
            # name, atom_nr, class, strand, start, end
            seq_name = fields[0]
            atom_id = fields[2]  # using the 'class' field as the atom identifier
            # Look up the type for this sequence id from the FASTA file
            seq_type = seq_types.get(seq_name)
            if seq_type:  # Only consider explicitly marked sequences
                if atom_id not in atom_dict:
                    atom_dict[atom_id] = set()
                atom_dict[atom_id].add(seq_type)
    return atom_dict

def main(fasta_file, geese_file, output_file):
    seq_types = parse_fasta(fasta_file)
    atom_dict = parse_geese(geese_file, seq_types)
    # Identify atoms found in both plasmid and chromosome sequences.
    with open(output_file, 'w') as out:
        for atom_id, types in atom_dict.items():
            if 'chromosome' in types and 'plasmid' in types:
                out.write(atom_id + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <fasta_file> <geese_file> <output_file>")
        sys.exit(1)
    fasta_file = sys.argv[1]
    geese_file = sys.argv[2]
    output_file = sys.argv[3]
    main(fasta_file, geese_file, output_file)

