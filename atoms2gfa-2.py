#!/usr/bin/env python3
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(
        description='Create GFA from Geese output, excluding certain classes, and writing result to a .gfa file.'
    )
    parser.add_argument('atoms', help='Input file with atoms from Geese')
    parser.add_argument('fasta', help='Input file with FASTA sequence')
    parser.add_argument(
        '-c', '--cutoff',
        type=int,
        default=10000000,
        help='Exclude atoms whose class appears more times than this cutoff (default=10000000).'
    )
    parser.add_argument(
        '-o', '--output',
        default='out.gfa',
        help='Output file name for the GFA (default: out.gfa).'
    )
    parser.add_argument(
        '-x', '--exclude',
        help='File listing classes to exclude, one per line.'
    )
    args = parser.parse_args()

    # 1) Read excluded classes into a set
    excluded_classes = set()
    if args.exclude:
        with open(args.exclude, 'r') as excl:
            for line in excl:
                cls = line.strip()
                if cls:
                    excluded_classes.add(cls)

    # 2) Prepare to read FASTA
    fachrom = {}
    faplas = {}
    sequences = read_fasta(args.fasta, fachrom, faplas)

    # 3) Read the atom definitions
    segments = read_atoms(args.atoms)

    # 4) Count occurrences of each class (for cutoff filtering)
    atomstats = {}
    # Track if class appears in chromosome or plasmid
    atchrom = {}
    atplas = {}

    for seg in segments:
        cls = seg['class']
        atomstats[cls] = atomstats.get(cls, 0) + 1
        if cls not in atchrom:
            atchrom[cls] = False
        if cls not in atplas:
            atplas[cls] = False

        # Mark whether the sequence for this segment is chromosome or plasmid
        if fachrom.get(seg['name'], False):
            atchrom[cls] = True
        if faplas.get(seg['name'], False):
            atplas[cls] = True

    # 5) Build a histogram of how many classes appear X times
    hist = []
    for cls, count in atomstats.items():
        while len(hist) <= count:
            hist.append(0)
        hist[count] += 1

    # Print histogram to stderr
    print("Histogram of atom counts per class:", file=sys.stderr)
    print(hist, file=sys.stderr)

    # 6) We'll keep track of which classes we've actually output as segments
    used_classes = set()

    # For bridging: remember the last included atom *per sequence* 
    # so we can link from it to the next included atom on the same seq.
    last_included_per_seq = {}

    # 7) Write GFA to file
    with open(args.output, 'w') as gfa_file:
        for atom in segments:
            cls = atom['class']
            seqname = atom['name']

            # Skip if:
            # - the class is over the cutoff, or
            # - the class is in the exclude list
            if atomstats[cls] > args.cutoff or cls in excluded_classes:
                # We do NOT update last_included_per_seq[seqname].
                continue

            # This atom is "included" in the GFA, so we see if we need 
            # to create an 'S' line for its class
            if cls not in used_classes:
                used_classes.add(cls)
                subseq = sequences[seqname][atom['start']:atom['end']]
                if atom['strand'] == '-':
                    subseq = reverse_complement(subseq)

                color_tags = ''
                if atchrom[cls]:
                    if atplas[cls]:
                        color_tags = '\tCL:z:#aaaa00\tC2:z:#aaaa00'
                    else:
                        color_tags = '\tCL:z:#00aa00\tC2:z:#00aa00'
                elif atplas[cls]:
                    color_tags = '\tCL:z:#aa0000\tC2:z:#aa0000'

                # Print the segment line: S <id> <sequence> [tags]
                gfa_file.write(f"S\t{cls}\t{subseq}{color_tags}\n")

            # If there was a previous included atom in the same seq, 
            # link that atom to this atom
            if seqname in last_included_per_seq:
                prev_atom = last_included_per_seq[seqname]
                # L <class1> <strand1> <class2> <strand2> 0M
                gfa_file.write("\t".join([
                    'L',
                    prev_atom["class"], prev_atom['strand'],
                    cls, atom['strand'],
                    '0M'
                ]) + "\n")

            # Update the last included atom for this sequence
            last_included_per_seq[seqname] = atom


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    seq = seq.upper()
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def parse_atom_line(line):
    """
    Parse a single line from the atoms file.
    Expected format: 
      #name  atom_nr  class  strand  start  end
      a1     1        1      +       0      129517
    """
    parts = line.split()
    return {
        'name': parts[0],
        'class': parts[2],
        'strand': parts[3],
        'start': int(parts[4]),
        'end': int(parts[5])
    }


def read_fasta(fasta_file, ischrom, isplas):
    """
    Read sequences from a FASTA file into a dict: { seq_name: sequence_string }
    Also fill ischrom[name], isplas[name] with booleans based on the header line.
    """
    sequences = {}
    current_seq_lines = []
    current_name = None

    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Store the previous sequence
                if current_name is not None:
                    sequences[current_name] = ''.join(current_seq_lines)

                current_seq_lines = []
                # Parse the header
                parts = line.split()
                current_name = parts[0][1:]  # remove '>'

                if 'chromosome=true' in line:
                    ischrom[current_name] = True
                else:
                    ischrom[current_name] = False

                if 'plasmid=true' in line:
                    isplas[current_name] = True
                else:
                    isplas[current_name] = False

            else:
                current_seq_lines.append(line.strip())

        # Last sequence
        if current_name is not None:
            sequences[current_name] = ''.join(current_seq_lines)

    return sequences


def read_atoms(atoms_file):
    """
    Read lines from an atoms file, ignoring #comment lines.
    Return a list of dicts (one per atom).
    """
    segments = []
    with open(atoms_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            atom = parse_atom_line(line)
            segments.append(atom)
    return segments


if __name__ == '__main__':
    main()

