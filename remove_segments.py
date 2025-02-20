#!/usr/bin/env python3
import argparse

def load_segments_to_remove(segments_file):
    """
    Load the segment names (one per line) from the given file.
    Returns a set of segment names.
    """
    segments_to_remove = set()
    with open(segments_file, 'r') as f:
        for line in f:
            seg = line.strip()
            if seg:
                segments_to_remove.add(seg)
    return segments_to_remove

def process_gfa(original_gfa, removed_segments):
    """
    Reads the original GFA file and returns a list of lines to be kept.
    
    - S lines: Only kept if the segment name is NOT in removed_segments.
    - L lines: Only kept if neither endpoint is in removed_segments.
    - Other lines: Kept unchanged.
    """
    kept_lines = []
    with open(original_gfa, 'r') as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("S"):
                # S line format: S <segment> <sequence> [other tags...]
                fields = line.split("\t")
                if len(fields) >= 2:
                    seg = fields[1]
                    if seg in removed_segments:
                        continue  # skip this S line
                kept_lines.append(line)
            elif line.startswith("L"):
                # L line format: L <seg1> <ori1> <seg2> <ori2> <overlap> [other tags...]
                fields = line.split("\t")
                if len(fields) >= 6:
                    seg1 = fields[1]
                    seg2 = fields[3]
                    if seg1 in removed_segments or seg2 in removed_segments:
                        continue  # skip this L line
                kept_lines.append(line)
            else:
                # Keep any other lines (headers, comments, etc.)
                kept_lines.append(line)
    return kept_lines

def main():
    parser = argparse.ArgumentParser(
        description="Remove segments (and corresponding edges) from a GFA file given a list of segments."
    )
    parser.add_argument("segments_file", help="Path to the file containing segment names to remove (one per line).")
    parser.add_argument("original_gfa", help="Path to the original GFA file.")
    parser.add_argument("output_gfa", help="Path to the output GFA file (with the segments removed).")
    args = parser.parse_args()

    # Load segments to remove.
    removed_segments = load_segments_to_remove(args.segments_file)
    print(f"Loaded {len(removed_segments)} segments to remove.")

    # Process the GFA file.
    kept_lines = process_gfa(args.original_gfa, removed_segments)
    print(f"Kept {len(kept_lines)} lines out of the original GFA.")

    # Write the output GFA file.
    with open(args.output_gfa, 'w') as out:
        for line in kept_lines:
            out.write(line + "\n")
    print(f"Output GFA file written to '{args.output_gfa}'.")

if __name__ == "__main__":
    main()

