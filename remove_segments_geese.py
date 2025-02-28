#!/usr/bin/env python3
import argparse
import csv

def load_segments_to_remove(segments_file):
    """
    Načíta názvy segmentov (jeden názov na riadok) zo zadaného súboru.
    Vráti množinu názvov segmentov.
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
    Načíta pôvodný GFA súbor a vyfiltruje S riadky (segmenty),
    ktoré nie sú v množine removed_segments. L riadky sa ignorujú, 
    pretože budú neskôr nahradené novými spojeniami založenými na .geese súbore.
    Ďalšie riadky (hlavičky, komentáre) sú ponechané.
    
    Vráti trojicu: (filtered S riadky, other riadky)
    """
    s_lines = []
    other_lines = []
    with open(original_gfa, 'r') as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("S"):
                fields = line.split("\t")
                if len(fields) >= 2:
                    seg = fields[1]
                    if seg in removed_segments:
                        continue  # preskočíme tento segment
                s_lines.append(line)
            elif line.startswith("L"):
                # L riadky budú ignorované (nahradíme ich novými)
                continue
            else:
                other_lines.append(line)
    return s_lines, other_lines

def build_links_from_geese(geese_file, removed_segments, valid_segments):
    """
    Načíta .geese súbor a pre každú vzorku (sample) zoradí záznamy podľa atom_nr.
    Záznamy, ktorých hodnota 'class' je v removed_segments alebo ktoré nie sú v valid_segments,
    sa ignorujú. Pre každú vzorku skript prejde zoradené záznamy a pre každú dvojicu susediacich
    segmentov vytvorí nový L riadok vo formáte:
       L <seg1> <ori1> <seg2> <ori2> 0M

    Vráti zoznam nových L riadkov.
    """
    # Dictionary: sample -> list of tuples (atom_nr, class, strand)
    sample_order = {}
    with open(geese_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        # Očakávané názvy stĺpcov: "#name" (alebo "name"), "atom_nr", "class", "strand", "start", "end"
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
            # Ignorujeme záznamy, ktoré nie sú v aktuálnom GFA (valid_segments)
            if seg_class not in valid_segments:
                continue
            if sample not in sample_order:
                sample_order[sample] = []
            sample_order[sample].append((atom_nr, seg_class, strand))

    new_l_lines_set = set()
    # Pre každú vzorku zoradíme záznamy podľa atom_nr a vytvoríme spojenia medzi susediacimi segmentmi
    for sample, entries in sample_order.items():
        sorted_entries = sorted(entries, key=lambda x: x[0])
        # Vyfiltrujeme záznamy, ktoré patria do removed_segments
        filtered_entries = [
            (atom_nr, seg_class, strand)
            for (atom_nr, seg_class, strand) in sorted_entries
            if seg_class not in removed_segments
        ]
        for i in range(len(filtered_entries) - 1):
            seg1 = filtered_entries[i][1]
            ori1 = filtered_entries[i][2]
            seg2 = filtered_entries[i+1][1]
            ori2 = filtered_entries[i+1][2]
            new_line = f"L\t{seg1}\t{ori1}\t{seg2}\t{ori2}\t0M"
            new_l_lines_set.add(new_line)
    return list(new_l_lines_set)

def main():
    parser = argparse.ArgumentParser(
        description="Odstráni segmenty zo súboru GFA podľa zadaného zoznamu a následne znovu prepojí zostávajúce segmenty "
                    "na základe poradia zo súboru .geese."
    )
    parser.add_argument("segments_file", help="Cesta k súboru s názvami segmentov na odstránenie (jeden názov na riadok).")
    parser.add_argument("original_gfa", help="Cesta k pôvodnému GFA súboru.")
    parser.add_argument("output_gfa", help="Cesta k výstupnému GFA súboru (s odstránenými segmentmi a novými spojeniami).")
    parser.add_argument("--geese-file", required=True,
                        help="Cesta k .geese súboru s poradím segmentov (hlavička: #name, atom_nr, class, strand, start, end).")
    args = parser.parse_args()

    # Načíta zoznam segmentov na odstránenie.
    removed_segments = load_segments_to_remove(args.segments_file)
    print(f"Nájdených {len(removed_segments)} segmentov na odstránenie.")

    # Načíta pôvodný GFA súbor a vyfiltruje S riadky.
    s_lines, other_lines = process_gfa(args.original_gfa, removed_segments)
    print(f"Z pôvodného GFA zostalo {len(s_lines)} segmentov (S riadkov).")

    # Získame množinu platných segmentov (tých, ktoré ostali v GFA).
    valid_segments = set()
    for line in s_lines:
        fields = line.split("\t")
        if len(fields) >= 2:
            valid_segments.add(fields[1])
    
    # Vytvoríme nové L riadky (spojenia) na základe poradia zo súboru .geese.
    new_l_lines = build_links_from_geese(args.geese_file, removed_segments, valid_segments)
    print(f"Vytvorených {len(new_l_lines)} nových L riadkov (spojení) na základe .geese súboru.")

    # Zapíšeme výsledný GFA súbor: pôvodné ostatné riadky, vyfiltrované S riadky a nové L riadky.
    with open(args.output_gfa, 'w') as out:
        for line in other_lines:
            out.write(line + "\n")
        for line in s_lines:
            out.write(line + "\n")
        for line in new_l_lines:
            out.write(line + "\n")
    print(f"Výstupný GFA súbor bol zapísaný do '{args.output_gfa}'.")

if __name__ == "__main__":
    main()

