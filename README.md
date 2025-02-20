# bc_senderak
Dokumentácia pre skripty
Všetky skripty vyžadujú Python 3 a väčšina z nich používa knižnicu NetworkX.

1. gfa_vertex_table.py
Popis:
Tento skript načíta GFA súbor, vytvorí neorientovaný graf (používajúci iba jedinečné hrany) a pre každý uzol (segment) vypočíta jeho stupeň (degree) a betweenness centrality. Následne vytvorí tabuľku (vo formáte CSV s oddelením tabulátorom), ktorá obsahuje zoznam segmentov spolu s vypočítanými metrikami. Výstup možno zoradiť podľa stupňa alebo centrality.
Kľúčové parametre:
gfa_file: Cesta k vstupnému GFA súboru.
--output-file: Súbor, do ktorého sa uloží výstupná tabuľka (ak nie je zadané, výstup sa vypíše na štandardný výstup).
--sort-by: Stĺpec, podľa ktorého sa má tabuľka zoradiť (možnosti: degree alebo centrality).
--descending: Ak je nastavené, zoradí tabuľku zostupne.
2. remove_segments.py
Popis:
Tento skript odstráni z pôvodného GFA súboru špecifikované segmenty (a všetky hrany, ktoré s nimi súvisia). Výsledkom je nový GFA súbor, v ktorom už nie sú zahrnuté dané segmenty ani žiadne spojenia, ktoré na ne odkazujú.
Kľúčové parametre:
segments_file: Textový súbor, ktorý obsahuje zoznam mien segmentov na odstránenie (jeden názov na riadok).
original_gfa: Cesta k pôvodnému GFA súboru.
output_gfa: Cesta k novému GFA súboru, v ktorom budú dané segmenty odstránené.

3. plot_distribution.py
Účel skriptu
Skript slúži na vizualizáciu dvoch metrík (degree a betweenness centrality) z tabuľky, ktorá obsahuje stĺpce:
segment
degree
centrality
Pomocou argumentov je možné obmedziť hodnoty, ktoré sa zobrazia v histogramoch a boxplotoch (napr. --degree-min 5 --degree-max 100).
Vstup
Tabuľkový súbor (CSV/TSV) s tromi základnými stĺpcami:
segment (názov segmentu/uzla)
degree (počet susedov/unikátnych spojení)
centrality (betweenness centrality)
Voliteľné argumenty:
--delimiter: Špecifikuje oddeľovač v súbore (predvolene \t).
--bins: Počet binov v histograme (predvolene 50).
--degree-min a --degree-max: Minimálna a maximálna hodnota stupňa, ktorú chcete zobraziť.
--centrality-min a --centrality-max: Minimálna a maximálna hodnota betweenness centrality, ktorú chcete zobraziť.
--output-distributions: Názov výstupného súboru s 2×2 podgrafmi (histogram + KDE a boxplot pre degree a centrality).
--output-scatter: Názov výstupného súboru so scatter plotom (degree vs. centrality).
Výstup
Obrázok (.png alebo iný formát), ktorý obsahuje 4 podgrafy:
Histogram + KDE pre degree (obmedzený na rozsah --degree-min až --degree-max)
Boxplot pre degree (s rovnakým rozsahom na osi x)
Histogram + KDE pre betweenness centrality (obmedzený na --centrality-min až --centrality-max)
Boxplot pre betweenness centrality
Scatter plot (.png), ktorý znázorňuje vzťah medzi stupňom (degree) a betweenness centrality v zadanom rozsahu.

4. iterative_flag_te_candidates_with_geese.py
Účel skriptu
Tento skript slúži na:
Iteratívne odstraňovanie kandidátnych uzlov z GFA grafu.
Rekonštrukciu spojení (bridging)
Vstupné súbory
GFA súbor
.geese súbor
Logovanie:
Skript vypíše štatistiky o odstránených uzloch – počet uzlov odstránených na základe prahu stupňa, prahu betweenness centrality a počet nútených odstránení.
Príklad spustenia
./iterative_flag_te_candidates_with_geese.py original.gfa simplified.gfa \
  --geese-file ecol-train-1000.geese \
  --iterations 5 \
  --degree-threshold 10 \
  --bc-threshold 0.05 \
  --flagged-output-file flagged_details.txt \
  --plain-output-file flagged_names.txt

Tento príkaz:
Načíta original.gfa a ecol-train-1000.geese.
Vykoná 5 iterácií odstraňovania uzlov s prahmi 10 (degree) a 0.05 (betweenness centrality).
Zapíše výsledný zjednodušený GFA súbor do simplified.gfa.
Vytvorí súbor flagged_details.txt s podrobnými informáciami o odstránených uzloch a flagged_names.txt so zoznamom názvov odstránených uzlov.

5. score_gfa.py
Účel
Tento skript slúži na ohodnotenie zostavového grafu vo formáte GFA podľa toho, ako je "rozpletený" (untangled) a koľko segmentov zostáva. Hodnotí sa pomocou dvoch metrík:
Priemerný stupeň (average degree)
Pomer zostávajúcich segmentov (fraction remaining)
Výpočet skóre
Skóre sa počíta podľa vzorca:
  Score = (avg_degree + 1) / (fraction_remaining^0.5)
Kde:
avg_degree – priemerný stupeň grafu (počítaný ako 2 × počet hrán / počet uzlov).
fraction_remaining – pomer aktuálneho počtu segmentov ku pôvodnému počtu segmentov (parameter --original-nodes).
Nižšie skóre znamená, že graf je menej "rozpletený" a zachováva viac segmentov.
Vstupy
gfa_file
--original-nodes: Povinný parameter, ktorý udáva pôvodný počet segmentov.
--verbose: (Voliteľné) Prepínač, ktorý zapne podrobné výstupy.
Výstup
Skript vypíše:
Konečné skóre grafu (číselná hodnota).
Ak je zapnutý prepínač --verbose, skript tiež zobrazí podrobné metriky.
Príklad použitia
./score_gfa.py simplified.gfa --original-nodes 4262 --verbose
