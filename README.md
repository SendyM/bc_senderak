# **bc_senderak**
# Pangenomic Graph Filtering and Visualization – README

Tento projekt obsahuje skripty na filtrovanie pangenomických grafov vytvorených z výstupu nástroja GEESE a ich následnú vizualizáciu pomocou D3.js alebo Bandage.

## Prehľad skriptov

### 1. `iterative_geese_filter.py`
Hlavný skript pre **filtrovanie** `.geese` súboru:
- kombinuje globálne pravidlá (minimálna hĺbka, dĺžka, duplikáty, unikátne kontexty),
- iteruje viacnásobne podľa nastavení,
- následne aplikuje **kontextovo založené čistenie** atómov na základe ich výskytu v rôznych genómoch.

**Výstup:** nový `.geese` súbor s odfiltrovanými atómami.

---

### 2. `atoms2gfa.py`
Konvertuje `.geese` a FASTA súbor na GFA graf.

**Výstup:** `.gfa` súbor so segmentmi a hranami pripravený na ďalšie spracovanie.

---

### 3. `score_graph.py`
Porovnáva pôvodný a filtrovaný GFA graf a vypočíta **HI (Hairball Index)**

**Výstup:** textový výstup so štatistikami a skóre (nižšie = lepšie).

---

### 4. `gfa2json.py`
Prevádza GFA + `.geese` súbor do formátu JSON vhodného na vizualizáciu:
- každý uzol obsahuje počet genómov, v ktorých sa vyskytuje,
- farba uzla sa prenáša z GFA.

**Výstup:** `graph.json` súbor pre vizualizáciu.

---

### 5. `graph_D3.html`
HTML stránka využívajúca D3.js na **interaktívnu vizualizáciu grafu**:
- úprava veľkosti uzlov a hrán,
- možnosť zapnúť heatmapu podľa výskytu uzlov (usage).

---

## Ako reprodukovať finálny graf

1. Spusti filtrovanie `.geese` súboru:

```bash
python3 iterative_geese_filter.py ecol-train-1000.geese ecol-train.fa \
 -o filtered.geese \
 --min-depth 25 --remove-dup --iterations 4
```
2. Vytvor GFA z filtrovaného `.geese`:

```bash
python3 atoms2gfa.py filtered.geese ecol-train.fa -o filtered.gfa
```
3. Vyzualizácia pomocou nástroja Bandage:
```bash
Bandage load filtered.gfa
```

3. Vyzualizácia pomocou D3.js:
```bash
python3 gfa2json.py filtered.gfa filtered.geese graph.json

#Spusti lokálny HTTP server a otvor vizualizáciu:
python3 -m http.server 8000
http://localhost:8000/graph_D3.html
