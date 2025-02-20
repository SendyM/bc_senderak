# **bc_senderak**
## Dokumentácia pre skripty

Všetky skripty vyžadujú **Python 3** a väčšina z nich používa knižnicu **NetworkX**.

---

## **1. gfa_vertex_table.py**
### **Popis**  
Tento skript načíta **GFA súbor**, vytvorí **neorientovaný graf** (používajúci iba jedinečné hrany) a pre každý uzol (segment) vypočíta jeho:
- **Stupeň (degree)**
- **Betweenness centrality**  

Následne vytvorí **CSV tabuľku (s oddelením tabulátorom)** obsahujúcu segmenty a vypočítané metriky.  
Tabuľku možno zoradiť podľa stupňa alebo centrality.

### **Kľúčové parametre**  
- `gfa_file` – Cesta k vstupnému GFA súboru.  
- `--output-file` – Výstupný súbor (ak nie je zadané, výstup sa vypíše na štandardný výstup).  
- `--sort-by` – Stĺpec na zoradenie tabuľky (`degree` alebo `centrality`).  
- `--descending` – Ak je nastavené, zoradí tabuľku zostupne.  

---

## **2. remove_segments.py**
### **Popis**  
Tento skript odstráni z pôvodného **GFA súboru** špecifikované segmenty a všetky súvisiace hrany.  
Výsledkom je nový GFA súbor bez daných segmentov.

### **Kľúčové parametre**  
- `segments_file` – Textový súbor s menami segmentov na odstránenie (jeden názov na riadok).  
- `original_gfa` – Cesta k pôvodnému GFA súboru.  
- `output_gfa` – Cesta k novému GFA súboru.  

---

## **3. plot_distribution.py**
### **Účel skriptu**  
Skript vizualizuje **degree a betweenness centrality** z tabuľky so stĺpcami:  
- `segment`  
- `degree`  
- `centrality`  

Možno obmedziť hodnoty v histograme a boxplotoch (napr. `--degree-min 5 --degree-max 100`).

### **Vstupné dáta**  
Tabuľkový súbor (`CSV/TSV`) s tromi základnými stĺpcami:
- **segment** – názov segmentu/uzla  
- **degree** – počet susedov/unikátnych spojení  
- **centrality** – betweenness centrality  

### **Voliteľné argumenty**  
- `--delimiter` – Oddeľovač v súbore (predvolene `\t`).  
- `--bins` – Počet binov v histograme (predvolene `50`).  
- `--degree-min`, `--degree-max` – Minimálna/maximálna hodnota stupňa.  
- `--centrality-min`, `--centrality-max` – Minimálna/maximálna hodnota centrality.  
- `--output-distributions` – Výstupný súbor s **2×2 podgrafmi** (histogram + KDE a boxplot pre **degree a centrality**).  
- `--output-scatter` – Výstupný **scatter plot** (*degree vs. centrality*).  

### **Výstup**  
Obrázok `.png` obsahujúci **4 podgrafy**:
1. Histogram + KDE pre **degree**
2. Boxplot pre **degree**
3. Histogram + KDE pre **centrality**
4. Boxplot pre **centrality**  

**Scatter plot (.png)** znázorňuje vzťah medzi **degree a betweenness centrality**.

---

## **4. iterative_flag_te_candidates_with_geese.py**
### **Účel skriptu**  
Tento skript slúži na:
- **Iteratívne odstraňovanie kandidátnych uzlov z GFA grafu**.
- **Rekonštrukciu spojení (bridging)**.

### **Vstupné súbory**  
- **GFA súbor**  
- **.geese súbor**  

### **Logovanie**  
Skript vypíše štatistiky o odstránených uzloch:
- Počet uzlov odstránených na základe **prahu stupňa**.
- Počet uzlov odstránených na základe **prahu betweenness centrality**.
- Počet **nútených odstránení**.

### **Príklad spustenia**  
```sh
./iterative_flag_te_candidates_with_geese.py original.gfa simplified.gfa \
  --geese-file ecol-train-1000.geese \
  --iterations 5 \
  --degree-threshold 10 \
  --bc-threshold 0.05 \
  --flagged-output-file flagged_details.txt \
  --plain-output-file flagged_names.txt
```

---

## **5. score_gfa.py**
### **Účel**  
Skript ohodnotí zostavový **GFA graf** podľa:  
1. **Priemerného stupňa (average degree)**  
2. **Pomeru zostávajúcich segmentov (fraction remaining)**  

### **Výpočet skóre**  
Vzorec:
```math
Score = (avg_degree + 1) / (fraction_remaining^0.5)
```
Kde:  
- **avg_degree** = 2 × počet hrán / počet uzlov  
- **fraction_remaining** = aktuálny počet segmentov / pôvodný počet segmentov (`--original-nodes`).  

**Nižšie skóre znamená, že graf je menej "rozpletený" a zachováva viac segmentov.**

### **Vstupné dáta**  
- `gfa_file`  
- `--original-nodes` (Povinné) – Udáva pôvodný počet segmentov.  
- `--verbose` (voliteľné) – Zapne podrobné výstupy.  

### **Výstup**  
**Konečné skóre grafu** (číselná hodnota).  
Pri `--verbose` aj podrobné metriky.

### **Príklad použitia**  
```sh
./score_gfa.py simplified.gfa --original-nodes 4262 --verbose
```
