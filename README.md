# make_astral_map 

Generates a **map file** (TSV `taxon<TAB>group`) for **ASTRAL** from hundreds (or thousands) of **NEXUS** and **FASTA** alignments. No external dependencies, simple and portable.

---

##  What it does

- Recursively scans a directory and reads alignment files:
  - **NEXUS** (`.nex`, `.nexus`): extracts taxa from `TAXLABELS` or, if missing, from the `MATRIX` block (supports interleaved).
  - **FASTA** (`.fasta`, `.fa`, `.fas`): extracts the first token after each `>` header.
- Deduplicates all labels found.
- Writes a **map file** in the format expected by ASTRAL:
  ```
  taxon<TAB>group
  ```
- Optionally consumes a CSV/TSV (`taxon,group`) to apply your **classification** (family, clade, etc.).
- Can also save a unique list of detected taxa for inspection.

---

##  Quick download

- Script: [`make_astral_map.py`](make_astral_map.py)
- README: this file (`README_make_astral_map_EN.md`)

> If you are seeing this in ChatGPT, use the download links provided in the conversation to save the files locally.

---

## Basic usage

```bash
python make_astral_map.py --input ./alignments --out-map astral.map
```

This creates `astral.map` with **each taxon mapped to itself** (default `--default-group species`).

### With your own classification (CSV/TSV)

```bash
python make_astral_map.py \
  --input ./alignments \
  --groups my_groups.csv \
  --out-map astral.map \
  --out-taxa taxa_list.txt \
  --verbose
```

- `my_groups.csv` (with header):
  ```csv
  taxon,group
  Homo_sapiens,Primates
  Pan_troglodytes,Primates
  Mus_musculus,Rodents
  ```
- Or **without** header, just two columns: `taxon,group`.

---

## 🔧 Parameters

- `--input DIR` **(required)**: directory with your alignments.
- `--out-map FILE` **(required)**: path to the output map file (TSV).
- `--groups FILE` (optional): CSV/TSV with columns `taxon,group`.
- `--out-taxa FILE` (optional): saves a unique list of detected taxa.
- `--pattern` (optional): comma-separated glob patterns. Default:
  ```
  "*.nex,*.nexus,*.fasta,*.fa,*.fas"
  ```
- `--default-group {species,NA,none}` (optional, default: `species`):
  - `species`: maps each taxon to itself (useful when ASTRAL expects “species = taxon name”).
  - `NA`: writes `NA` when no mapping is found in the CSV.
  - `none`: writes an empty string for the group.
- `--strict` (optional): fails if any file has an unknown/invalid format.
- `--verbose` (optional): progress messages on `stderr`.

---

##  How the parser works

- **NEXUS**
  - Removes comments `[ ... ]`.
  - Searches for `TAXLABELS` and tokenizes respecting `'quotes'` and `"quotes"`.
  - If `TAXLABELS` is missing, extracts the **first token of each line** in the `MATRIX` block until `;`.
- **FASTA**
  - For each line starting with `>`, uses the **first token** of the header as the taxon name.

> Tip: If your FASTA headers have spaces (e.g., `>Homo_sapiens mitogenome`), only `Homo_sapiens` will be used. Standardize your labels according to ASTRAL needs.

---

##  Example commands

- “Quick map, no external classification”
  ```bash
  python make_astral_map.py --input ./nexus_fasta --out-map astral.map
  ```

- “Using my taxonomy and inspecting detected taxa”
  ```bash
  python make_astral_map.py \
    --input ./nexus_fasta \
    --groups taxonomy.csv \
    --out-map astral.map \
    --out-taxa detected_taxa.txt \
    --verbose
  ```

- “Only match certain file patterns”
  ```bash
  python make_astral_map.py \
    --input ./alignments \
    --pattern "*.nex,*.fas" \
    --out-map astral.map
  ```

---

##  Best practices

- **Consistent labels**: avoid spaces and exotic characters in names (the script supports quotes in NEXUS, but it’s safer to standardize).
- **Group CSV**: keep a single “source of truth” file for your taxonomy. This avoids divergence across projects.
- **Quick verification**: generate `--out-taxa` and check if all taxa were captured as expected.

---

## 🩹 Common issues

- **“No taxa found”**: check that the directory and patterns (`--pattern`) are correct.
- **Long FASTA labels**: only the first token after `>` is used by default. Adjust your headers or edit the script if you want another logic.
- **Unsupported formats**: PHYLIP and others are not parsed *in this tool*. Convert to FASTA/NEXUS or extend the script.

---

## Map file output

- TSV format: `taxon<TAB>group`, sorted alphabetically by taxon.
- Examples:
  - Without group CSV (`--default-group species`):
    ```text
    Homo_sapiens\tHomo_sapiens
    Pan_troglodytes\tPan_troglodytes
    ```
  - With group CSV:
    ```text
    Homo_sapiens\tPrimates
    Pan_troglodytes\tPrimates
    ```

---

## Tiago Belintani 2025 #brave the sun
