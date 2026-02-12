# 16S rRNA amplicon end-to-end pipeline (NCBI FASTQ → QIIME 2 → mbX → diversity → ANCOM-BC & MaAsLin2)

This repo is an **end-to-end, reproducible** workflow for paired-end 16S rRNA gene amplicon data:

1) Download FASTQs from NCBI SRA (via SRA Toolkit)  
2) Build a QIIME 2 manifest from a FASTQ directory  
3) Install QIIME 2 Amplicon **2026.1** (Conda)  
4) Run a QIIME 2 pipeline to produce taxa barplots + `level-5.csv`, `level-6.csv`, `level-7.csv`  
5) Run **mbX** in R for cleaning, visualization, and stats (always using `level-7.csv`)  
6) Run alpha/beta diversity into a dedicated directory  
7) Run differential abundance with **ANCOM-BC** and **MaAsLin2** into separate directories  

---

## Assumptions (explicit)

These are *defaults* in the scripts (edit via CLI flags or by editing the scripts):

- **Paired-end** reads (R1/R2) with **Phred33** qualities.
- DADA2 parameters (can be changed in `scripts/04_qiime2_taxa_barplot_pipeline.sh`):
  - `--p-trunc-len-f 248`
  - `--p-trunc-len-r 233`
  - `--p-trim-left-f 18`
  - `--p-trim-left-r 22`
- Taxonomic table exports:
  - `level-5.csv` = QIIME 2 **taxonomic level 5** (Family in Greengenes/SILVA)
  - `level-6.csv` = level 6 (Genus)
  - `level-7.csv` = level 7 (Species)
- A **pre-trained classifier artifact** exists (e.g. `classifier.qza`) and matches your amplicon region/primers.

---

## Repo layout

```
.
├── config/
│   └── srr_ids.txt                 ###IMP: you need to manually give the ids in the file, one id (like SRR31972796) in each line
├── scripts/
│   ├── 01_download_fastq_ncbi.sh
│   ├── 02_make_manifest.R
│   ├── 03_install_qiime2_2026.1.sh
│   ├── 04_qiime2_taxa_barplot_pipeline.sh
│   ├── 05_run_mbx.R
│   ├── 06_diversity.sh
│   ├── 07_ancombc.R
│   └── 08_maaslin2.R
└── README.md
```

---

## Quickstart (Please use this order, similar to what we did in our video)

### 0) Put your SRR accessions in a file

Create `config/srr_ids.txt` with one SRR per line, e.g.

```
SRR31972796
SRR31972797
SRR31972798
SRR31972799
```


---

### 1) Download FASTQs from NCBI

This downloads with `fasterq-dump`, compresses with `pigz` if available, and renames to **Casava 1.8** format:

```bash
bash scripts/01_download_fastq_ncbi.sh \
  --srr_file config/srr_ids.txt \
  --out_dir data/fastq \
  --threads 8
```

Output: `data/fastq/*.fastq.gz`

---

### 2) Create a manifest from the FASTQ directory

```bash
Rscript scripts/02_make_manifest.R \
  --fastq_dir data/fastq \
  --read_type paired \
  --out manifests/manifest_paired.txt
```

Output: `manifests/manifest_paired.txt` (tab-delimited)

---

### 3) Install QIIME 2 Amplicon 2026.1

This script creates a conda env named `qiime2-amplicon-2026.1` using the **official QIIME 2 Library** install YAMLs.

```bash
bash scripts/03_install_qiime2_2026.1.sh
```

Then:

```bash
conda activate qiime2-amplicon-2026.1
qiime info
```

---

### 4) Run QIIME 2 pipeline (taxa bar plots + level tables)

You must provide:
- `--manifest` from step 2
- `--metadata` sample metadata file
- `--classifier` a pre-trained `classifier.qza`

```bash
bash scripts/04_qiime2_taxa_barplot_pipeline.sh \
  --manifest manifests/manifest_paired.txt \
  --metadata metadata.txt \
  --classifier classifier.qza \
  --out_dir qiime2_out
```

Outputs (in `qiime2_out/`):
- `feature_table.qza`, `representative_sequences.qza`, `dada2_stats.qza`
- `feature_table_summary.qzv`, `representative_sequences_summary.qzv`
- `taxonomy.qza`, `taxonomy_summary.qzv`
- `taxa_bar_plots.qzv`
- `level-5.csv`, `level-6.csv`, `level-7.csv`

And the script also copies `level-5.csv/6/7` into the **repo root** for downstream R scripts.

---

### 5) Run mbX (clean / viz / stats) on level-7.csv

This always uses `level-7.csv`. 

```bash
Rscript scripts/05_run_mbx.R \
  --microbiome level-7.csv \
  --level g
  --metadata metadata.txt \
  --group_var BMIClass \
  --top_taxa 10
```

Output: `mbx_out/` with cleaned tables, plots, and stats (plus `sessionInfo.txt`).

---

### 6) Alpha & Beta diversity (to a separate directory)

This step requires a **sampling depth**

```bash
bash scripts/06_diversity.sh \
  --table qiime2_out/feature_table.qza \
  --rep_seqs qiime2_out/representative_sequences.qza \
  --metadata metadata.txt \
  --sampling_depth 10000 \
  --out_dir diversity_out
```

Outputs: `diversity_out/` (core metrics, alpha, beta distance matrices, PCoA)

---

### 7) Differential abundance: ANCOM-BC and MaAsLin2

Both scripts take:
- `level-7.csv` (counts / frequencies at a taxonomic level)
- metadata file
- a metadata variable to test (you choose)

#### ANCOM-BC

```bash
Rscript scripts/07_ancombc.R \
  --microbiome level-7.csv \
  --metadata metadata.txt \
  --group_var BMIClass \
  --out_dir ancombc_out
```

#### MaAsLin2

```bash
Rscript scripts/08_maaslin2.R \
  --microbiome level-7.csv \
  --metadata metadata.txt \
  --group_var BMIClass \
  --out_dir maaslin2_out
```

Each writes results and `sessionInfo.txt` into its output directory.

---

## Notes on reproducibility

- QIIME 2 is pinned to **Amplicon 2026.1** via the official `qiime2/distributions` YAML references.
- Each R script writes out `sessionInfo()` to record package versions.
- Scripts use explicit, deterministic filenames so downstream steps can rely on upstream outputs.

