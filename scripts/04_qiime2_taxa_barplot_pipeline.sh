#!/usr/bin/env bash
set -euo pipefail

# 04_qiime2_taxa_barplot_pipeline.sh
#
# Inputs:
# - manifest file (from scripts/02_make_manifest.R)
# - metadata file (QIIME2 compatible)
# - classifier.qza (pre-trained)
#
# Outputs:
# - QIIME2 artifacts and visualizations
# - level-5.csv, level-6.csv, level-7.csv (collapsed feature tables exported as CSV)
#
# IMPORTANT: This script assumes a QIIME2 conda env exists and is activatable.
# Default env name: qiime2-amplicon-2026.1 (from scripts/03_install_qiime2_2026.1.sh)

usage() {
  cat <<'EOF'
Usage:
  bash scripts/04_qiime2_taxa_barplot_pipeline.sh \
    --manifest manifests/manifest_paired.txt \
    --metadata metadata.txt \
    --classifier classifier.qza \
    --out_dir qiime2_out

Required:
  --manifest    Paired-end manifest (tab-delimited; sample-id, forward-absolute-filepath, reverse-absolute-filepath)
  --metadata    Sample metadata file (QIIME2 TSV)
  --classifier  Pre-trained classifier artifact (classifier.qza)

Optional:
  --out_dir     Output directory (default: qiime2_out)
  --qiime_env   Conda env name (default: qiime2-amplicon-2026.1)

  # DADA2 params (defaults match your requested pipeline)
  --trunc_f     (default: 248)
  --trunc_r     (default: 233)
  --trim_f      (default: 18)
  --trim_r      (default: 22)

  --threads     DADA2 threads (default: 0 = all cores)
EOF
}

MANIFEST=""
METADATA=""
CLASSIFIER=""
OUT_DIR="qiime2_out"
QIIME_ENV="qiime2-amplicon-2026.1"

TRUNC_F=248
TRUNC_R=233
TRIM_F=18
TRIM_R=22
THREADS=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2 ;;
    --metadata) METADATA="$2"; shift 2 ;;
    --classifier) CLASSIFIER="$2"; shift 2 ;;
    --out_dir) OUT_DIR="$2"; shift 2 ;;
    --qiime_env) QIIME_ENV="$2"; shift 2 ;;
    --trunc_f) TRUNC_F="$2"; shift 2 ;;
    --trunc_r) TRUNC_R="$2"; shift 2 ;;
    --trim_f) TRIM_F="$2"; shift 2 ;;
    --trim_r) TRIM_R="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "${MANIFEST}" || -z "${METADATA}" || -z "${CLASSIFIER}" ]]; then
  echo "ERROR: --manifest, --metadata, and --classifier are required."
  usage
  exit 1
fi
if [[ ! -f "${MANIFEST}" ]]; then echo "ERROR: manifest not found: ${MANIFEST}"; exit 1; fi
if [[ ! -f "${METADATA}" ]]; then echo "ERROR: metadata not found: ${METADATA}"; exit 1; fi
if [[ ! -f "${CLASSIFIER}" ]]; then echo "ERROR: classifier not found: ${CLASSIFIER}"; exit 1; fi

# Activate conda env
if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found. Install conda and QIIME2 first."
  exit 1
fi
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${QIIME_ENV}"

mkdir -p "${OUT_DIR}"
RUN_ROOT="$(pwd)"

# Work inside out_dir to keep outputs contained
pushd "${OUT_DIR}" >/dev/null

cp -f "${MANIFEST}" ./manifest.txt
cp -f "${METADATA}" ./metadata.txt
cp -f "${CLASSIFIER}" ./classifier.qza

# Detect paired vs single by header columns in manifest
header="$(head -n 1 manifest.txt | tr -d '\r')"
ncol="$(awk -F'\t' 'NR==1{print NF}' manifest.txt)"

if [[ "${ncol}" -eq 3 ]]; then
  echo "[INFO] Detected paired-end manifest (3 columns). Importing as PairedEndSequencesWithQuality..."
  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path manifest.txt \
    --output-path PE_seqs.qza \
    --input-format PairedEndFastqManifestPhred33V2

  echo "[INFO] Running DADA2 denoise-paired..."
  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs PE_seqs.qza \
    --p-trunc-len-f "${TRUNC_F}" \
    --p-trunc-len-r "${TRUNC_R}" \
    --p-trim-left-f "${TRIM_F}" \
    --p-trim-left-r "${TRIM_R}" \
    --p-n-threads "${THREADS}" \
    --o-table feature_table.qza \
    --o-representative-sequences representative_sequences.qza \
    --o-denoising-stats dada2_stats.qza
else
  echo "ERROR: manifest.txt does not look like a paired-end manifest (expected 3 tab-separated columns)."
  echo "Header: ${header}"
  exit 1
fi

echo "[INFO] Summarizing feature table..."
qiime feature-table summarize \
  --i-table feature_table.qza \
  --m-sample-metadata-file metadata.txt \
  --o-visualization feature_table_summary.qzv

echo "[INFO] Tabulating representative sequences..."
qiime feature-table tabulate-seqs \
  --i-data representative_sequences.qza \
  --o-visualization representative_sequences_summary.qzv

echo "[INFO] Taxonomic classification (classify-sklearn)..."
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads representative_sequences.qza \
  --o-classification taxonomy.qza

echo "[INFO] Taxonomy summary..."
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy_summary.qzv

echo "[INFO] Taxa barplot..."
qiime taxa barplot \
  --i-table feature_table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metadata.txt \
  --o-visualization taxa_bar_plots.qzv

# ---- Export collapsed tables at levels 5, 6, 7 into CSV ----
export_level() {
  local LVL="$1"
  local OUTCSV="$2"
  local TMPDIR="export_level_${LVL}"

  echo "[INFO] Collapsing to taxonomy level ${LVL} ..."
  qiime taxa collapse \
    --i-table feature_table.qza \
    --i-taxonomy taxonomy.qza \
    --p-level "${LVL}" \
    --o-collapsed-table "table_L${LVL}.qza"

  rm -rf "${TMPDIR}"
  mkdir -p "${TMPDIR}"
  qiime tools export \
    --input-path "table_L${LVL}.qza" \
    --output-path "${TMPDIR}"

  # Convert biom → TSV → CSV
  if ! command -v biom >/dev/null 2>&1; then
    echo "ERROR: biom command not found inside QIIME2 environment."
    echo "Install biom-format or use an environment that includes it."
    exit 1
  fi

  biom convert \
    -i "${TMPDIR}/feature-table.biom" \
    -o "${TMPDIR}/feature-table.tsv" \
    --to-tsv

  # Remove the first comment line from biom TSV
  tail -n +2 "${TMPDIR}/feature-table.tsv" > "${TMPDIR}/feature-table.nohdr.tsv"

  # Convert tab → comma (keeps first column as taxon IDs)
  python - <<'PY'
import csv, sys
from pathlib import Path

lvl = sys.argv[1]
in_tsv = Path(sys.argv[2])
out_csv = Path(sys.argv[3])

with in_tsv.open("r", newline="") as f_in, out_csv.open("w", newline="") as f_out:
    reader = csv.reader(f_in, delimiter="\t")
    writer = csv.writer(f_out)
    for i, row in enumerate(reader):
        if i == 0 and row and row[0].strip() == "#OTU ID":
            row[0] = "Taxon"
        writer.writerow(row)
PY "${LVL}" "${TMPDIR}/feature-table.nohdr.tsv" "${OUTCSV}"

  echo "[DONE] Wrote ${OUTCSV}"
}

export_level 5 "level-5.csv"
export_level 6 "level-6.csv"
export_level 7 "level-7.csv"

# Copy level CSVs to run root for downstream R scripts
cp -f level-5.csv "${RUN_ROOT}/level-5.csv"
cp -f level-6.csv "${RUN_ROOT}/level-6.csv"
cp -f level-7.csv "${RUN_ROOT}/level-7.csv"

echo "[DONE] QIIME2 pipeline complete. Outputs in: ${OUT_DIR}"
echo "[DONE] level-5/6/7.csv copied to repo root."

popd >/dev/null
