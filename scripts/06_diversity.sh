#!/usr/bin/env bash
set -euo pipefail

# 06_diversity.sh
# Runs QIIME2 phylogeny + core-metrics-phylogenetic to compute alpha and beta diversity outputs.
# Writes everything into a separate directory (default: diversity_out).

usage() {
  cat <<'EOF'
Usage:
  bash scripts/06_diversity.sh \
    --table qiime2_out/feature_table.qza \
    --rep_seqs qiime2_out/representative_sequences.qza \
    --metadata metadata.txt \
    --sampling_depth 10000 \
    [--out_dir diversity_out] \
    [--qiime_env qiime2-amplicon-2026.1]

Required:
  --table           feature_table.qza
  --rep_seqs        representative_sequences.qza
  --metadata        metadata file
  --sampling_depth  integer sampling depth for rarefaction (choose based on your data)

Optional:
  --out_dir         output directory (default: diversity_out)
  --qiime_env       conda env name (default: qiime2-amplicon-2026.1)
EOF
}

TABLE=""
REP_SEQS=""
METADATA=""
DEPTH=""
OUT_DIR="diversity_out"
QIIME_ENV="qiime2-amplicon-2026.1"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --table) TABLE="$2"; shift 2 ;;
    --rep_seqs) REP_SEQS="$2"; shift 2 ;;
    --metadata) METADATA="$2"; shift 2 ;;
    --sampling_depth) DEPTH="$2"; shift 2 ;;
    --out_dir) OUT_DIR="$2"; shift 2 ;;
    --qiime_env) QIIME_ENV="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "${TABLE}" || -z "${REP_SEQS}" || -z "${METADATA}" || -z "${DEPTH}" ]]; then
  echo "ERROR: --table, --rep_seqs, --metadata, --sampling_depth are required."
  usage
  exit 1
fi
[[ -f "${TABLE}" ]] || { echo "ERROR: table not found: ${TABLE}"; exit 1; }
[[ -f "${REP_SEQS}" ]] || { echo "ERROR: rep_seqs not found: ${REP_SEQS}"; exit 1; }
[[ -f "${METADATA}" ]] || { echo "ERROR: metadata not found: ${METADATA}"; exit 1; }

# Activate QIIME2
if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found. Install conda and QIIME2 first."
  exit 1
fi
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${QIIME_ENV}"

mkdir -p "${OUT_DIR}"
RUN_ROOT="$(pwd)"
pushd "${OUT_DIR}" >/dev/null

cp -f "${TABLE}" ./feature_table.qza
cp -f "${REP_SEQS}" ./representative_sequences.qza
cp -f "${METADATA}" ./metadata.txt

echo "[INFO] Building phylogenetic tree..."
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences representative_sequences.qza \
  --o-alignment aligned_rep_seqs.qza \
  --o-masked-alignment masked_aligned_rep_seqs.qza \
  --o-tree unrooted_tree.qza \
  --o-rooted-tree rooted_tree.qza

echo "[INFO] Running core-metrics-phylogenetic (sampling depth = ${DEPTH})..."
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted_tree.qza \
  --i-table feature_table.qza \
  --p-sampling-depth "${DEPTH}" \
  --m-metadata-file metadata.txt \
  --output-dir core_metrics

# Export key results to plain-text formats
mkdir -p exports

echo "[INFO] Exporting alpha diversity vectors..."
for a in faith_pd_vector.qza evenness_vector.qza shannon_vector.qza observed_features_vector.qza; do
  if [[ -f "core_metrics/${a}" ]]; then
    outd="exports/${a%.qza}"
    mkdir -p "${outd}"
    qiime tools export --input-path "core_metrics/${a}" --output-path "${outd}"
  fi
done

echo "[INFO] Exporting beta diversity distance matrices..."
for b in unweighted_unifrac_distance_matrix.qza weighted_unifrac_distance_matrix.qza jaccard_distance_matrix.qza bray_curtis_distance_matrix.qza; do
  if [[ -f "core_metrics/${b}" ]]; then
    outd="exports/${b%.qza}"
    mkdir -p "${outd}"
    qiime tools export --input-path "core_metrics/${b}" --output-path "${outd}"
  fi
done

echo "[DONE] Diversity outputs in: ${OUT_DIR}"
popd >/dev/null
