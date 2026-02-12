#!/usr/bin/env bash
set -euo pipefail

##    01_download_fastq_ncbi.sh #this will create data dir as mentioned in the video
# Downloads paired-end FASTQ from NCBI SRA using SRA Toolkit (fasterq-dump),
# compresses with pigz (if available), and renames to Casava 1.8 format.
# Based on "NCBI FASTQ Download Cheatsheet.pdf" (SRA Toolkit v3.2.1, fasterq-dump, pigz, Casava renaming).

usage() {
  cat <<'EOF'
Usage:
  bash scripts/01_download_fastq_ncbi.sh --srr_file config/srr_ids.txt [--out_dir data/fastq] [--threads 8]

Required:
  --srr_file   Path to a text file with one SRR accession per line (e.g., SRR123...)

Optional:
  --out_dir    Output directory for FASTQs (default: data/fastq)
  --threads    Threads for fasterq-dump and pigz (default: 8)
  --sra_ver    SRA toolkit version to auto-install if needed (default: 3.2.1)
  --install_dir Where to install SRA toolkit if fasterq-dump is missing (default: tools)
  --no_rename  Do not rename into Casava format (keeps *_1.fastq.gz and *_2.fastq.gz)
EOF
}

SRR_FILE=""
OUT_DIR="data/fastq"
THREADS="8"
SRA_VER="3.2.1"
INSTALL_DIR="tools"
DO_RENAME="1"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --srr_file) SRR_FILE="$2"; shift 2 ;;
    --out_dir) OUT_DIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --sra_ver) SRA_VER="$2"; shift 2 ;;
    --install_dir) INSTALL_DIR="$2"; shift 2 ;;
    --no_rename) DO_RENAME="0"; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "${SRR_FILE}" ]]; then
  echo "ERROR: --srr_file is required"
  usage
  exit 1
fi
if [[ ! -f "${SRR_FILE}" ]]; then
  echo "ERROR: SRR file not found: ${SRR_FILE}"
  exit 1
fi

mkdir -p "${OUT_DIR}" "${INSTALL_DIR}"

# --- Ensure fasterq-dump exists (auto-install SRA Toolkit if missing) ---
if ! command -v fasterq-dump >/dev/null 2>&1; then
  echo "[INFO] fasterq-dump not found. Installing SRA Toolkit ${SRA_VER} (user-local) ..."
  ARCH="$(uname -m)"
  OS="$(uname -s)"
  FILE=""
  URL=""

  if [[ "${OS}" == "Darwin" ]]; then
    if [[ "${ARCH}" == "arm64" ]]; then
      FILE="sratoolkit.${SRA_VER}-mac-arm64.tar.gz"
    else
      FILE="sratoolkit.${SRA_VER}-mac-x86_64.tar.gz"
    fi
    URL="https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VER}/${FILE}"
  elif [[ "${OS}" == "Linux" ]]; then
    # Cheatsheet specifies Ubuntu64 tarball name:
    FILE="sratoolkit.${SRA_VER}-ubuntu64.tar.gz"
    URL="https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VER}/${FILE}"
  else
    echo "ERROR: Unsupported OS for auto-install: ${OS}"
    echo "Install SRA Toolkit manually and re-run."
    exit 1
  fi

  curl -L -o "${INSTALL_DIR}/${FILE}" "${URL}"
  tar -xvzf "${INSTALL_DIR}/${FILE}" -C "${INSTALL_DIR}"
  # Add to PATH for this script run
  export PATH="$(pwd)/${INSTALL_DIR}/sratoolkit.${SRA_VER}-"*"/bin:${PATH}"
fi

echo "[INFO] SRA Toolkit version:"
fasterq-dump --version || true

# --- Download and process each SRR ---
mapfile -t SRRS < <(grep -v '^[[:space:]]*$' "${SRR_FILE}" | tr -d '\r')

if [[ "${#SRRS[@]}" -eq 0 ]]; then
  echo "ERROR: No SRR accessions found in ${SRR_FILE}"
  exit 1
fi

echo "[INFO] Downloading ${#SRRS[@]} SRR accessions into: ${OUT_DIR}"
echo "[INFO] Threads: ${THREADS}"

pushd "${OUT_DIR}" >/dev/null

for id in "${SRRS[@]}"; do
  echo "[INFO] Downloading: ${id}"
  fasterq-dump "${id}" --split-files -e "${THREADS}" -O ./

  # Compress (prefer pigz; fallback to gzip)
  if command -v pigz >/dev/null 2>&1; then
    [[ -f "${id}_1.fastq" ]] && pigz -p "${THREADS}" "${id}_1.fastq"
    [[ -f "${id}_2.fastq" ]] && pigz -p "${THREADS}" "${id}_2.fastq"
  else
    echo "[WARN] pigz not found; using gzip (slower)."
    [[ -f "${id}_1.fastq" ]] && gzip -f "${id}_1.fastq"
    [[ -f "${id}_2.fastq" ]] && gzip -f "${id}_2.fastq"
  fi

  if [[ "${DO_RENAME}" == "1" ]]; then
    # Rename to Casava 1.8 format (sample-id derived from prefix before _S1_)
    # forward
    if [[ -f "${id}_1.fastq.gz" ]]; then mv "${id}_1.fastq.gz" "${id}_S1_L001_R1_001.fastq.gz"; fi
    if [[ -f "${id}_1.fastq" ]]; then mv "${id}_1.fastq" "${id}_S1_L001_R1_001.fastq"; fi
    # reverse
    if [[ -f "${id}_2.fastq.gz" ]]; then mv "${id}_2.fastq.gz" "${id}_S1_L001_R2_001.fastq.gz"; fi
    if [[ -f "${id}_2.fastq" ]]; then mv "${id}_2.fastq" "${id}_S1_L001_R2_001.fastq"; fi
  fi
done

popd >/dev/null

echo "[DONE] FASTQs are in: ${OUT_DIR}"
