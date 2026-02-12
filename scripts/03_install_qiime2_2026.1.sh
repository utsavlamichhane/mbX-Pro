#!/usr/bin/env bash
set -euo pipefail

# 03_install_qiime2_2026.1.sh
# One-shot installer for:
# - Miniconda (if conda is missing)
# - QIIME 2 Amplicon 2026.1 conda environment
#
# This follows your INSTALLATION.md approach (Miniconda → conda update → conda env create from QIIME2 distributions YAML),
# but upgrades the QIIME 2 version to 2026.1 and makes it non-interactive.

usage() {
  cat <<'EOF'
Usage:
  bash scripts/03_install_qiime2_2026.1.sh [--env qiime2-amplicon-2026.1] [--conda_prefix $HOME/miniconda3]

Optional:
  --env          Conda environment name (default: qiime2-amplicon-2026.1)
  --conda_prefix Where to install Miniconda if needed (default: $HOME/miniconda3)
EOF
}

ENV_NAME="qiime2-amplicon-2026.1"
CONDA_PREFIX="${HOME}/miniconda3"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --env) ENV_NAME="$2"; shift 2 ;;
    --conda_prefix) CONDA_PREFIX="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
done

OS="$(uname -s)"
ARCH="$(uname -m)"

# ---- Install Miniconda if conda is missing ----
if ! command -v conda >/dev/null 2>&1; then
  echo "[INFO] conda not found. Installing Miniconda into: ${CONDA_PREFIX}"

  mkdir -p "${CONDA_PREFIX}"
  INSTALLER=""
  URL=""

  if [[ "${OS}" == "Darwin" && "${ARCH}" == "x86_64" ]]; then
    INSTALLER="Miniconda3-latest-MacOSX-x86_64.sh"
    URL="https://repo.anaconda.com/miniconda/${INSTALLER}"
  elif [[ "${OS}" == "Linux" && "${ARCH}" == "x86_64" ]]; then
    INSTALLER="Miniconda3-latest-Linux-x86_64.sh"
    URL="https://repo.anaconda.com/miniconda/${INSTALLER}"
  else
    echo "ERROR: Auto-install only supports Intel macOS (x86_64) and Linux x86_64."
    echo "You are on OS=${OS} ARCH=${ARCH}."
    echo "Install conda (Miniconda/Miniforge) manually, then re-run this script."
    exit 1
  fi

  curl -L "${URL}" -o "${CONDA_PREFIX}/miniconda.sh"
  bash "${CONDA_PREFIX}/miniconda.sh" -b -u -p "${CONDA_PREFIX}"
  rm -f "${CONDA_PREFIX}/miniconda.sh"

  export PATH="${CONDA_PREFIX}/bin:${PATH}"
fi

# Enable conda activation in this script
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"

echo "[INFO] Updating conda..."
conda update -y conda

# ---- QIIME2 distributions YAMLs for 2026.1 (Amplicon) ----
YML_LINUX="https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.1/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml"
YML_MACOS="https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.1/amplicon/released/qiime2-amplicon-macos-latest-conda.yml"

echo "[INFO] Installing QIIME 2 2026.1 env: ${ENV_NAME}"
echo "[INFO] Detected OS=${OS} ARCH=${ARCH}"

# If env already exists, skip create
if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "[INFO] Conda env already exists: ${ENV_NAME} (skipping creation)"
else
  if [[ "${OS}" == "Linux" ]]; then
    conda env create --name "${ENV_NAME}" --file "${YML_LINUX}"
  elif [[ "${OS}" == "Darwin" ]]; then
    if [[ "${ARCH}" == "arm64" ]]; then
      echo "ERROR: Apple Silicon detected (arm64)."
      echo "QIIME 2 on Apple Silicon typically installs via osx-64 (Rosetta) or via a compatible conda stack."
      echo "Install conda first, then create env using the macOS YAML with osx-64 subdir."
      exit 1
    fi
    conda env create --name "${ENV_NAME}" --file "${YML_MACOS}"
  else
    echo "ERROR: Unsupported OS: ${OS}"
    exit 1
  fi
fi

echo "[INFO] Testing install..."
conda deactivate || true
conda activate "${ENV_NAME}"
qiime info

echo "[DONE] QIIME 2 installed in conda env: ${ENV_NAME}"
echo "Next: conda activate ${ENV_NAME}"
