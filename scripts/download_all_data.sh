#!/usr/bin/env bash
# Download the fixed raw/reference inputs for this project.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

echo "[download] Creating directory structure"
bash "${SCRIPT_DIR}/setup_directories.sh"

echo "[download] Reference genomes and SILVA classifier"
bash "${REPO_ROOT}/code/00_data_retrieval/download_reference_data.sh"

echo "[download] RNA-seq FASTQs from ENA"
bash "${REPO_ROOT}/code/00_data_retrieval/download_ena_rnaseq.sh"

echo "[download] Microbiome FASTQs from ENA"
bash "${REPO_ROOT}/code/00_data_retrieval/download_ena_microbiome.sh"

echo "[download] Done"
