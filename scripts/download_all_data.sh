#!/usr/bin/env bash
# Convenience wrapper to download all raw/reference data using bundled scripts.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Centralized configuration -------------------------------------------------
: "${ENA_PROJECT:=PRJEB102830}"
: "${RNASEQ_METADATA:=${REPO_ROOT}/data/metadata/rnaseq.tsv}"
: "${MICROBIOME_METADATA:=${REPO_ROOT}/data/metadata/microbiome.tsv}"
: "${RNASEQ_OUTPUT:=${REPO_ROOT}/data/raw/rnaseq}"
: "${MICROBIOME_OUTPUT:=${REPO_ROOT}/data/raw/microbiome}"
: "${ENA_REPORT_PATH:=${REPO_ROOT}/data/metadata/ena_file_report.tsv}"
: "${REFERENCE_OUTPUT:=${REPO_ROOT}/data/external/reference/GRCm39}"
: "${SILVA_CLASSIFIER_PATH:=${REPO_ROOT}/data/external/reference/silva-138-99-nb-classifier.qza}"

export ENA_PROJECT \
       RNASEQ_METADATA \
       MICROBIOME_METADATA \
       RNASEQ_OUTPUT \
       MICROBIOME_OUTPUT \
       ENA_REPORT_PATH \
       REFERENCE_OUTPUT \
       SILVA_CLASSIFIER_PATH

echo "[download] Creating directory structure"
bash "${SCRIPT_DIR}/setup_directories.sh"

echo "[download] Reference genomes and SILVA classifier"
bash "${REPO_ROOT}/code/00_data_retrieval/download_reference_data.sh"

echo "[download] RNA-seq FASTQs from ENA"
bash "${REPO_ROOT}/code/00_data_retrieval/download_ena_rnaseq.sh"

echo "[download] Microbiome FASTQs from ENA"
bash "${REPO_ROOT}/code/00_data_retrieval/download_ena_microbiome.sh"

echo "[download] Run verification"
Rscript "${REPO_ROOT}/code/00_data_retrieval/verify_downloads.R"

echo "[download] Done"
