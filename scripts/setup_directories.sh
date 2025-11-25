#!/usr/bin/env bash
# Create minimal directory structure for data, results, figures, tables.

set -euo pipefail

BASE_DIR="$(pwd)"

mkdir -p \
  "${BASE_DIR}/data/raw/rnaseq/colon" \
  "${BASE_DIR}/data/raw/rnaseq/spleen" \
  "${BASE_DIR}/data/raw/microbiome" \
  "${BASE_DIR}/data/processed/rnaseq/rsem_counts" \
  "${BASE_DIR}/data/processed/rnaseq/nfcore/star_rsem" \
  "${BASE_DIR}/data/processed/microbiome" \
  "${BASE_DIR}/data/external/reference/GRCm39" \
  "${BASE_DIR}/data/metadata" \
  "${BASE_DIR}/results/rnaseq/deseq2_results" \
  "${BASE_DIR}/results/rnaseq/deg_lists" \
  "${BASE_DIR}/results/rnaseq/quality_metrics" \
  "${BASE_DIR}/figures/rnaseq/pca_plots" \
  "${BASE_DIR}/tables"

echo "Created directory skeleton under ${BASE_DIR}"
