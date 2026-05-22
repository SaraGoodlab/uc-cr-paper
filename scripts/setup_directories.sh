#!/usr/bin/env bash
# Create the fixed directory skeleton used by this project.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

mkdir -p \
  "${REPO_ROOT}/data/raw/rnaseq/colon" \
  "${REPO_ROOT}/data/raw/rnaseq/spleen" \
  "${REPO_ROOT}/data/raw/microbiome" \
  "${REPO_ROOT}/data/processed/rnaseq/nfcore/star_rsem" \
  "${REPO_ROOT}/data/processed/rnaseq/rsem_counts" \
  "${REPO_ROOT}/data/processed/microbiome/tree" \
  "${REPO_ROOT}/data/external/reference/GRCm39" \
  "${REPO_ROOT}/data/metadata" \
  "${REPO_ROOT}/results/rnaseq/deseq2_results" \
  "${REPO_ROOT}/results/rnaseq/deg_lists" \
  "${REPO_ROOT}/results/rnaseq/quality_metrics" \
  "${REPO_ROOT}/results/pathway_analysis/mouse_gsea" \
  "${REPO_ROOT}/results/microbiome/diversity" \
  "${REPO_ROOT}/results/microbiome/prevalence" \
  "${REPO_ROOT}/results/microbiome/differential_abundance" \
  "${REPO_ROOT}/figures/rnaseq/pca_plots" \
  "${REPO_ROOT}/figures/pathway_analysis/gsea_plots" \
  "${REPO_ROOT}/figures/microbiome/diversity_plots" \
  "${REPO_ROOT}/figures/microbiome/prevalence" \
  "${REPO_ROOT}/figures/microbiome/abundance_plots" \
  "${REPO_ROOT}/tables/rnaseq" \
  "${REPO_ROOT}/tables/microbiome/prevalence" \
  "${REPO_ROOT}/tables/microbiome/abundance" \
  "${REPO_ROOT}/logs/01_rnaseq" \
  "${REPO_ROOT}/logs/04_microbiome"

echo "Created directory skeleton under ${REPO_ROOT}"
