#!/usr/bin/env bash
# Download reference genome resources and SILVA classifier.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
REFERENCE_DIR="${REPO_ROOT}/data/external/reference"

# Configuration
REFERENCE_OUTPUT="${REFERENCE_OUTPUT:-${REFERENCE_DIR}/GRCm39}"
SILVA_CLASSIFIER_PATH="${SILVA_CLASSIFIER_PATH:-${REFERENCE_DIR}/silva-138-99-nb-classifier.qza}"
FASTA_URL="${FASTA_URL:-ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz}"
TRANSCRIPTOME_URL="${TRANSCRIPTOME_URL:-ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_rna.fna.gz}"
GTF_URL="${GTF_URL:-https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gtf.gz}"
SILVA_CLASSIFIER_URL="${SILVA_CLASSIFIER_URL:-https://data.qiime2.org/2023.5/common/silva-138-99-nb-classifier.qza}"

mkdir -p "${REFERENCE_OUTPUT}" "$(dirname "${SILVA_CLASSIFIER_PATH}")"

download_and_extract() {
  local url="$1"
  local dest_dir="$2"
  local outfile="${dest_dir}/$(basename "${url}" .gz)"

  if [[ -f "${outfile}" ]]; then
    echo "[reference] ${outfile} already exists, skipping."
    return
  fi

  echo "[reference] Downloading $(basename "${outfile}")"
  wget -q -O - "${url}" | gunzip -c > "${outfile}"
}

echo "[reference] Downloading genome resources into ${REFERENCE_OUTPUT}"
download_and_extract "${FASTA_URL}" "${REFERENCE_OUTPUT}"
download_and_extract "${TRANSCRIPTOME_URL}" "${REFERENCE_OUTPUT}"
download_and_extract "${GTF_URL}" "${REFERENCE_OUTPUT}"

if [[ -f "${SILVA_CLASSIFIER_PATH}" ]]; then
  echo "[reference] SILVA classifier already present at ${SILVA_CLASSIFIER_PATH}, skipping."
else
  echo "[reference] Downloading SILVA classifier to ${SILVA_CLASSIFIER_PATH}"
  wget -q -O "${SILVA_CLASSIFIER_PATH}" "${SILVA_CLASSIFIER_URL}"
fi

echo "[reference] Reference downloads complete."
