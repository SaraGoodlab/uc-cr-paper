#!/usr/bin/env bash
# Run nf-core/rnaseq for the the UC-CR mouse RNA-seq dataset.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

PIPELINE="nf-core/rnaseq"
REVISION="3.14.0"
PROFILE="singularity"
SAMPLESHEET="${REPO_ROOT}/data/metadata/nfcore_samplesheet.csv"
PARAMS_FILE="${REPO_ROOT}/config/nfcore_rnaseq_params.json"
OUTDIR="${REPO_ROOT}/data/processed/rnaseq/nfcore"
WORKDIR="${OUTDIR}/work"

command -v nextflow >/dev/null 2>&1 || { echo "nextflow not found in PATH"; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "python3 not found in PATH"; exit 1; }

mkdir -p "${OUTDIR}" "${WORKDIR}" "$(dirname "${SAMPLESHEET}")"
python3 "${REPO_ROOT}/scripts/build_nfcore_samplesheet.py"

cd "${REPO_ROOT}"

nextflow run "${PIPELINE}" \
  -r "${REVISION}" \
  -name mice_rnaseq \
  -profile "${PROFILE}" \
  -work-dir "${WORKDIR}" \
  -params-file "${PARAMS_FILE}"

echo "[nf-core] Finished. STAR/RSEM outputs expected in ${OUTDIR}/star_rsem/"
