#!/usr/bin/env bash
# Helper script to run nf-core/rnaseq outside of Snakemake

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
BASE_DIR="${BASE_DIR:-${REPO_ROOT}}"

RNASEQ_PIPELINE="${RNASEQ_PIPELINE:-nf-core/rnaseq}"
METADATA="${METADATA:-${BASE_DIR}/data/metadata/rnaseq.tsv}"
FASTQ_DIR="${FASTQ_DIR:-${BASE_DIR}/data/raw/rnaseq}"
SAMPLESHEET="${SAMPLESHEET:-${BASE_DIR}/data/metadata/nfcore_samplesheet.csv}"
OUTDIR="${OUTDIR:-${BASE_DIR}/data/processed/rnaseq/nfcore}"
WORKDIR="${WORKDIR:-${OUTDIR}/work}"
PARAMS_FILE="${PARAMS_FILE:-${BASE_DIR}/config/nfcore_rnaseq_params.json}"
PROFILE="${PROFILE:-singularity}"
REVISION="${REVISION:-3.14.0}"
ALIGNER="${ALIGNER:-star}"
QUANTIFIER="${QUANTIFIER:-rsem}"
FASTA_PATH="${FASTA_PATH:-${BASE_DIR}/data/external/reference/GRCm39/GCF_000001635.27_GRCm39_genomic.fna}"
GTF_PATH="${GTF_PATH:-${BASE_DIR}/data/external/reference/GRCm39/GCF_000001635.27_GRCm39_genomic.gtf}"
TRANSCRIPTOME_PATH="${TRANSCRIPTOME_PATH:-${BASE_DIR}/data/external/reference/GRCm39/GCF_000001635.27_GRCm39_rna.fna}"
STRANDEDNESS="${STRANDEDNESS:-auto}"
NEXTFLOW_NAME="${NEXTFLOW_NAME:-mice_rnaseq}"
EXTRA_ARGS="${EXTRA_ARGS:-}"
NFCORE_FILTER_ARGS=()

NFCORE_ALIGNER="${ALIGNER}"
if [[ "${ALIGNER}" == "star" ]]; then
  case "${QUANTIFIER}" in
    rsem|salmon) NFCORE_ALIGNER="star_${QUANTIFIER}" ;;
  esac
fi

if [[ ! -f "${METADATA}" ]]; then
  echo "Metadata file not found: ${METADATA}"
  exit 1
fi

if [[ ! -d "${FASTQ_DIR}" ]]; then
  echo "FASTQ directory not found: ${FASTQ_DIR}"
  exit 1
fi

mkdir -p "${OUTDIR}" "${WORKDIR}" "$(dirname "${SAMPLESHEET}")"

echo "[nf-core] Repository root: ${BASE_DIR}"
echo "[nf-core] FASTQ dir: ${FASTQ_DIR}"
echo "[nf-core] Metadata: ${METADATA}"
echo "[nf-core] Samplesheet: ${SAMPLESHEET}"
echo "[nf-core] Output dir: ${OUTDIR}"
echo "[nf-core] Work dir: ${WORKDIR}"
echo "[nf-core] Pipeline: ${RNASEQ_PIPELINE}"
echo "[nf-core] Params template: ${PARAMS_FILE}"
echo "[nf-core] Profile: ${PROFILE}"
echo "[nf-core] Revision: ${REVISION}"
echo "[nf-core] Aligner: ${NFCORE_ALIGNER}"

REVISION_ARGS=()
if [[ ! -d "${RNASEQ_PIPELINE}" ]]; then
  REVISION_ARGS=(-r "${REVISION}")
fi

command -v nextflow >/dev/null 2>&1 || { echo "nextflow not found in PATH"; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "python3 not found in PATH"; exit 1; }

python3 "${BASE_DIR}/scripts/build_nfcore_samplesheet.py" \
  --metadata "${METADATA}" \
  --fastq_dir "${FASTQ_DIR}" \
  --output "${SAMPLESHEET}" \
  --strandedness "${STRANDEDNESS}"

EFFECTIVE_GTF="${GTF_PATH}"
FILTER_SCRIPT="${RNASEQ_PIPELINE}/workflow/bin/filter_gtf.py"

if [[ -d "${RNASEQ_PIPELINE}" ]]; then
  if [[ -f "${FILTER_SCRIPT}" ]]; then
    FILTER_PREFIX="${WORKDIR}/gtf_filtered"
    python3 "${FILTER_SCRIPT}" \
      --gtf "${GTF_PATH}" \
      --fasta "${FASTA_PATH}" \
      --prefix "${FILTER_PREFIX}"
    EFFECTIVE_GTF="${FILTER_PREFIX}.filtered.gtf"
    NFCORE_FILTER_ARGS+=(--skip_gtf_filter)
    echo "[nf-core] Using filtered GTF: ${EFFECTIVE_GTF}"
  else
    echo "[nf-core] Warning: ${FILTER_SCRIPT} not found; skipping pre-run GTF filtering."
  fi
fi

NFCORE_PARAMS_FILE="$(mktemp -p "${WORKDIR}" nfcore_params.XXXXXX.json)"
cleanup_params() {
  rm -f "${NFCORE_PARAMS_FILE}"
}
trap cleanup_params EXIT

echo "[nf-core] Writing parameter overrides to ${NFCORE_PARAMS_FILE}"
python3 "${BASE_DIR}/scripts/make_nfcore_params_json.py" \
  --template "${PARAMS_FILE}" \
  --output "${NFCORE_PARAMS_FILE}" \
  --samplesheet "${SAMPLESHEET}" \
  --outdir "${OUTDIR}" \
  --fasta "${FASTA_PATH}" \
  --gtf "${EFFECTIVE_GTF}" \
  --transcriptome "${TRANSCRIPTOME_PATH}" \
  --aligner "${NFCORE_ALIGNER}"

echo "[nf-core] Launching nf-core/rnaseq"
nextflow run "${RNASEQ_PIPELINE}" \
  "${REVISION_ARGS[@]}" \
  -name "${NEXTFLOW_NAME}" \
  -profile "${PROFILE}" \
  -work-dir "${WORKDIR}" \
  -params-file "${NFCORE_PARAMS_FILE}" \
  "${NFCORE_FILTER_ARGS[@]}" \
  "${EXTRA_ARGS}"

echo "[nf-core] Finished. STAR/RSEM outputs expected in ${OUTDIR}/star_rsem/"
