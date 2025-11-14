#!/usr/bin/env bash
# Download 16S rRNA microbiome FASTQ files from ENA

set -euo pipefail

# Configuration
ENA_PROJECT="${ENA_PROJECT:-PRJEB102830}"
OUTPUT_DIR="${OUTPUT_DIR:-../../data/raw/microbiome}"
SAMPLE_MAPPING="${SAMPLE_MAPPING:-../../data/metadata/microbiome.tsv}"
REPORT_PATH="../../data/metadata/ena_file_report.tsv"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "Downloading 16S rRNA microbiome data from ENA project: ${ENA_PROJECT}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Metadata file: ${SAMPLE_MAPPING}"

REPORT_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ENA_PROJECT}&result=read_run&fields=run_accession,sample_accession,submitted_ftp&format=tsv&download=true"

if [[ ! -f "${REPORT_PATH}" ]]; then
    echo "Fetching file report from ENA..."
    wget -q -O "${REPORT_PATH}" "${REPORT_URL}"
fi

declare -A SAMPLE_MAP

echo "Loading sample mapping..."
while IFS=$'\t' read -r sample_accession sample_id _; do
    # Skip header or empty lines
    [[ -z "$sample_accession" ]] && continue
    [[ "$sample_accession" == "sample_accession" ]] && continue

    SAMPLE_MAP["$sample_accession"]="$sample_id"
done < "${SAMPLE_MAPPING}"

# Filter ENA report to rows with a single file (no ';' in submitted_ftp) and download each file

DOWNLOADED=0
SKIPPED=0
ERRORS=0

while IFS=$'\t' read -r run_accession sample_accession submitted_ftp; do

    if [[ -z "${submitted_ftp}" ]]; then
        echo "  [SKIP] no submitted_ftp for run ${run_accession}"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # Look up internal sample_id (handle missing key safely)
    sample_id="${SAMPLE_MAP[$sample_accession]:-}"

    if [[ -z "${sample_id}" ]]; then
        echo "  [SKIP] no mapping for ENA sample ${sample_accession}"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    remote_name="$(basename "${submitted_ftp}")"
    out_file="${sample_id}.fastq.gz"

    if [[ -f "${OUTPUT_DIR}/${out_file}" ]]; then
        echo "  [SKIP] ${out_file} already exists"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    echo "  [DOWNLOAD] ${remote_name} -> ${out_file}..."
    if wget -q --show-progress -O "${OUTPUT_DIR}/${out_file}" "${submitted_ftp}"; then
        DOWNLOADED=$((DOWNLOADED + 1))
    else
        echo "  [ERROR] Failed to download ${out_file}"
        [[ -f "${OUTPUT_DIR}/${out_file}" ]] && rm -f "${OUTPUT_DIR}/${out_file}"
        ERRORS=$((ERRORS + 1))
    fi

done < <(awk -F'\t' 'NR > 1 && $3 !~ /;/' "${REPORT_PATH}")

echo ""
echo "Download Summary"
echo "Files downloaded: ${DOWNLOADED}"
echo "Files skipped: ${SKIPPED}"
echo "Errors: ${ERRORS}"
echo ""

# Verify downloads
TOTAL_FILES=$(find "${OUTPUT_DIR}" -name "*.fastq.gz" 2>/dev/null | wc -l)
echo "Total FASTQ files: ${TOTAL_FILES}"

if [[ ${ERRORS} -gt 0 ]]; then
    echo ""
    echo "WARNING: ${ERRORS} errors occurred during download. Please review the output above."
    exit 1
fi
