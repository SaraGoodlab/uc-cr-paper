#!/bin/bash
# Download RNA-seq FASTQ files from European Nucleotide Archive (ENA) using a sample mapping file

set -euo pipefail 

# Configuration
ENA_PROJECT="${ENA_PROJECT:-PRJEB102830}"
OUTPUT_DIR="${OUTPUT_DIR:-../../data/raw/rnaseq}"
SAMPLE_MAPPING="${SAMPLE_MAPPING:-../../data/metadata/rnaseq.tsv}"  # Sample mapping file

# Create output directories
COLON_DIR="${OUTPUT_DIR}/colon"
SPLEEN_DIR="${OUTPUT_DIR}/spleen"
mkdir -p "${COLON_DIR}"
mkdir -p "${SPLEEN_DIR}"

echo "ENA RNA-seq Data Download"
echo "ENA Project: ${ENA_PROJECT}"
echo "Output directory: ${OUTPUT_DIR}"

# Download file report from ENA
REPORT_URL="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ENA_PROJECT}&result=read_run&fields=run_accession,sample_accession,submitted_ftp&format=tsv&download=true"
    
if [[ ! -f "../../data/metadata/ena_file_report.tsv" ]]; then
    echo "Fetching file report from ENA..."
    wget -q -O "../../data/metadata/ena_file_report.tsv" "${REPORT_URL}"
fi

if [[ ! -f "../../data/metadata/ena_file_report.tsv" ]]; then
    echo "Error: Failed to download ENA file report"
    exit 1
fi

# Create sample mappings
declare -A SAMPLE_TO_TISSUE
declare -A SAMPLE_TO_SAMPLE_ID

echo "Loading sample mapping from ${SAMPLE_MAPPING}..."

# Column positions (0-based indexing for array)
SAMPLE_ACC_COL=0
TISSUE_COL=1
TREATMENT_COL=2
SAMPLE_ID_COL=3

# Read mapping file
while IFS=$'\t' read -r -a FIELDS; do
    # Skip header
    [[ "${FIELDS[0]:-}" == "sample_accession" ]] && continue
    [[ -z "${FIELDS[0]:-}" ]] && continue
    
    # Check if we have enough fields
    if [[ ${#FIELDS[@]} -lt 4 ]]; then
        echo "Warning: Skipping line with insufficient fields (expected 4, got ${#FIELDS[@]})"
        continue
    fi
    
    # Extract fields (0-based array indexing, with default empty string if unset)
    sample_accession="${FIELDS[$SAMPLE_ACC_COL]:-}"
    tissue="${FIELDS[$TISSUE_COL]:-}"
    sample_id="${FIELDS[$SAMPLE_ID_COL]:-}"
    
    # Validate required fields
    if [[ -z "$sample_accession" ]] || [[ -z "$tissue" ]] || [[ -z "$sample_id" ]]; then
        echo "Warning: Skipping line with missing required fields (sample_accession, tissue, or sample_id)"
        continue
    fi
    
    # Normalize tissue name (lowercase, remove spaces)
    tissue=$(echo "$tissue" | tr '[:upper:]' '[:lower:]' | tr -d ' ')
    
    # Store mappings
    case "$tissue" in
        colon)
            SAMPLE_TO_TISSUE["$sample_accession"]="colon"
            SAMPLE_TO_SAMPLE_ID["$sample_accession"]="$sample_id"
            ;;
        spleen)
            SAMPLE_TO_TISSUE["$sample_accession"]="spleen"
            SAMPLE_TO_SAMPLE_ID["$sample_accession"]="$sample_id"
            ;;
        *)
            echo "Warning: Unknown tissue type '$tissue' for sample $sample_accession, skipping..."
            continue
            ;;
    esac
done < "${SAMPLE_MAPPING}"

# Parse file report and download FASTQ files
DOWNLOADED=0
SKIPPED=0
ERRORS=0

echo "Processing ENA file report..."
while IFS=$'\t' read -r run_accession sample_accession submitted_ftp; do
    # Skip header
    [[ "$run_accession" == "run_accession" ]] && continue
    [[ -z "$run_accession" ]] && continue
    
    # Determine tissue type
    TARGET_DIR=""
    tissue="${SAMPLE_TO_TISSUE[$sample_accession]:-}"
    if [[ "$tissue" == "colon" ]]; then
        TARGET_DIR="${COLON_DIR}"
    elif [[ "$tissue" == "spleen" ]]; then
        TARGET_DIR="${SPLEEN_DIR}"
    else
        echo "Warning: No tissue mapping found for sample ${sample_accession}, skipping..."
        SKIPPED=$((SKIPPED + 1))
        continue
    fi
    
    # Download paired-end files
    if [[ -z "$submitted_ftp" ]]; then
        echo "Warning: No FASTQ FTP URL for ${run_accession}, skipping..."
        SKIPPED=$((SKIPPED + 1))
        continue
    fi
    
    # Get sample_id for renaming
    sample_id="${SAMPLE_TO_SAMPLE_ID[$sample_accession]:-}"
    if [[ -z "$sample_id" ]]; then
        echo "Warning: No sample_id found for sample ${sample_accession}, using original filename"
        sample_id="${sample_accession}"
    fi
    
    echo "[${run_accession}] Downloading to ${TARGET_DIR} (sample_id: ${sample_id})..."
    IFS=';' read -ra FILES <<< "$submitted_ftp"
    
    # Determine read number for paired-end files
    for i in "${!FILES[@]}"; do
        file_url="${FILES[$i]}"
        original_filename=$(basename "${file_url}")
        temp_file_path="${TARGET_DIR}/.tmp_${original_filename}"
        final_filename=""
        read_num=""
        
        # Use index to assign read order (first file = R1, second = R2)
        read_num=$((i + 1))
        
        final_filename="${sample_id}_R${read_num}.fastq.gz"
        final_file_path="${TARGET_DIR}/${final_filename}"
        
        # Skip if file already exists
        if [[ -f "${final_file_path}" ]]; then
            echo "  [SKIP] ${final_filename} already exists"
            continue
        fi
        
        # Download file to temporary location
        echo "  [DOWNLOAD] ${original_filename} -> ${final_filename}..."
        if wget -q --show-progress -O "${temp_file_path}" "${file_url}"; then
            # Move to final location with new name
            mv "${temp_file_path}" "${final_file_path}"
            DOWNLOADED=$((DOWNLOADED + 1))
        else
            echo "  [ERROR] Failed to download ${original_filename}"
            [[ -f "${temp_file_path}" ]] && rm -f "${temp_file_path}"
            ERRORS=$((ERRORS + 1))
        fi
    done
done < "../../data/metadata/ena_file_report.tsv"

echo ""
echo "Download Summary"
echo "Files downloaded: ${DOWNLOADED}"
echo "Files skipped: ${SKIPPED}"
echo "Errors: ${ERRORS}"
echo ""
echo "Files are located in:"
echo "  Colon: ${COLON_DIR}"
echo "  Spleen: ${SPLEEN_DIR}"
echo ""

# Verify downloads
TOTAL_FILES=$(find "${OUTPUT_DIR}" -name "*.fastq.gz" 2>/dev/null | wc -l)
echo "Total FASTQ files: ${TOTAL_FILES}"

if [[ ${ERRORS} -gt 0 ]]; then
    echo ""
    echo "WARNING: ${ERRORS} errors occurred during download. Please review the output above."
    exit 1
fi
