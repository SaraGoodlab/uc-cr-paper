"""
Shared workflow constants.
"""

WORKFLOW_DIR = PROJECT_DIR / "workflow"
RAW_DATA_DIR = DATA_DIR / "raw"
PROCESSED_DATA_DIR = DATA_DIR / "processed"
EXTERNAL_REFERENCE_DIR = DATA_DIR / "external" / "reference"
METADATA_DIR = DATA_DIR / "metadata"
CONTAINERS_DIR = WORKFLOW_DIR / "containers"

R_ANALYSIS_CONTAINER = str(CONTAINERS_DIR / "r-analysis" / "r-analysis.sif")
QIIME2_CONTAINER = str(CONTAINERS_DIR / "microbiome-qiime2" / "qiime2.sif")
SHELL_PREAMBLE = "set -euo pipefail"
