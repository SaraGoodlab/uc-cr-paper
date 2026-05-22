#!/usr/bin/env python3
# Create the nf-core/rnaseq samplesheet for this dataset.

import csv
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
METADATA = REPO_ROOT / "data" / "metadata" / "rnaseq.tsv"
FASTQ_DIR = REPO_ROOT / "data" / "raw" / "rnaseq"
OUTPUT = REPO_ROOT / "data" / "metadata" / "nfcore_samplesheet.csv"
STRANDEDNESS = "auto"


def main() -> None:
    rows = []
    missing_fastqs = []

    with METADATA.open("r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample_id = row["sample_id"].strip()
            tissue = row["tissue"].strip().lower()
            fastq1 = FASTQ_DIR / tissue / f"{sample_id}_R1.fastq.gz"
            fastq2 = FASTQ_DIR / tissue / f"{sample_id}_R2.fastq.gz"

            if not fastq1.exists() or not fastq2.exists():
                missing_fastqs.append(sample_id)
                continue

            rows.append({
                "sample": sample_id,
                "fastq_1": str(fastq1),
                "fastq_2": str(fastq2),
                "strandedness": STRANDEDNESS,
            })

    if missing_fastqs:
        raise SystemExit("Missing FASTQ pairs for samples: " + ", ".join(sorted(missing_fastqs)))
    if not rows:
        raise SystemExit("No RNA-seq samples were written")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "fastq_1", "fastq_2", "strandedness"])
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
