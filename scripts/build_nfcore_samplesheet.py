#!/usr/bin/env python3
"""Create nf-core/rnaseq sample sheet from rnaseq.tsv metadata."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, type=Path, help="rnaseq.tsv metadata file")
    parser.add_argument("--fastq_dir", required=True, type=Path, help="Root FASTQ directory (containing colon/ and spleen/)")
    parser.add_argument("--output", required=True, type=Path, help="Output samplesheet CSV for nf-core")
    parser.add_argument(
        "--strandedness",
        default="auto",
        choices=["auto", "forward", "reverse", "unstranded"],
        help="Library strandedness reported to nf-core/rnaseq",
    )
    return parser.parse_args()


def detect_fastqs(sample_id: str, tissue: str, fastq_root: Path) -> tuple[Path, Path]:
    tissue_dir = fastq_root / tissue.lower()
    r1 = tissue_dir / f"{sample_id}_R1.fastq.gz"
    r2 = tissue_dir / f"{sample_id}_R2.fastq.gz"
    return r1, r2


def main() -> None:
    args = parse_args()
    rows = []
    missing_fastqs: list[str] = []
    with args.metadata.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sample_id = row.get("sample_id")
            tissue = (row.get("tissue") or "").strip().lower()
            if not sample_id or not tissue:
                continue
            fastq1, fastq2 = detect_fastqs(sample_id, tissue, args.fastq_dir)
            if not fastq1.exists() or not fastq2.exists():
                missing_fastqs.append(sample_id)
                continue
            rows.append(
                {
                    "sample": sample_id,
                    "fastq_1": str(fastq1),
                    "fastq_2": str(fastq2),
                    "strandedness": args.strandedness,
                }
            )

    if not rows:
        raise SystemExit("No samples found; ensure metadata includes sample_id and tissue columns.")
    if missing_fastqs:
        missing = ", ".join(sorted(missing_fastqs))
        raise SystemExit(f"Missing FASTQ pairs for samples: {missing}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "fastq_1", "fastq_2", "strandedness"])
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
