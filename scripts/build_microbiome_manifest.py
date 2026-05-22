#!/usr/bin/env python3

# Generate the QIIME2 manifest for the 16S dataset.

import argparse
from pathlib import Path

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastq_dir", required=True, help="Directory containing {sample_id}.fastq.gz files")
    parser.add_argument("--output", required=True, help="Path to manifest TSV to write")
    return parser.parse_args()

def main() -> None:
  args = parse_args()
  fastq_dir = Path(args.fastq_dir)
  output_path = Path(args.output)
  fastqs = sorted(fastq_dir.glob("*.fastq.gz"))

    if not fastqs:
        raise SystemExit(f"No .fastq.gz files found in {fastq_dir}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("sample-id\tabsolute-filepath\tdirection\n")
        for fastq in fastqs:
            sample_id = fastq.name.removesuffix(".fastq.gz")
            handle.write(f"{sample_id}\t{fastq.resolve()}\tforward\n")

if __name__ == "__main__":
    main()
