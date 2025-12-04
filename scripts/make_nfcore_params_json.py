#!/usr/bin/env python3
"""Apply runtime overrides to an nf-core/rnaseq params JSON template."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", required=True, type=Path, help="Base nf-core params JSON")
    parser.add_argument("--output", required=True, type=Path, help="Path to write rendered params JSON")
    parser.add_argument("--samplesheet", required=True, type=Path, help="nf-core samplesheet CSV")
    parser.add_argument("--outdir", required=True, type=Path, help="nf-core output directory")
    parser.add_argument("--fasta", required=True, type=Path, help="Reference genome FASTA")
    parser.add_argument("--gtf", required=True, type=Path, help="Reference annotation GTF")
    parser.add_argument("--transcriptome", required=True, type=Path, help="Reference transcriptome FASTA")
    parser.add_argument("--aligner", required=True, help="Aligner value to set (e.g., star_rsem)")
    return parser.parse_args()


def normalize(path_value: Path) -> str:
    return str(path_value.expanduser().resolve(strict=False))


def main() -> None:
    args = parse_args()

    with args.template.open("r", encoding="utf-8") as handle:
        params = json.load(handle)

    overrides = {
        "input": normalize(args.samplesheet),
        "outdir": normalize(args.outdir),
        "fasta": normalize(args.fasta),
        "gtf": normalize(args.gtf),
        "transcript_fasta": normalize(args.transcriptome),
        "aligner": args.aligner,
    }

    params.update({k: v for k, v in overrides.items() if v})

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        json.dump(params, handle, indent=4)
        handle.write("\n")


if __name__ == "__main__":
    main()
