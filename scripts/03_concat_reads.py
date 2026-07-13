#!/usr/bin/env python3
"""Merge primer-orientation FASTQs and concatenate read1 with reverse-complemented read2."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def load_env(path: Path) -> dict[str, str]:
    env: dict[str, str] = {}
    if not path.exists():
        return env
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        env[key] = value.strip().strip('"')
    return env


def read_samples(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def concatenate_reads(read1_path: Path, read2_path: Path, output_path: Path) -> None:
    """Append reverse-complemented read2 to read1 and preserve quality scores."""
    with read1_path.open() as handle_r1, read2_path.open() as handle_r2, output_path.open("w") as handle_out:
        read1_records = SeqIO.parse(handle_r1, "fastq")
        read2_records = SeqIO.parse(handle_r2, "fastq")
        for r1, r2 in zip(read1_records, read2_records):
            rc_r2_seq = r2.seq.reverse_complement()
            concatenated_seq = r1.seq + rc_r2_seq
            concatenated_qual = r1.letter_annotations["phred_quality"] + r2.letter_annotations["phred_quality"][::-1]
            record = SeqRecord(
                concatenated_seq,
                id=r1.id,
                description=r1.description,
                letter_annotations={"phred_quality": concatenated_qual},
            )
            SeqIO.write(record, handle_out, "fastq")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", default="config/samples.tsv")
    parser.add_argument("--env", default="config/pipeline.env")
    args = parser.parse_args()

    env = load_env(Path(args.env))
    results_dir = Path(env.get("RESULTS_DIR", "results"))
    preprocess_dir = env.get("PREPROCESS_DIR", "preprocessing")

    for sample in read_samples(Path(args.samples)):
        sample_id = sample["sample_id"]
        sample_dir = results_dir / sample_id / preprocess_dir
        prefix = sample_dir / f"{sample_id}_output"

        read1_paths = [
            prefix.with_name(prefix.name + "_PRIMER_F_R_read1_poly8A_trimmed.fastq"),
            prefix.with_name(prefix.name + "_PRIMER_R_F_read1_poly8A_trimmed.fastq"),
        ]
        read2_paths = [
            prefix.with_name(prefix.name + "_PRIMER_F_R_read2_poly8A.fastq"),
            prefix.with_name(prefix.name + "_PRIMER_R_F_read2_poly8A.fastq"),
        ]
        combined_read1 = sample_dir / f"{sample_id}_trimpolyA_combined_read1.fastq"
        combined_read2 = sample_dir / f"{sample_id}_trimpolyA_combined_read2.fastq"

        # Merge both valid primer-orientation groups before concatenation.
        with combined_read1.open("w") as out_handle:
            for path in read1_paths:
                out_handle.write(path.read_text())
        with combined_read2.open("w") as out_handle:
            for path in read2_paths:
                out_handle.write(path.read_text())

        output_path = sample_dir / f"{sample_id}_concatenated.fastq"
        concatenate_reads(combined_read1, combined_read2, output_path)
        print(f"Saved concatenated reads for {sample_id}: {output_path}")


if __name__ == "__main__":
    main()
