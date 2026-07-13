#!/usr/bin/env python3
"""Sort primer orientations, filter read pairs by polyA signal, and trim read1."""

from __future__ import annotations

import argparse
import csv
import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq


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


def core_id(seq_id: str) -> str:
    return seq_id.split("/")[0]


def primer_type(record) -> str:
    match = re.search(r"\|PRIMER=([^_\s|]+)", record.description)
    return match.group(1) if match else "unknown"


def max_consecutive_a(sequence: str) -> int:
    runs = re.findall(r"A+", sequence)
    return max(map(len, runs), default=0)


def has_polya_signal(sequence: str, min_run: int, search_bp: int) -> bool:
    return re.search(f"A{{{min_run},}}", sequence[:search_bp]) is not None


def trim_after_last_polya(sequence: str, qualities: list[int], min_run: int, search_bp: int) -> tuple[str, list[int]]:
    matches = list(re.finditer(f"A{{{min_run},}}", sequence[:search_bp]))
    if not matches:
        return sequence, qualities
    trim_pos = matches[-1].end()
    return sequence[trim_pos:], qualities[trim_pos:]


def write_histogram(values: list[int], output_pdf: Path, title: str) -> None:
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(8, 5))
    bins = range(0, 41)
    plt.hist(values, bins=bins, color="#4C78A8", edgecolor="black", align="left")
    plt.title(title)
    plt.xlabel("Maximum consecutive A bases")
    plt.ylabel("Read count")
    plt.xticks(range(0, 41, 5))
    plt.xlim(0, 40)
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_pdf)
    plt.close()


def sort_by_primer_orientation(read1_path: Path, read2_path: Path, output_prefix: Path) -> dict[str, Path]:
    """Split paired reads according to Fprimer/Rprimer orientation tags."""
    reads1 = {core_id(record.id): record for record in SeqIO.parse(read1_path, "fastq")}
    reads2 = {core_id(record.id): record for record in SeqIO.parse(read2_path, "fastq")}

    outputs = {
        "rf_r1": output_prefix.with_name(output_prefix.name + "_PRIMER_R_F_read1.fastq"),
        "rf_r2": output_prefix.with_name(output_prefix.name + "_PRIMER_R_F_read2.fastq"),
        "fr_r1": output_prefix.with_name(output_prefix.name + "_PRIMER_F_R_read1.fastq"),
        "fr_r2": output_prefix.with_name(output_prefix.name + "_PRIMER_F_R_read2.fastq"),
        "both_f": output_prefix.with_name(output_prefix.name + "_both_F.fastq"),
        "both_r": output_prefix.with_name(output_prefix.name + "_both_R.fastq"),
        "no_match": output_prefix.with_name(output_prefix.name + "_no_match.fastq"),
    }

    # Keep the original handle mapping from the analysis notes. The R_F output
    # handles look reversed by name, but the downstream file pairing depends on
    # this exact write order.
    with open(outputs["rf_r1"], "w") as out_R_F2, \
         open(outputs["rf_r2"], "w") as out_R_F1, \
         open(outputs["fr_r1"], "w") as out_F_R1, \
         open(outputs["fr_r2"], "w") as out_F_R2, \
         open(outputs["both_f"], "w") as out_both_F, \
         open(outputs["both_r"], "w") as out_both_R, \
         open(outputs["no_match"], "w") as out_no_match:
        for seq_id, record1 in reads1.items():
            record2 = reads2.get(seq_id)
            if record2 is None:
                SeqIO.write(record1, out_no_match, "fastq")
                continue

            p1 = primer_type(record1)
            p2 = primer_type(record2)
            if p1 == "Rprimer" and p2 == "Fprimer":
                SeqIO.write(record1, out_R_F1, "fastq")
                SeqIO.write(record2, out_R_F2, "fastq")
            elif p1 == "Fprimer" and p2 == "Rprimer":
                SeqIO.write(record1, out_F_R1, "fastq")
                SeqIO.write(record2, out_F_R2, "fastq")
            elif p1 == "Fprimer" and p2 == "Fprimer":
                SeqIO.write(record1, out_both_F, "fastq")
                SeqIO.write(record2, out_both_F, "fastq")
            elif p1 == "Rprimer" and p2 == "Rprimer":
                SeqIO.write(record1, out_both_R, "fastq")
                SeqIO.write(record2, out_both_R, "fastq")

        for seq_id, record2 in reads2.items():
            if seq_id not in reads1:
                SeqIO.write(record2, out_no_match, "fastq")

    return outputs


def filter_read1_polya(read1_path: Path, read2_path: Path, output_prefix: Path, min_run: int, search_bp: int) -> tuple[Path, Path]:
    """Keep read pairs whose read1 contains the polyA signal, then trim read1."""
    kept_ids: set[str] = set()
    all_counts: list[int] = []
    kept_counts: list[int] = []
    out_read1 = read1_path.with_name(read1_path.stem + "_poly8A.fastq")
    out_read2 = read2_path.with_name(read2_path.stem + "_poly8A.fastq")
    trimmed_read1 = read1_path.with_name(read1_path.stem + "_poly8A_trimmed.fastq")

    # The polyA signal is searched on read1. Only matching read IDs are kept
    # from read2, preserving the paired-read relationship.
    with out_read1.open("w") as handle:
        for record in SeqIO.parse(read1_path, "fastq"):
            sequence = str(record.seq)
            count = max_consecutive_a(sequence)
            all_counts.append(count)
            if has_polya_signal(sequence, min_run, search_bp):
                kept_ids.add(core_id(record.id))
                kept_counts.append(count)
                SeqIO.write(record, handle, "fastq")

    with out_read2.open("w") as handle:
        for record in SeqIO.parse(read2_path, "fastq"):
            if core_id(record.id) in kept_ids:
                SeqIO.write(record, handle, "fastq")

    # Trim read1 after the last qualifying polyA run within the search window;
    # this removes the polyA prefix before downstream read concatenation.
    with trimmed_read1.open("w") as handle:
        for record in SeqIO.parse(out_read1, "fastq"):
            sequence = str(record.seq)
            qualities = record.letter_annotations["phred_quality"]
            trimmed_sequence, trimmed_quality = trim_after_last_polya(sequence, qualities, min_run, search_bp)
            record.letter_annotations = {}
            record.seq = Seq(trimmed_sequence)
            record.letter_annotations["phred_quality"] = trimmed_quality
            SeqIO.write(record, handle, "fastq")

    write_histogram(all_counts, output_prefix.with_name(output_prefix.name + "_polya_all.pdf"), f"polyA distribution: {read1_path.name}")
    write_histogram(kept_counts, output_prefix.with_name(output_prefix.name + "_polya_kept.pdf"), f"polyA+ reads: {read1_path.name}")
    return trimmed_read1, out_read2


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", default="config/samples.tsv")
    parser.add_argument("--env", default="config/pipeline.env")
    args = parser.parse_args()

    env = load_env(Path(args.env))
    results_dir = Path(env.get("RESULTS_DIR", "results"))
    preprocess_dir = env.get("PREPROCESS_DIR", "preprocessing")
    min_run = int(env.get("POLYA_MIN_RUN", "8"))
    search_bp = int(env.get("POLYA_SEARCH_BP", "80"))

    for sample in read_samples(Path(args.samples)):
        sample_id = sample["sample_id"]
        sample_dir = results_dir / sample_id / preprocess_dir
        output_prefix = sample_dir / f"{sample_id}_output"
        read1_pass = sample_dir / (Path(sample["read1"]).name.replace(".gz", "").replace(".fastq", "").replace(".fq", "") + "_primers-pass.fastq")
        read2_pass = sample_dir / (Path(sample["read2"]).name.replace(".gz", "").replace(".fastq", "").replace(".fq", "") + "_primers-pass.fastq")

        orientation_files = sort_by_primer_orientation(read1_pass, read2_pass, output_prefix)
        for label in ("rf", "fr"):
            filter_read1_polya(
                orientation_files[f"{label}_r1"],
                orientation_files[f"{label}_r2"],
                output_prefix.with_name(f"{output_prefix.name}_{label}"),
                min_run,
                search_bp,
            )
        print(f"Finished primer sorting and polyA trimming for {sample_id}")


if __name__ == "__main__":
    main()
