# SpaCir Workflow

## 1. Primer And Quality Filtering

`scripts/01_filter_primers.sh` runs pRESTO quality filtering and primer masking for each row in `config/samples.tsv`.
The `primer_fasta` column should point to one of the included primer files, or to a user-supplied replacement:

- `primers/IGH_primer_all.fasta`
- `primers/TRB_primer_all.fasta`

Key outputs:

- `{sample_id}_1_quality-pass.fastq`
- `{sample_id}_2_quality-pass.fastq`
- `{sample_id}_1_primers-pass.fastq`
- `{sample_id}_2_primers-pass.fastq`

## 2. Primer Orientation Sorting And PolyA Trimming

`scripts/02_sort_trim_polya.py` sorts reads by primer orientation, filters read pairs where read1 contains a polyA signal, trims read1 after the last polyA run, and writes orientation-specific FASTQ files.

## 3. Read Concatenation

`scripts/03_concat_reads.py` combines primer-orientation outputs and concatenates read1 with the reverse complement of read2.

## 4. MiXCR Annotation

`scripts/04_run_mixcr.sh` runs MiXCR alignment, assembly, and alignment export.
The retained manuscript workflow uses `generic-amplicon`, `--floating-left-alignment-boundary VEnd`, and `--floating-right-alignment-boundary C`.

## 5. Alignment Filtering

`scripts/05_filter_alignments.sh` filters exported MiXCR alignments by chain and clone assignment:

- keeps rows where `chains` starts with the sample chain (`TRB` or `IGH`)
- removes rows where `cloneId == -1`
- writes `{sample_id}-filtered_alignments_poly8A.tsv`

## 6. Spatial Barcode And UMI Processing

`scripts/06_spatial_barcode_umi_pipeline.sh` matches spatial barcodes from MiXCR alignments, extracts UMIs, filters by Levenshtein distance and UMI count, then produces clonotype, matrix, diversity, mutation-rate, and VDJC summary files.

The large MiXCR alignment TSV is read in chunks to reduce memory use.
