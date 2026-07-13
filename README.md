# SpaCir

SpaCir contains scripts used to preprocess spatial immune-receptor sequencing data and generate barcode-level clonotype summaries.

The current repository version includes the upstream workflow:

1. primer and quality filtering with pRESTO
2. primer-orientation sorting and polyA trimming
3. paired-read merging and read concatenation
4. MiXCR annotation
5. MiXCR alignment filtering
6. spatial barcode matching, UMI filtering, clonotype matrix generation, and VDJC summaries

Downstream plotting scripts are organized under `analysis/figure1/`, `analysis/figure2/`, `analysis/figure3/`, `analysis/figure4/`, and `analysis/figure5/`.

## Repository Layout

```text
SpaCir/
├── config/
│   ├── pipeline.env.example
│   └── samples.tsv
├── primers/
│   ├── IGH_primer_all.fasta
│   └── TRB_primer_all.fasta
├── scripts/
│   ├── 01_filter_primers.sh
│   ├── 02_sort_trim_polya.py
│   ├── 03_concat_reads.py
│   ├── 04_run_mixcr.sh
│   ├── 05_filter_alignments.sh
│   └── 06_spatial_barcode_umi_pipeline.sh
├── docs/
│   └── workflow.md
├── analysis/
│   ├── figure1/
│   │   ├── R/
│   │   └── config/
│   ├── figure2/
│   │   ├── R/
│   │   ├── python/
│   │   └── config/
│   ├── figure3/
│   │   ├── R/
│   │   └── config/
│   ├── figure4/
│   │   ├── R/
│   │   └── config/
│   └── figure5/
│       ├── R/
│       └── config/
├── notebooks/
├── requirements.txt
└── .gitignore
```

## Quick Start

Copy the example configuration and edit it for your local environment:

```bash
cp config/pipeline.env.example config/pipeline.env
```

Edit `config/samples.tsv` to list sample IDs, receptor chains, paired FASTQ files, and primer FASTA files.
The repository includes the primer files used by the workflow:

- `primers/IGH_primer_all.fasta` for IGH libraries
- `primers/TRB_primer_all.fasta` for TRB libraries

MiXCR is run with the manuscript workflow: `generic-amplicon`, `VEnd` to `C` floating boundaries, and saved original reads.

Run the preprocessing steps from the repository root:

```bash
bash scripts/01_filter_primers.sh config/samples.tsv config/pipeline.env
python scripts/02_sort_trim_polya.py --samples config/samples.tsv --env config/pipeline.env
python scripts/03_concat_reads.py --samples config/samples.tsv --env config/pipeline.env
bash scripts/04_run_mixcr.sh config/samples.tsv config/pipeline.env
bash scripts/05_filter_alignments.sh config/samples.tsv config/pipeline.env
```

Run the spatial barcode and UMI workflow for one sample:

```bash
BASE_DIR=results SAMPLE_GLOB='*-TRB' CHUNKSIZE=50000 MAX_WORKERS=4 \
  bash scripts/06_spatial_barcode_umi_pipeline.sh
```

Figure 1 and Figure S1 plotting scripts are organized by panel:

```bash
ls analysis/figure1/R
```

See `analysis/figure1/README.md` and `analysis/figure1/config/panel_map.csv` for the mapping between scripts and panels.

Figure 2 and Figure S2-S4 plotting scripts are organized the same way:

```bash
ls analysis/figure2/R
```

See `analysis/figure2/README.md` and `analysis/figure2/config/panel_map.csv` for the mapping between scripts and panels.

Figure 3 and Figure S5-S7 plotting scripts are organized by panel:

```bash
ls analysis/figure3/R
```

See `analysis/figure3/README.md` and `analysis/figure3/config/panel_map.csv` for the mapping between scripts and panels.

Figure 4 and Figure S8-S10 plotting scripts are organized by panel:

```bash
ls analysis/figure4/R
```

See `analysis/figure4/README.md` and `analysis/figure4/config/panel_map.csv` for the mapping between scripts and panels.

Figure 5 and Figure S11-S12 plotting scripts are organized by panel:

```bash
ls analysis/figure5/R
```

See `analysis/figure5/README.md` and `analysis/figure5/config/panel_map.csv` for the mapping between scripts and panels.

## Notes

- Local absolute paths and original sample identifiers were intentionally removed.
- Primer FASTA files are included in `primers/` so pRESTO can be run without private local paths.
- Large alignment files are processed in chunks in `05_filter_alignments.sh` and `06_spatial_barcode_umi_pipeline.sh` to reduce memory use.
- The default UMI threshold is `3`, and the default Levenshtein distance threshold is `1`, matching the current analysis workflow.
- Intermediate preprocessing files are written under `results/<sample_id>/preprocessing` by default. Set `PREPROCESS_DIR` in `config/pipeline.env` to use a different directory name.
