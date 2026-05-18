# Whole-Exome Sequencing (WES) Analysis Pipeline

A reproducible, end-to-end pipeline for Whole-Exome Sequencing analysis, from raw SRA reads through aligned BAM, base recalibration, variant calling, filtering, and functional annotation.

---

## Background

The human genome contains approximately 3 billion base pairs, but only ~1–2% of that DNA encodes proteins. This coding region — the **exome** — is where ~85% of known disease-causing mutations occur. Whole-Exome Sequencing (WES) captures and sequences only these protein-coding regions, giving:

- **Clinical relevance**: Most pathogenic variants are in coding regions
- **Cost efficiency**: Sequencing ~1% of the genome at significantly lower cost
- **Higher depth**: For the same budget, WES achieves 100× coverage vs ~30× for WGS
- **Manageable data size**: WES files are far smaller than whole-genome data

This pipeline implements **GATK Best Practices** for germline short variant discovery on GRCh38.

---

## Pipeline Overview

```
Raw FASTQ (SRA)
    │
    ▼
Pre-alignment QC        ← FastQC + MultiQC
    │
    ▼
Adapter trimming        ← fastp
    │
    ▼
Post-trim QC            ← FastQC + MultiQC
    │
    ▼
Alignment               ← BWA-MEM → sorted BAM
    │
    ▼
Mark Duplicates         ← Picard MarkDuplicates
    │
    ▼
Base Recalibration      ← GATK BQSR (BaseRecalibrator + ApplyBQSR)
    │
    ▼
Coverage & QC           ← samtools depth, flagstat, CollectHsMetrics
    │
    ▼
Variant Calling         ← GATK HaplotypeCaller → GVCF
    │
    ▼
Genotyping              ← GATK GenotypeGVCFs → raw VCF
    │
    ▼
Variant Filtering       ← Hard-filter SNPs & INDELs → filtered VCF
    │
    ▼
Variant Annotation      ← Ensembl VEP → annotated VCF
    │
    ▼
Checksums               ← MD5 for all key outputs
```

---

## Repository Structure

```
Whole-Exome-Analysis/
├── .github/
│   └── workflows/
│       └── ci.yml              # GitHub Actions CI (lint, dry-run, integration test)
├── config/
│   └── config.yaml             # All parameters: paths, threads, tool settings
├── envs/
│   └── environment.yml         # Pinned conda environment with all tool versions
├── workflow/
│   └── Snakefile               # Snakemake pipeline rules (11 steps, 13 rules)
├── scripts/
│   ├── fix_bed_chr_names.awk   # Rename NCBI chr names to UCSC style in BED files
│   ├── coverage_summary.py     # Multi-threshold coverage report generator
│   └── prepare_test_data.sh    # Download chr20 subset for CI/testing
├── test_data/
│   ├── config_test.yaml        # Test config pointing at chr20 subset
│   ├── reads/                  # test_1.fastq.gz, test_2.fastq.gz (chr20 subset)
│   ├── ref/                    # chr20.fa + indices
│   └── bed/                    # chr20_capture.bed
├── resources/                  # (gitignored large files — download separately)
│   ├── GRCh38/                 # Reference FASTA + indices
│   ├── beds/                   # Capture BED files
│   └── known_sites/            # dbSNP, Mills, 1000G VCFs for BQSR
├── results/                    # (gitignored) Pipeline outputs
├── logs/                       # (gitignored) Per-rule log files
├── Dockerfile                  # Container for fully reproducible runs
├── .gitignore
└── README.md
```

---

## Prerequisites

- [Conda / Mamba](https://github.com/conda-forge/miniforge) (recommended: Mambaforge)
- [Snakemake ≥ 8.0](https://snakemake.readthedocs.io/)
- Docker (optional, for container-based runs)

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/suraj-chauhan-21/Whole-Exome-Analysis-.git
cd Whole-Exome-Analysis-
```

### 2. Create the conda environment

```bash
mamba env create -f envs/environment.yml
conda activate wes-pipeline
```

### 3. Download reference resources

```bash
# Reference genome (GRCh38)
mkdir -p resources/GRCh38
datasets download genome accession GCF_000001405.40 \
    --include genome,gff3 --filename resources/GRCh38/human_GRCh38.zip
unzip resources/GRCh38/human_GRCh38.zip -d resources/GRCh38/

# Index the reference
bwa index resources/GRCh38/GRCh38.p14.fa
samtools faidx resources/GRCh38/GRCh38.p14.fa
gatk CreateSequenceDictionary -R resources/GRCh38/GRCh38.p14.fa

# GATK resource bundle (known variant sites for BQSR)
# Download from: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0
mkdir -p resources/known_sites
# Place dbsnp_146.hg38.vcf.gz, Mills_and_1000G_gold_standard.indels.hg38.vcf.gz,
# and Homo_sapiens_assembly38.known_indels.vcf.gz into resources/known_sites/
```

### 4. Set up capture BED file

```bash
# Download Agilent SureSelect V6 BED from:
# https://earray.chem.agilent.com/suredesign/
mkdir -p resources/beds

# Fix chromosome names (NCBI → UCSC format)
awk -f scripts/fix_bed_chr_names.awk \
    resources/beds/chr_map.txt \
    resources/beds/Agilent_SureSelect_V6.bed \
    > resources/beds/Agilent_SureSelect_V6_chr_fixed.bed
```

### 5. Edit the config

Open `config/config.yaml` and update:
- `samples`: your SRA accession IDs
- `reference.fasta`: path to your GRCh38 FASTA
- `capture_bed_fixed`: path to your fixed BED file
- `known_sites.*`: paths to BQSR VCFs
- `threads.*` and `memory.picard_jvm` to match your machine

### 6. Run the pipeline

```bash
# Dry-run first (no jobs executed)
snakemake --cores 8 --use-conda --dry-run

# Full run
snakemake --cores 8 --use-conda

# On a cluster (SLURM example)
snakemake --cores 8 --use-conda \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "sbatch --mem=32G --cpus-per-task={threads}"
```

---

## Running Tests (chr20 subset)

A small chr20 test dataset is included for quick validation:

```bash
# Prepare test data (one-time setup, ~5 min)
bash scripts/prepare_test_data.sh

# Run pipeline on test data
snakemake --snakefile workflow/Snakefile \
          --configfile test_data/config_test.yaml \
          --cores 4 --use-conda

# Expected outputs in test_results/
ls test_results/metrics/
```

---

## Running with Docker

```bash
# Build the image
docker build -t wes-pipeline:1.0 .

# Dry-run
docker run --rm \
    -v $(pwd):/workspace -w /workspace \
    wes-pipeline:1.0 \
    snakemake --cores 4 --dry-run --configfile test_data/config_test.yaml

# Full run
docker run --rm \
    -v $(pwd):/workspace -w /workspace \
    wes-pipeline:1.0 \
    snakemake --cores 8
```

---

## Key Outputs

| File | Description |
|------|-------------|
| `results/qc/raw_multiqc/multiqc_report.html` | Pre-trim QC summary |
| `results/qc/trimmed_multiqc/multiqc_report.html` | Post-trim QC summary |
| `results/aligned/{sample}.dedup.bam` | Deduplicated sorted BAM |
| `results/bqsr/{sample}.bqsr.bam` | Base-recalibrated BAM |
| `results/metrics/{sample}.hs_metrics.txt` | WES capture efficiency metrics |
| `results/metrics/{sample}.coverage_summary.tsv` | Coverage at 1×, 10×, 20×, 30×, 50×, 100× |
| `results/metrics/{sample}.flagstat.txt` | Alignment rate and read statistics |
| `results/variants/{sample}.filtered.vcf.gz` | Hard-filtered SNP + INDEL calls |
| `results/annotated/{sample}.vep.vcf.gz` | VEP-annotated variant calls |
| `results/metrics/{sample}.md5sums.txt` | Checksums for output verification |

---

## Tool Versions

All versions are pinned in `envs/environment.yml`. Key tools:

| Tool | Version | Purpose |
|------|---------|---------|
| bwa | 0.7.17 | Read alignment |
| samtools | 1.19.2 | BAM manipulation |
| fastp | 0.23.4 | Adapter trimming |
| fastqc | 0.12.1 | Read QC |
| multiqc | 1.21 | Aggregate QC report |
| picard | 3.1.1 | Duplicate marking |
| gatk4 | 4.5.0.0 | BQSR, variant calling, filtering |
| ensembl-vep | 111.0 | Variant annotation |
| snakemake | 8.5.3 | Workflow management |

---

## Sample Data

This pipeline was developed and validated using **SRR22317682** (Illumina NovaSeq WES data, GEO accession GSE220768).

Results summary for SRR22317682:

| Metric | Value |
|--------|-------|
| Read pairs | 52,928,161 |
| % duplication | 0.56% |
| Estimated library size | ~4.7 billion molecules |
| Covered bases (≥1×) | 60,197,453 |
| Mean depth | 123× |

> Due to large file sizes (FASTQ, BAM, reference genome), only the pipeline code and documentation are stored in this repository. Raw data is available from NCBI SRA under accession SRR22317682. Output files are available upon request.

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **BWA**: Li & Durbin, *Bioinformatics* 2009
- **GATK**: Van der Auwera et al., *Current Protocols in Bioinformatics* 2013
- **Picard**: Broad Institute, https://broadinstitute.github.io/picard/
- **VEP**: McLaren et al., *Genome Biology* 2016
- **Snakemake**: Mölder et al., *F1000Research* 2021

---

## License

MIT License. See `LICENSE` for details.
