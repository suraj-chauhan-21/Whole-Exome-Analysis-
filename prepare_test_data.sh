#!/usr/bin/env bash
# =============================================================================
# prepare_test_data.sh
# =============================================================================
# Downloads a tiny chr20 subset of SRR22317682 and the matching reference
# slice for use in CI / quick validation runs.
#
# Runtime: ~3-5 min on a standard connection
# Disk usage: ~500 MB
#
# Usage:
#   bash scripts/prepare_test_data.sh
#
# Outputs (under test_data/):
#   test_data/reads/test_1.fastq.gz
#   test_data/reads/test_2.fastq.gz
#   test_data/ref/chr20.fa
#   test_data/ref/chr20.fa.fai
#   test_data/bed/chr20_capture.bed
# =============================================================================

set -euo pipefail

OUTDIR="test_data"
READS_DIR="${OUTDIR}/reads"
REF_DIR="${OUTDIR}/ref"
BED_DIR="${OUTDIR}/bed"

mkdir -p "${READS_DIR}" "${REF_DIR}" "${BED_DIR}"

echo "=== Step 1: Download chr20 FASTA from NCBI ==="
# GRCh38 chr20 (NC_000020.11), ~64 Mb compressed
datasets download genome accession GCF_000001405.40 \
    --chromosomes 20 \
    --include genome \
    --filename "${REF_DIR}/chr20.zip"

unzip -o "${REF_DIR}/chr20.zip" -d "${REF_DIR}/"
# The FASTA will land inside ncbi_dataset/data/...
NCBI_FA=$(find "${REF_DIR}" -name "*.fna" | head -1)

# Rename to a simple path
cp "${NCBI_FA}" "${REF_DIR}/chr20.fa"
samtools faidx "${REF_DIR}/chr20.fa"

echo "=== Step 2: Build BWA index on chr20 FASTA ==="
bwa index "${REF_DIR}/chr20.fa"

echo "=== Step 3: Subset SRR22317682 reads mapping to chr20 ==="
# Fetch a small portion of the SRA run, align to chr20, extract paired reads
# We use a partial download (first 500k read pairs) for speed
fastq-dump SRR22317682 \
    --split-files \
    --maxSpotId 500000 \
    --outdir /tmp/srr_subset/

# Align to chr20 only
bwa mem -t 4 "${REF_DIR}/chr20.fa" \
    /tmp/srr_subset/SRR22317682_1.fastq \
    /tmp/srr_subset/SRR22317682_2.fastq \
    | samtools view -b -F 4 - \
    | samtools sort -o /tmp/chr20_subset.bam

samtools index /tmp/chr20_subset.bam

# Extract FASTQs from chr20-mapped reads
samtools fastq -1 "${READS_DIR}/test_1.fastq" \
               -2 "${READS_DIR}/test_2.fastq" \
               -0 /dev/null -s /dev/null \
               /tmp/chr20_subset.bam

# Compress
pigz -f "${READS_DIR}/test_1.fastq"
pigz -f "${READS_DIR}/test_2.fastq"

echo "=== Step 4: Create a chr20 capture BED (subset of Agilent V6) ==="
# If the full BED is already present, subset it; otherwise create a minimal one
FULL_BED="resources/beds/Agilent_SureSelect_V6_chr_fixed.bed"
if [[ -f "${FULL_BED}" ]]; then
    grep "^chr20" "${FULL_BED}" > "${BED_DIR}/chr20_capture.bed"
    echo "  Subsetted ${FULL_BED} → ${BED_DIR}/chr20_capture.bed"
else
    echo "  Full BED not found — creating a small synthetic chr20 BED"
    cat > "${BED_DIR}/chr20_capture.bed" << 'EOF'
chr20	34590000	34620000
chr20	44000000	44030000
chr20	55000000	55030000
EOF
fi

echo "=== Step 5: Create test config ==="
cat > "${OUTDIR}/config_test.yaml" << EOF
# Test configuration — uses chr20 subset only
samples:
  - test

reference:
  fasta: "test_data/ref/chr20.fa"
  fai:   "test_data/ref/chr20.fa.fai"
  dict:  "test_data/ref/chr20.dict"

capture_bed_fixed: "test_data/bed/chr20_capture.bed"

# No BQSR in test mode (known sites not included)
known_sites:
  dbsnp:  ""
  mills:  ""
  indels: ""

threads:
  bwa: 4
  samtools: 2
  fastp: 2
  gatk: 2
  fastqc: 2

memory:
  picard_jvm: "4g"

fastp:
  min_length: 36
  qualified_quality: 20
  unqualified_pct: 40

gatk:
  haplotypecaller_mode: "EMIT_VARIANTS_ONLY"
  stand_call_conf: 10

hard_filter:
  snp:   {QD: 2.0, MQ: 40.0, FS: 60.0, SOR: 3.0, MQRankSum: -12.5, ReadPosRankSum: -8.0}
  indel: {QD: 2.0, FS: 200.0, SOR: 10.0, ReadPosRankSum: -20.0}

downsample_fraction: "123.0813"
depth_thresholds: [1, 10, 20, 30]

dirs:
  data:     "test_data/reads"
  aligned:  "test_results/aligned"
  qc:       "test_results/qc"
  metrics:  "test_results/metrics"
  bqsr:     "test_results/bqsr"
  variants: "test_results/variants"
  annotated: "test_results/annotated"
  logs:     "test_results/logs"
EOF

echo ""
echo "=== Test data ready ==="
echo "  Reads  : ${READS_DIR}/test_1.fastq.gz"
echo "          : ${READS_DIR}/test_2.fastq.gz"
echo "  Ref    : ${REF_DIR}/chr20.fa"
echo "  BED    : ${BED_DIR}/chr20_capture.bed"
echo "  Config : ${OUTDIR}/config_test.yaml"
echo ""
echo "Run the pipeline on test data with:"
echo "  snakemake --cores 4 --use-conda --configfile test_data/config_test.yaml"
