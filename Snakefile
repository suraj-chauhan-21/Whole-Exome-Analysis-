"""
Whole-Exome Sequencing (WES) Analysis Pipeline
===============================================
Snakemake workflow covering:
  1.  SRA download
  2.  Pre-alignment QC  (FastQC + MultiQC)
  3.  Adapter trimming  (fastp)
  4.  Post-trim QC      (FastQC + MultiQC)
  5.  Alignment         (BWA-MEM)
  6.  Duplicate marking (Picard)
  7.  BQSR              (GATK BaseRecalibrator + ApplyBQSR)
  8.  Coverage & QC     (samtools depth, flagstat, CollectHsMetrics)
  9.  Variant calling   (GATK HaplotypeCaller → GenomicsDBImport → GenotypeGVCFs)
  10. Variant filtering (hard-filter SNP + INDEL)
  11. Variant annotation (VEP)

Usage:
  snakemake --cores 8 --use-conda
  snakemake --cores 8 --use-conda --configfile config/config.yaml
"""

import os

configfile: "config/config.yaml"

SAMPLES    = config["samples"]
REF        = config["reference"]["fasta"]
BED        = config["capture_bed_fixed"]
KNOWN      = config["known_sites"]
THREADS    = config["threads"]
DIRS       = config["dirs"]

# ─────────────────────────────────────────────────────────────────────────────
# Target rule — defines the final outputs Snakemake aims to produce
# ─────────────────────────────────────────────────────────────────────────────
rule all:
    input:
        # MultiQC reports
        expand("{qc}/raw_multiqc/multiqc_report.html",     qc=DIRS["qc"]),
        expand("{qc}/trimmed_multiqc/multiqc_report.html", qc=DIRS["qc"]),
        # Final annotated VCF
        expand("{ann}/{sample}.vep.vcf.gz", ann=DIRS["annotated"], sample=SAMPLES),
        # Coverage summary
        expand("{met}/{sample}.hs_metrics.txt", met=DIRS["metrics"], sample=SAMPLES),
        # Checksums
        expand("{met}/{sample}.md5sums.txt",    met=DIRS["metrics"], sample=SAMPLES),


# ─────────────────────────────────────────────────────────────────────────────
# RULE 1: Download SRA data
# ─────────────────────────────────────────────────────────────────────────────
rule sra_download:
    output:
        r1 = "{data}/{sample}/{sample}_1.fastq",
        r2 = "{data}/{sample}/{sample}_2.fastq",
    params:
        outdir = lambda wc: f"{DIRS['data']}/{wc.sample}",
    log:
        "logs/{sample}/sra_download.log"
    threads: 4
    shell:
        """
        mkdir -p {params.outdir}
        prefetch -O {params.outdir} --progress {wildcards.sample} 2>{log}
        fasterq-dump {params.outdir}/{wildcards.sample}.sra \
            --split-files --threads {threads} \
            --outdir {params.outdir} 2>>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# RULE 2a: FastQC on raw reads
# ─────────────────────────────────────────────────────────────────────────────
rule fastqc_raw:
    input:
        r1 = rules.sra_download.output.r1,
        r2 = rules.sra_download.output.r2,
    output:
        html_r1 = "{qc}/raw/{sample}_1_fastqc.html".format(qc=DIRS["qc"], sample="{sample}"),
        html_r2 = "{qc}/raw/{sample}_2_fastqc.html".format(qc=DIRS["qc"], sample="{sample}"),
        zip_r1  = "{qc}/raw/{sample}_1_fastqc.zip".format( qc=DIRS["qc"], sample="{sample}"),
        zip_r2  = "{qc}/raw/{sample}_2_fastqc.zip".format( qc=DIRS["qc"], sample="{sample}"),
    params:
        outdir = f"{DIRS['qc']}/raw",
    log:
        "logs/{sample}/fastqc_raw.log"
    threads: config["threads"]["fastqc"]
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2} 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# RULE 2b: MultiQC on raw FastQC reports
# ─────────────────────────────────────────────────────────────────────────────
rule multiqc_raw:
    input:
        expand("{qc}/raw/{sample}_{read}_fastqc.zip",
               qc=DIRS["qc"], sample=SAMPLES, read=["1","2"]),
    output:
        report = f"{DIRS['qc']}/raw_multiqc/multiqc_report.html",
    params:
        indir  = f"{DIRS['qc']}/raw",
        outdir = f"{DIRS['qc']}/raw_multiqc",
    log:
        "logs/multiqc_raw.log"
    shell:
        "multiqc {params.indir} -o {params.outdir} --force 2>{log}"


# ─────────────────────────────────────────────────────────────────────────────
# RULE 3: Adapter & quality trimming with fastp
# ─────────────────────────────────────────────────────────────────────────────
rule fastp_trim:
    input:
        r1 = rules.sra_download.output.r1,
        r2 = rules.sra_download.output.r2,
    output:
        r1      = f"{DIRS['data']}/{{sample}}/{{sample}}_1.trimmed.fastq.gz",
        r2      = f"{DIRS['data']}/{{sample}}/{{sample}}_2.trimmed.fastq.gz",
        html    = f"{DIRS['qc']}/fastp/{{sample}}_fastp.html",
        json    = f"{DIRS['qc']}/fastp/{{sample}}_fastp.json",
    params:
        min_len = config["fastp"]["min_length"],
        qual    = config["fastp"]["qualified_quality"],
        unqual  = config["fastp"]["unqualified_pct"],
    log:
        "logs/{sample}/fastp.log"
    threads: config["threads"]["fastp"]
    shell:
        """
        mkdir -p {DIRS[qc]}/fastp
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --html {output.html} --json {output.json} \
            --length_required {params.min_len} \
            --qualified_quality_phred {params.qual} \
            --unqualified_percent_limit {params.unqual} \
            --detect_adapter_for_pe \
            --thread {threads} 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# RULE 4a: FastQC on trimmed reads
# ─────────────────────────────────────────────────────────────────────────────
rule fastqc_trimmed:
    input:
        r1 = rules.fastp_trim.output.r1,
        r2 = rules.fastp_trim.output.r2,
    output:
        html_r1 = f"{DIRS['qc']}/trimmed/{{sample}}_1.trimmed_fastqc.html",
        html_r2 = f"{DIRS['qc']}/trimmed/{{sample}}_2.trimmed_fastqc.html",
        zip_r1  = f"{DIRS['qc']}/trimmed/{{sample}}_1.trimmed_fastqc.zip",
        zip_r2  = f"{DIRS['qc']}/trimmed/{{sample}}_2.trimmed_fastqc.zip",
    params:
        outdir = f"{DIRS['qc']}/trimmed",
    log:
        "logs/{sample}/fastqc_trimmed.log"
    threads: config["threads"]["fastqc"]
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2} 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# RULE 4b: MultiQC on trimmed FastQC + fastp JSON reports
# ─────────────────────────────────────────────────────────────────────────────
rule multiqc_trimmed:
    input:
        expand(f"{DIRS['qc']}/trimmed/{{sample}}_{{read}}.trimmed_fastqc.zip",
               sample=SAMPLES, read=["1","2"]),
        expand(f"{DIRS['qc']}/fastp/{{sample}}_fastp.json", sample=SAMPLES),
    output:
        report = f"{DIRS['qc']}/trimmed_multiqc/multiqc_report.html",
    params:
        indir  = f"{DIRS['qc']}",
        outdir = f"{DIRS['qc']}/trimmed_multiqc",
    log:
        "logs/multiqc_trimmed.log"
    shell:
        "multiqc {params.indir}/trimmed {params.indir}/fastp -o {params.outdir} --force 2>{log}"


# ─────────────────────────────────────────────────────────────────────────────
# RULE 5: Alignment with BWA-MEM
# ─────────────────────────────────────────────────────────────────────────────
rule bwa_align:
    input:
        r1  = rules.fastp_trim.output.r1,
        r2  = rules.fastp_trim.output.r2,
        ref = REF,
    output:
        bam = f"{DIRS['aligned']}/{{sample}}.sorted.bam",
        bai = f"{DIRS['aligned']}/{{sample}}.sorted.bam.bai",
    params:
        rg = "@RG\\tID:{sample}\\tSM:{sample}_WES\\tPL:ILLUMINA\\tLB:WES_Lib\\tPU:NovaSeq",
    log:
        "logs/{sample}/bwa_align.log"
    threads: config["threads"]["bwa"]
    shell:
        """
        mkdir -p {DIRS[aligned]}
        bwa mem -t {threads} -M \
            -R "{params.rg}" \
            {input.ref} {input.r1} {input.r2} \
            | samtools view -Sb - \
            | samtools sort -o {output.bam} 2>{log}
        samtools index {output.bam} 2>>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# RULE 6: Mark duplicates with Picard
# ─────────────────────────────────────────────────────────────────────────────
rule mark_duplicates:
    input:
        bam = rules.bwa_align.output.bam,
    output:
        bam     = f"{DIRS['aligned']}/{{sample}}.dedup.bam",
        bai     = f"{DIRS['aligned']}/{{sample}}.dedup.bam.bai",
        metrics = f"{DIRS['metrics']}/{{sample}}.dedup_metrics.txt",
    log:
        "logs/{sample}/mark_duplicates.log"
    params:
        jvm = config["memory"]["picard_jvm"],
    shell:
        """
        picard -Xmx{params.jvm} MarkDuplicates \
            I={input.bam} \
            O={output.bam} \
            M={output.metrics} \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=SILENT 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# RULE 7a: GATK BaseRecalibrator
# ─────────────────────────────────────────────────────────────────────────────
rule base_recalibrator:
    input:
        bam    = rules.mark_duplicates.output.bam,
        ref    = REF,
        dbsnp  = KNOWN["dbsnp"],
        mills  = KNOWN["mills"],
        indels = KNOWN["indels"],
    output:
        table = f"{DIRS['bqsr']}/{{sample}}.recal.table",
    log:
        "logs/{sample}/base_recalibrator.log"
    threads: config["threads"]["gatk"]
    shell:
        """
        mkdir -p {DIRS[bqsr]}
        gatk --java-options "-Xmx{config[memory][picard_jvm]}" BaseRecalibrator \
            -I {input.bam} \
            -R {input.ref} \
            --known-sites {input.dbsnp} \
            --known-sites {input.mills} \
            --known-sites {input.indels} \
            -O {output.table} 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# RULE 7b: GATK ApplyBQSR
# ─────────────────────────────────────────────────────────────────────────────
rule apply_bqsr:
    input:
        bam   = rules.mark_duplicates.output.bam,
        table = rules.base_recalibrator.output.table,
        ref   = REF,
    output:
        bam = f"{DIRS['bqsr']}/{{sample}}.bqsr.bam",
        bai = f"{DIRS['bqsr']}/{{sample}}.bqsr.bam.bai",
    log:
        "logs/{sample}/apply_bqsr.log"
    shell:
        """
        gatk --java-options "-Xmx{config[memory][picard_jvm]}" ApplyBQSR \
            -I {input.bam} \
            -R {input.ref} \
            --bqsr-recal-file {input.table} \
            -O {output.bam} 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# RULE 8a: samtools flagstat
# ─────────────────────────────────────────────────────────────────────────────
rule flagstat:
    input:
        bam = rules.apply_bqsr.output.bam,
    output:
        f"{DIRS['metrics']}/{{sample}}.flagstat.txt",
    log:
        "logs/{sample}/flagstat.log"
    shell:
        "samtools flagstat {input.bam} > {output} 2>{log}"


# ─────────────────────────────────────────────────────────────────────────────
# RULE 8b: Picard CollectHsMetrics (WES-specific on-target metrics)
# ─────────────────────────────────────────────────────────────────────────────
rule collect_hs_metrics:
    input:
        bam = rules.apply_bqsr.output.bam,
        ref = REF,
        bed = BED,
    output:
        f"{DIRS['metrics']}/{{sample}}.hs_metrics.txt",
    params:
        jvm = config["memory"]["picard_jvm"],
    log:
        "logs/{sample}/hs_metrics.log"
    shell:
        """
        picard -Xmx{params.jvm} CollectHsMetrics \
            I={input.bam} \
            O={output} \
            R={input.ref} \
            BAIT_INTERVALS={input.bed} \
            TARGET_INTERVALS={input.bed} 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# RULE 8c: samtools depth with multiple thresholds
# ─────────────────────────────────────────────────────────────────────────────
rule depth_coverage:
    input:
        bam = rules.apply_bqsr.output.bam,
        bed = BED,
    output:
        depth   = f"{DIRS['metrics']}/{{sample}}.depth.txt",
        summary = f"{DIRS['metrics']}/{{sample}}.coverage_summary.txt",
    params:
        thresholds = config["depth_thresholds"],
    log:
        "logs/{sample}/depth.log"
    shell:
        """
        samtools depth -b {input.bed} {input.bam} > {output.depth} 2>{log}
        echo "threshold\tcovered_bases\tmean_depth" > {output.summary}
        for t in {params.thresholds}; do
            cov=$(awk -v t=$t '$3>=t{{n++}} END{{print n+0}}' {output.depth})
            avg=$(awk -v t=$t '$3>=t{{s+=$3;n++}} END{{if(n>0) printf "%.2f",s/n; else print "NA"}}' {output.depth})
            echo -e "$t\\t$cov\\t$avg" >> {output.summary}
        done
        """


# ─────────────────────────────────────────────────────────────────────────────
# RULE 9: GATK HaplotypeCaller → per-sample GVCF
# ─────────────────────────────────────────────────────────────────────────────
rule haplotype_caller:
    input:
        bam = rules.apply_bqsr.output.bam,
        ref = REF,
        bed = BED,
    output:
        gvcf = f"{DIRS['variants']}/{{sample}}.g.vcf.gz",
        tbi  = f"{DIRS['variants']}/{{sample}}.g.vcf.gz.tbi",
    params:
        stand_call_conf = config["gatk"]["stand_call_conf"],
    log:
        "logs/{sample}/haplotype_caller.log"
    threads: config["threads"]["gatk"]
    shell:
        """
        mkdir -p {DIRS[variants]}
        gatk --java-options "-Xmx{config[memory][picard_jvm]}" HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -L {input.bed} \
            -O {output.gvcf} \
            -ERC GVCF \
            --stand-call-conf {params.stand_call_conf} 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# RULE 10: Genotype GVCFs (single-sample direct path)
# ─────────────────────────────────────────────────────────────────────────────
rule genotype_gvcfs:
    input:
        gvcf = rules.haplotype_caller.output.gvcf,
        ref  = REF,
    output:
        vcf = f"{DIRS['variants']}/{{sample}}.raw.vcf.gz",
        tbi = f"{DIRS['variants']}/{{sample}}.raw.vcf.gz.tbi",
    log:
        "logs/{sample}/genotype_gvcfs.log"
    shell:
        """
        gatk --java-options "-Xmx{config[memory][picard_jvm]}" GenotypeGVCFs \
            -R {input.ref} \
            -V {input.gvcf} \
            -O {output.vcf} 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# RULE 11a: Hard-filter SNPs
# ─────────────────────────────────────────────────────────────────────────────
rule filter_snps:
    input:
        vcf = rules.genotype_gvcfs.output.vcf,
        ref = REF,
    output:
        vcf = f"{DIRS['variants']}/{{sample}}.snps.filtered.vcf.gz",
    params:
        hf = config["hard_filter"]["snp"],
    log:
        "logs/{sample}/filter_snps.log"
    shell:
        """
        gatk SelectVariants -R {input.ref} -V {input.vcf} --select-type-to-include SNP \
            -O /tmp/{wildcards.sample}.snps.vcf.gz 2>{log}
        gatk VariantFiltration -R {input.ref} -V /tmp/{wildcards.sample}.snps.vcf.gz \
            --filter-expression "QD < {params.hf[QD]}"          --filter-name "QD_filter" \
            --filter-expression "MQ < {params.hf[MQ]}"          --filter-name "MQ_filter" \
            --filter-expression "FS > {params.hf[FS]}"          --filter-name "FS_filter" \
            --filter-expression "SOR > {params.hf[SOR]}"        --filter-name "SOR_filter" \
            -O {output.vcf} 2>>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# RULE 11b: Hard-filter INDELs
# ─────────────────────────────────────────────────────────────────────────────
rule filter_indels:
    input:
        vcf = rules.genotype_gvcfs.output.vcf,
        ref = REF,
    output:
        vcf = f"{DIRS['variants']}/{{sample}}.indels.filtered.vcf.gz",
    params:
        hf = config["hard_filter"]["indel"],
    log:
        "logs/{sample}/filter_indels.log"
    shell:
        """
        gatk SelectVariants -R {input.ref} -V {input.vcf} --select-type-to-include INDEL \
            -O /tmp/{wildcards.sample}.indels.vcf.gz 2>{log}
        gatk VariantFiltration -R {input.ref} -V /tmp/{wildcards.sample}.indels.vcf.gz \
            --filter-expression "QD < {params.hf[QD]}"          --filter-name "QD_filter" \
            --filter-expression "FS > {params.hf[FS]}"          --filter-name "FS_filter" \
            --filter-expression "SOR > {params.hf[SOR]}"        --filter-name "SOR_filter" \
            -O {output.vcf} 2>>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# RULE 11c: Merge filtered SNPs and INDELs
# ─────────────────────────────────────────────────────────────────────────────
rule merge_filtered_vcf:
    input:
        snps   = rules.filter_snps.output.vcf,
        indels = rules.filter_indels.output.vcf,
        ref    = REF,
    output:
        vcf = f"{DIRS['variants']}/{{sample}}.filtered.vcf.gz",
        tbi = f"{DIRS['variants']}/{{sample}}.filtered.vcf.gz.tbi",
    log:
        "logs/{sample}/merge_filtered.log"
    shell:
        """
        gatk MergeVcfs \
            I={input.snps} I={input.indels} \
            O={output.vcf} 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# RULE 12: Variant annotation with Ensembl VEP
# ─────────────────────────────────────────────────────────────────────────────
rule vep_annotate:
    input:
        vcf = rules.merge_filtered_vcf.output.vcf,
    output:
        vcf  = f"{DIRS['annotated']}/{{sample}}.vep.vcf.gz",
        html = f"{DIRS['annotated']}/{{sample}}.vep_summary.html",
    log:
        "logs/{sample}/vep.log"
    shell:
        """
        mkdir -p {DIRS[annotated]}
        vep \
            --input_file {input.vcf} \
            --output_file {output.vcf} \
            --stats_file {output.html} \
            --format vcf \
            --vcf \
            --compress_output bgzip \
            --species homo_sapiens \
            --assembly GRCh38 \
            --cache \
            --offline \
            --everything \
            --fork 4 2>{log}
        """


# ─────────────────────────────────────────────────────────────────────────────
# RULE 13: Generate MD5 checksums for key output files
# ─────────────────────────────────────────────────────────────────────────────
rule checksums:
    input:
        bam     = rules.apply_bqsr.output.bam,
        vcf     = rules.merge_filtered_vcf.output.vcf,
        metrics = rules.collect_hs_metrics.output,
    output:
        f"{DIRS['metrics']}/{{sample}}.md5sums.txt",
    shell:
        """
        md5sum {input.bam} {input.vcf} {input.metrics} > {output}
        """
