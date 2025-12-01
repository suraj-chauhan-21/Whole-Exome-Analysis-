# Whole-Exome-Analysis-
The human genome consists of approximately 3 billion base pairs, but only ~1-2% of this DNA codes for proteins. This coding region is known as the Exome.  Whole Exome Sequencing is a targeted sequencing strategy that focuses exclusively on these protein-coding regions.


1. The Genome vs. The Exome

    The Genome: The entire human DNA sequence (approx. 3 billion base pairs). It contains everything: genes, regulatory sequences, and a vast amount of non-coding DNA (often called "junk DNA," though it does have functions).

    The Exome: The specific part of the genome formed by Exons (the coding regions). These are the parts of the DNA that actually get translated into proteins.

   <img width="532" height="673" alt="image" src="https://github.com/user-attachments/assets/d3fadd7d-03c5-4795-8060-a8e504dd4a43" />


The exome makes up only about 1% to 2% of the whole genome.


Why do we need it? (The Advantages) If we can sequence the whole genome (WGS), why bother with just the exome?


**Clinical Relevance**: About 85% of known disease-causing mutations occur in the exome. If you are looking for a genetic cause of a disease (like cancer or a hereditary disorder), the answer is almost certainly in the exome.

**Cost Efficiency**: Sequencing 1% of the genome is significantly cheaper than sequencing 100% of it.

**Higher Depth (Accuracy):** Because you are sequencing a smaller target, you can sequence it "deeper" (more times) for the same cost.

Example: For $500, you might get 30× coverage on a Whole Genome, or 100× coverage on an Exome. Higher coverage means you can be much more confident that a mutation you found is real and not a machine error.

**Data Management**: WGS files are massive (hundreds of Gigabytes). WES files are smaller and easier to store and analyze on standard computers.






# **Whole Exome Sequencing (WES) Analysis – SRR22317682**


Short, clean workflow containing only essential commands + brief purpose notes.

*Note: 

“10× depth” below is only an example. Real QC usually evaluates **10×, 20×, 30×, 50×, 100×**, etc.*

> Due to the large size of sequencing and alignment files (FASTQ, BAM, reference genomes), only the pipeline structure, code, and essential README documentation are uploaded here.  

> The full dataset and output files are stored locally and can be provided upon request.

---

================================================================================
                  WHOLE EXOME SEQUENCING (WES) PIPELINE FLOWCHART
================================================================================

   [ INPUTS ]                                          [ PROCESS ]
       |
       |--------------------------------------------------+
       |                                                  |
  (Raw FastQ Files)                               (Reference Genome)
 SRR22317682_1.fastq                                   hg38.fa
 SRR22317682_2.fastq                                      |
       |                                                  v
       |                                             [ Indexing ]
       |                                           (bwa index/faidx)
       |                                                  |
       |                                                  |
       v                                                  v
  +--------------------------------------------------------------------------+
  |  1. ALIGNMENT (BWA-MEM)                                                  |
  |     Map raw reads to reference genome.                                   |
  +--------------------------------------------------------------------------+
       |
       v
  (Raw SAM/BAM)
       |
       v
  +--------------------------------------------------------------------------+
  |  2. SORTING & INDEXING (Samtools)                                        |
  |     Organize reads by genomic coordinate (chr1, chr2...).                |
  +--------------------------------------------------------------------------+
       |
       v
  (Sorted BAM)
       |
       v
  +--------------------------------------------------------------------------+
  |  3. MARK DUPLICATES (Picard)                                             |
  |     Flag PCR duplicates (optical/chemical artifacts).                    |
  +--------------------------------------------------------------------------+
       |                                          |
       v                                          |---> [ QC METRICS ]
  (Dedup BAM)                                     |     - Insert Size
       |                                          |     - Library Complexity (ROI)
       |                                          |     - Depth of Coverage
       v
  +--------------------------------------------------------------------------+
  |  4. BQSR (GATK BaseRecalibrator)                                         |
  |     Fix sequencing machine errors using known sites (dbSNP).             |
  +--------------------------------------------------------------------------+
       |
       v
  (Recalibrated BAM) <--- *FINAL CLEAN BAM*
       |
       v
  +--------------------------------------------------------------------------+
  |  5. VARIANT CALLING (GATK HaplotypeCaller)                               |
  |     Identify SNPs and Indels in Exome regions (BED file).                |
  +--------------------------------------------------------------------------+
       |
       v
  (Raw VCF)
       |
       v
  +--------------------------------------------------------------------------+
  |  6. HARD FILTERING (GATK VariantFiltration)                              |
  |     Remove low-confidence variants (Quality < 30, Strand Bias, etc.).    |
  +--------------------------------------------------------------------------+
       |
       v
  (Filtered VCF)
       |
       v
  +--------------------------------------------------------------------------+
  |  7. ANNOTATION (SnpEff)                                                  |
  |     Add biological context (Gene name, Missense/Silent, Impact).         |
  +--------------------------------------------------------------------------+
       |
       v
================================================================================
                              FINAL OUTPUT
                     SRR22317682.annotated.vcf.gz
================================================================================



## **1. Download SRA Data**


> **Purpose:** Obtain raw FASTQ reads for WES sample SRR22317682.


```bash

apt install sra-toolkit


mkdir -p data && cd data

prefetch -O ./ --progress SRR22317682

fasterq-dump SRR22317682.sra --split-files --threads 4

```


---


## **2. Download Reference Genome (GRCh38)**


> **Purpose:** Required for alignment and downstream analysis.


```bash

conda create -n ncbi_datasets

conda activate ncbi_datasets

conda install -c conda-forge ncbi-datasets-cli


mkdir HG38 && cd HG38

datasets download genome accession GCF_000001405.40 --include genome,gff3 --filename human_GRCh38.zip

unzip human_GRCh38.zip

```


---


## **3. Index Reference**


> **Purpose:** Enable alignment and fast genomic coordinate lookup.


```bash

conda activate ngs1


bwa index ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna

samtools faidx ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna

```


---


## **4. Alignment with BWA-MEM**


> **Purpose:** Map reads to the human genome and produce sorted BAM.


```bash

REF="HG38/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna"

RUN="SRR22317682"

R1="data/${RUN}/${RUN}_1.fastq"

R2="data/${RUN}/${RUN}_2.fastq"

mkdir -p aligned

OUT="aligned/${RUN}.sorted.bam"


RG="@RG\tID:${RUN}\tSM:${RUN}_WES\tPL:ILLUMINA\tLB:WES_Lib\tPU:NovaSeq"


bwa mem -t 8 -M -R "$RG" "$REF" "$R1" "$R2" \

  | samtools view -Sb - \

  | samtools sort -o "$OUT"


samtools index "$OUT"

```


---


## **5. Coverage Calculation Using BED File**


> **Why BED file is used?**

> A BED file defines **only the target exome regions** (Agilent V6).

> Coverage should be calculated *only* in these captured regions, not the whole genome.


### **Fix BED chromosome names (if needed)**


```bash

awk 'NR==FNR{a[$1]=$2;next}{if($1 in a)$1=a[$1];print}' \

  HG38/chr_map.txt \

  HG38/Exome-Agilent_V6.bed \

  > HG38/Exome-Agilent_V6_fixed.bed

```


### **Calculate depth**


```bash

mkdir -p metrics


samtools depth -b HG38/Exome-Agilent_V6_fixed.bed \

  aligned/SRR22317682.sorted.bam > metrics/depth.txt


# Covered bases (≥1×)

awk '$3>=1' metrics/depth.txt | wc -l


# Average depth

awk '$3>=1 {sum+=$3; n++} END {print sum/n}' metrics/depth.txt

```


---


## **6. Downsample to ~10× (Example Only)**


> **Purpose:** Demonstrate how depth affects QC (10× is only an educational example).


```bash

samtools view -s 123.0813 -b aligned/SRR22317682.sorted.bam \

  -o aligned/SRR22317682.subsampled.bam

samtools index aligned/SRR22317682.subsampled.bam


samtools depth -b HG38/Exome-Agilent_V6_fixed.bed \

  aligned/SRR22317682.subsampled.bam > metrics/depth_subsampled.txt

```


---


## **7. Mark Duplicates (Picard)**


> **Purpose:** Remove PCR duplicates for accurate coverage & variant calling.


```bash

conda install -c bioconda picard java-jdk=8.0.112


picard MarkDuplicates \

  I=aligned/SRR22317682.sorted.bam \

  O=aligned/SRR22317682.dedup.bam \

  M=metrics/SRR22317682.dedup_metrics.txt \

  CREATE_INDEX=true

```


### **Subsampled BAM**


```bash

picard MarkDuplicates \

  I=aligned/SRR22317682.subsampled.bam \

  O=aligned/SRR22317682.sub.dedup.bam \

  M=metrics/SRR22317682.sub.dedup_metrics.txt \

  CREATE_INDEX=true

```


---


## **8. Insert Size Metrics**


> **Purpose:** Evaluate library preparation and fragment size distribution.


```bash

picard CollectInsertSizeMetrics \

  I=aligned/SRR22317682.dedup.bam \

  O=metrics/SRR22317682.insert_metrics.txt \

  H=metrics/SRR22317682.insert_hist.pdf

```


---


## **9. Library Complexity**


> **Purpose:** Estimate unique molecules vs duplicates.


```bash

picard EstimateLibraryComplexity \

  I=aligned/SRR22317682.dedup.bam \

  O=metrics/SRR22317682.library_complexity.txt \

  TMP_DIR=/tmp

```


---


## **10. ROI (Return on Investment)**


> **Purpose:** How many reads remain after removing duplicates.


```bash

awk 'NR==8 {print ($3-$7)/$3}' metrics/SRR22317682.library_complexity.txt

```

## 📌 Coverage & Depth Summary


| Sample         | Covered Bases (≥1×) | Average Depth |

| -------------- | ------------------- | ------------- |

| **Original**   | 60,197,453          | 123×          |

| **Subsampled** | 58,035,334          | 10.38×        |


**Interpretation:**

Subsample retains most target regions; depth matches expected coverage reduction.


---


## 📌 Library Complexity & ROI Summary


| Sample         | Read Pairs Examined | Duplicates | % Duplication | Estimated Library Size | ROI    |

| -------------- | ------------------- | ---------- | ------------- | ---------------------- | ------ |

| **Original**   | 52,928,161          | 295,987    | 0.5592%       | 4,714,626,709          | 0.994  |

| **Subsampled** | 4,305,140           | 2,059      | 0.0478%       | 4,499,349,301          | 0.9995 |


**Interpretation:**


* Both samples show *extremely low duplication* → high library complexity.

* Subsampling preserved the diversity and complexity of the original library.


---


## 📌 Original Picard Complexity Output 


### **Original Sample**


```

EstimatedLibrarySize = 4714626709

ReadPairsExamined    = 52928161

ReadPairDuplicates   = 295987

PercentDuplication   = 0.005592

```


### **Subsampled Sample**


```

EstimatedLibrarySize = 4499349301

ReadPairsExamined    = 4305140

ReadPairDuplicates   = 2059

PercentDuplication   = 0.000478

```

## **Summary**


* Alignment, coverage, and duplicate processing follow standard WES best practices.

* BED file ensures evaluation strictly on targeted exome regions.

* 10× subsampling is only an example, not a QC standard.

* Real WES quality checks commonly include **10×, 20×, 30×, 50×, 100×** depth summaries.

I am planning to add this in on my Github can you first explain the what is Whole exome analysis, why even we need this 



