# TriLabs Snakemake Pipelines

This repository contains a collection of Snakemake pipelines for processing various types of sequencing data. All pipelines are optimized for running on the Biowulf HPC cluster at NIH.

## Installation

Note: This should be done only once.

### 1- Start an interactive session in Biowulf and go to your working directory (WD)

### 2- Download the TriLabs Snakemake Pipelines to your working directory by running the following command

```bash
git clone https://github.com/TriLab-bioinf/TriLabs_pipelines.git && cd TriLabs_pipelines
```

#### Now you should be ready to run the pipelines.

## Available Pipelines

### 1. RNA-Seq Pipeline

The RNA-Seq pipeline is designed for processing bulk RNA sequencing data. It performs comprehensive analysis from raw sequencing reads to gene expression quantification and differential expression analysis.

**Key Features:**

- Handling of single- or paired-end reads
- Compatible with bacterial
- Quality control of raw reads (fastp)
- Read trimming and adapter removal (fastp)
- Alignment to reference genome (STAR)
- Duplicate removal (GATK's Picard)
- Gene expression quantification (featureCounts)
- Generation of expression profiles in BigWig format (deeptools::bamCoverage)
- Pipeline QC report (MultiQC)

**Typical Use Cases:**

- Differential gene expression analysis
- Transcriptome profiling
- Gene expression quantification across samples

**Location:** `RNA-Seq/`

**Documentation:** See [RNA-Seq/README.md](RNA-Seq/README.md)

---

### 2. ChIP-Seq Pipeline

The ChIP-Seq pipeline processes chromatin immunoprecipitation sequencing data to identify protein-DNA binding sites and histone modification patterns across the genome.

**Key Features:**

- Quality control and read trimming (fastp)
- Alignment to reference genome (bowtie2)
- Duplicate removal (Picard)
- Narrow and broad peak calling (macs2)
- Generation of peak profiles in BigWig format (deeptools::bamCoverage)
- Pipeline QC report (MultiQC)

**Typical Use Cases:**

- Identification of transcription factor binding sites
- Mapping histone modification landscapes
- Comparative analysis of protein-DNA interactions

**Location:** `ChIP-Seq/`

**Documentation:** See [ChIP-Seq/README.md](ChIP-Seq/README.md)

---

### 3. ATAC-Seq Pipeline

The ATAC-Seq pipeline analyzes Assay for Transposase-Accessible Chromatin sequencing data to map open chromatin regions and nucleosome positioning across the genome.

**Key Features:**

- Quality control and read trimming (fastp)
- Alignment to reference genome (bowtie2)
- Duplicate removal (Picard)
- Peak calling for open chromatin regions (macs2)
- Generation of peak profiles in BigWig format (deeptools::bamCoverage)
- Pipeline QC report (MultiQC)

**Typical Use Cases:**

- Mapping chromatin accessibility landscapes
- Identification of regulatory elements
- Transcription factor footprinting
- Nucleosome positioning analysis

**Location:** `ATAC-Seq/`

**Documentation:** See [ATAC-Seq/README.md](ATAC-Seq/README.md)

---

### 4. GATK Pipeline (Human germline)

The GATK pipeline is specifically designed for processing human whole genome sequencing data following GATK best practices for variant calling. It identifies single nucleotide polymorphisms (SNPs) and insertions/deletions (INDELs) in human samples.

**Key Features:**

- Quality control and read trimming (fastp)
- Read alignment and preprocessing (bwa2)
- Duplicate removal (GATK MarkDuplicates)
- Base quality score recalibration (BQSR) (GATK BaseRecalibrator -> ApplyBQSR)
- Variant calling per sample using HaplotypeCaller (GVCF) (GATK HaplotypeCaller)
- Combine per-sample GVCF files (GATK CombineGVCFs)
- Joint genotyping across multiple samples (GATK GenotypeGVCFs)
- Variant quality score recalibration and filtering (VQSR) (VariantRecalibrator -> ApplyVQSR)
- Support for WGS data

_Note: Support for WES data is not yet implemented_

**Typical Use Cases:**

- Human germline variant discovery
- Population genetics studies
- Clinical sequencing analysis
- Identification of disease-associated variants

**Location:** `GATK/`

**Documentation:** See [GATK/README.md](GATK/README.md)

---

### 5. GATK Modified for Other Organisms (GATK_mod_org)

The GATK_mod_org pipeline is a flexible, organism-agnostic variant calling pipeline adapted from GATK best practices. It can process whole genome or whole exome sequencing data from any organism for SNP and INDEL identification.

**Key Features:**

- Organism-independent workflow
- Configurable reference genome and annotations
- Quality control and read trimming (fastp)
- Read alignment to reference (bwa2)
- Duplicate removal (GATK MarkDuplicates)
- Base quality score recalibration (BQSR) (GATK BaseRecalibrator -> ApplyBQSR)
- Variant calling per sample using HaplotypeCaller (GVCF) (GATK HaplotypeCaller)
- Implementation of GenomicsDB for efficient variant analysis of large cohorts (GATK GenomicsDBImport)
- Joint genotyping across multiple samples (GATK GenotypeGVCFs)
- Hard filtering for variants (alternative to VQSR for non-model organisms) (GATK VariantFiltration)
- Support for both WGS and WES data
- Custom interval targeting for exome data
- Parallel calling of variants per chromosome/contig

**Typical Use Cases:**

- Variant discovery in non-human organisms

**Location:** `GATK_mod_org/`

**Documentation:** See [GATK_mod_org/README.md](GATK_mod_org/README.md)

---

## General Pipeline Usage

All pipelines follow a similar structure and execution pattern:

1. **Configuration:** Edit `config/config.yaml` and `config/samplesheet.csv`
2. **Dry Run:** Test pipeline with `run_snakemake.sh -n`
3. **Execution:** Submit to cluster with `sbatch run_snakemake.sh`
4. **Results:** Find outputs in `results/` directory

## Command Extraction Utility

All pipelines include a Python utility (`parse_snakemake_commands.py`) that extracts and organizes shell commands from Snakemake dry-run output. This is useful for:

- Understanding pipeline execution flow
- Debugging specific pipeline steps
- Running individual commands manually
- Documentation and reproducibility

**Usage:**

```bash
snakemake --snakefile ./workflow/PIPELINE.smk -p -n --forceall > snakemake_output.txt
python parse_snakemake_commands.py -i snakemake_output.txt -o commands_ordered.txt
```

## Requirements

- Access to Biowulf HPC cluster (or similar SLURM-based system)
- Snakemake 7.32.4 or higher
- Appropriate reference genomes and annotations
- Required software modules (loaded via environment modules)

## Support

For questions or issues with specific pipelines, please refer to the individual README files in each pipeline directory or contact the TriLab Bioinformatics team.
