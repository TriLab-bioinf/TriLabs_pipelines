# TriLabs Snakemake Pipelines

This repository contains a collection of Snakemake pipelines for processing various types of sequencing data. All pipelines are optimized for running on the Biowulf HPC cluster at NIH.

## Available Pipelines

### 1. RNA-Seq Pipeline

The RNA-Seq pipeline is designed for processing bulk RNA sequencing data. It performs comprehensive analysis from raw sequencing reads to gene expression quantification and differential expression analysis.

**Key Features:**

- Handling of single- or paired-end reads
- Quality control of raw reads (fastp)
- Read trimming and adapter removal (fastp)
- Alignment to reference genome (STAR)
- Duplicate removal (GATK's Picard)
- Gene expression quantification (featureCounts)
- Generation of expression profiles in BigWig format (bamCoverage)
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

- Quality control and read trimming
- Alignment to reference genome
- Duplicate removal
- Peak calling for transcription factor binding sites or histone marks
- Peak visualization


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

- Quality control and adapter trimming
- Alignment to reference genome
- Duplicate removal
- Peak calling for open chromatin regions


**Typical Use Cases:**

- Mapping chromatin accessibility landscapes
- Identification of regulatory elements
- Transcription factor footprinting
- Nucleosome positioning analysis

**Location:** `ATAC-Seq/`

**Documentation:** See [ATAC-Seq/README.md](ATAC-Seq/README.md)

---

### 4. GATK Pipeline (Human)

The GATK pipeline is specifically designed for processing human whole genome sequencing data following GATK best practices for variant calling. It identifies single nucleotide polymorphisms (SNPs) and insertions/deletions (INDELs) in human samples.

**Key Features:**

- Read alignment and preprocessing
- Duplicate removal
- Base quality score recalibration (BQSR)
- Variant calling using HaplotypeCaller
- Variant quality score recalibration (VQSR)
- Joint genotyping across multiple samples
- Variant filtering
- Support for WGS data

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
- Read alignment and preprocessing
- Duplicate removal
- Variant calling with HaplotypeCaller
- Hard filtering for variants (alternative to VQSR for non-model organisms)
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
2. **Dry Run:** Test pipeline with `snakemake -n -p`
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
