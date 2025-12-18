# How to run the GATK pipeline 

![](dag.png)

**Figure 1:** Schematic representation of the Snakemake GATK pipeline

## A. Setup the pipeline for the first time:

Note: This should be done only once.

### 1- Start an interactive session in Biowulf and go to your working directory (WD).

### 2- Download the GATK pipeline by running the  following command in your WD:

```bash
git clone https://github.com/TriLab-bioinf/TriLabs_pipelines.git && cd TriLabs_pipelines/GATK
```

### 3- Download the biowulf snakemake profile from GitHub to GATK/config directory

```bash
git clone https://github.com/NIH-HPC/snakemake_profile.git ./config/snakemake_profile
```

### Now you should be ready to run the GATK pipeline 

## B. Running the GATK pipeline in Biowulf

### 1- Within the config directory do the following:

- Edit the [config.yaml](config/config.yaml) file with required information
- Edit [samplesheet.csv](config/samplesheet.csv) with sample data information

### 2- Whole Exome Sequencing (WES) data

If you are working with Whole Exome Sequencing (WES) data, you will need to provide the following two files:

- _chr_intervals.bed_: a bed file with the start and end coordinates of each chromosome or contig in the reference fasta file used for mapping. 
- _merged.exons.bed_: a bed file with the coordinates of all exons used for generating the WES data. If this information is unknown, you can generate this file from the annotation file in gtf/gff format.

There are two bash scripts that can be used to generate the two bed files above. To create _chr_intervals.bed_ you can run the following commands from the terminal:

```bash
# Using the fasta file for your reference genome run
# the following from the pipeline's working directory:

FASTA_GENOME=/path/to/fasta/file

./scripts/fasta_to_genome_interval.sh ${FASTA_GENOME} > ./data/intervals/chr_intervals.bed
```

_Note_: The script _fasta_to_genome_interval.sh_ makes use of the _samtools_ tool via the `load module` command, so make sure _samtools_ can be found via your environmental variable $PATH.

To generate the _merged.exons.bed_ file run the following commands in a terminal:

```bash
# Using the gtf/gff annotation file for your reference genome run
# the following from the pipeline's working directory:
GTF_FILE=/path/to/gtf_gff/file
FASTA_GENOME=/path/to/fasta/file

./scripts/gtf_to_exon_intervals.sh ${GTF_FILE} ${FASTA_GENOME} > ./data/intervals/merged.exons.bed
```

_Note_: The script _gtf_to_exon_intervals.sh_ makes use of the _bedtools_ tool via the `load module` command, so make sure _bedtools_ can be found via your environmental variable $PATH.

### 3- Load Snakemake module

```bash
module load snakemake/7.32.4
```

### 4- OPTIONAL: Activate conda environment (if running snakemake using 5.a or 5.b below)

```bash
source ~/bin/myconda
```

### 5.a- To run the Snakemake pipeline to process sequencing data locally (dry-run)

```bash
snakemake --profile ./config/snakemake_profile --snakefile ./workflow/GATK_pipeline.smk -p -n
```

### 5.b- To run the Snakemake pipeline to process sequencing data locally

```bash
snakemake --profile ./config/snakemake_profile --snakefile ./workflow/GATK_pipeline.smk -p
```

### 5.c- To run the Snakemake pipeline to process sequencing data in a cluster machine (best option)

```bash
sbatch run_snakemake.sh ./workflow/GATK_pipeline.smk
```

### 6- To fetch the shell commands ran by the snakemake pipeline use the following command

```bash
bash print_snakemake_shell_commands.sh -c module,fastp,bwa-mem2,samtools,picard,GATK,multiqc -p workflow/GATK_pipeline.smk
```

The command above will create the file [snakemake_shell_commands.txt](workflow/snakemake_shell_commands.txt) within the `workflow` directory with an unsorted print out of all shell commands executed by the snakemake pipeline.
