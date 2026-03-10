# BS-Seq processing pipeline paired-end workflow v1.0
# Hernan Lorenzi
# hernan.lorenzi@nih.gov
# Workflow requires to configure the config.yml file and samplesheet.csv file to include all parameters and metadatata required.

import os
import glob
import pandas as pd

# Read config file
configfile: "./config/config.yaml"
genome: str = config["reference"]["fasta"]
reads_dir: str = config["reads"]

# Read sample data from samplesheet and skip comments
metadata = pd.read_csv(config["metadata"], comment='#', sep=',', header=0, dtype=str)

# Create dictionaries from metadata for @RG line
# metadata is imported from ./config/samplesheet.csv file
fq_1: dict = {s:fq1 for s, fq1 in zip(metadata['sample_ID'], metadata['fastq_1'])}
fq_2: dict = {s:fq2 for s, fq2 in zip(metadata['sample_ID'], metadata['fastq_2'])} if config["paired"] == True else {}
samples: list = list(metadata['sample_ID'])

# Set what rules to run locally
localrules: all 

rule all:
    # IMPORTANT: output file for all rules has to match the name specified in the output file
    # and not include suffixes that the command use might add to it.
    input: "results/6-multiqc/multiqc_report.html"  

# 1- Trim and QC reads with fastp
include: "rules/trim_reads.smk"

# 2- Run bismark
include: "rules/bismark.smk"

# # 4- Merge bismark bam files
# include: "rules/merge_bam_files.smk"

# # 5- Deduplicate reads
include: "rules/deduplicate_reads.smk"

# # 6- Bismark_methylation_extractor
include: "rules/methylation_extractor.smk"

# # 7- Run multiqc
include: "rules/multiqc.smk"
