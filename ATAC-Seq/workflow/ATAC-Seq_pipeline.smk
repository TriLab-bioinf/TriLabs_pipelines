# vim: set ft=python:

# ATAC-Seq Single-end / Paired-end workflow v2.0
# Yuejun Wang
# yuejun.wang@nih.gov
# Workflow requires to configure the config.yml file accordingly to include all metadatata required.

import os
import glob
import pandas as pd

# Read config file
configfile: "config/config.yaml"
genome: str = config["reference"]["fasta"]
gsize = config["reference"]["gsize"]
bowtiedb_path: str = config["bowtiedb"]["path"]
adapters: str = config["adapters"]
idx: str = config["reference"]["idx"]
reads_dir: str = config["reads"]

# Read sample data from samplesheet and skip comments
metadata = pd.read_csv(config["metadata"], comment='#', sep=',', header=0, dtype=str)
#SAMPLES = metadata['sample'].unique().tolist()

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

if config["trimming"] == True:
    trim_output = "results/1-trimming"
    # 1- Trim reads with fastp
    include: "rules/trim_reads.smk"
else:
    trim_output = "data/reads"

# 2- Map trimmed reads to reference with bowtie2
include: "rules/map_reads.smk"

if config["flag_dup"] == True:
    # 3- Flag duplicated reads with GATK MarkDuplicates
    include: "rules/mark_duplicates.smk"

## 4- Call peaks with macs2
include: "rules/call_peaks.smk"

# 5- Make bigwig files from bam files
include: "rules/make_bigwig.smk"

# 6- Run multiqc report
include: "rules/multiqc.smk"
