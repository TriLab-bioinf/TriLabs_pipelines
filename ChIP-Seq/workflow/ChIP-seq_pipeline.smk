# vim: set ft=python:

# ChIP-Seq Single-end / Paired-end workflow v2.0
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

# Create dictionaries from metadata for @RG line
# metadata is imported from ./config/samplesheet.csv file
fq_1: dict = {s:fq1 for s, fq1 in zip(metadata['sample_ID'], metadata['fastq_1'])}
fq_2: dict = {s:fq2 for s, fq2 in zip(metadata['sample_ID'], metadata['fastq_2'])} if config["paired"] == True else {}
samples: metadata['sample_ID'].unique().tolist()

def get_fastq(sample, type, read):
    """
    sample : str (sample_ID)
    type  : str (treatment/control)
    read   : str ("fastq_1" or "fastq_2")
    """
    row = metadata[(metadata['sample_ID']==sample) & (metadata['type']==type)]

    # If single-end and read=="fastq_2", return None
    if read == "fastq_2" and (not config.get("paired", False)):
        return None

    return row.iloc[0][read]

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
