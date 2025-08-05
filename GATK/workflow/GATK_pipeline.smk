# vim: set ft=python:

# ChIP-Seq Single-end workflow v1.0
# Yuejun Wang
# yuejun.wang@nih.gov
# Workflow requires to configure the config.yml file accordingly to include all metadatata required.

import os
import glob
import pandas as pd

# Read config file
configfile: "config/config.yaml"
genome: str = config["reference"]["fasta"]
bwadb_path: str = config["bwadb"]["path"]
adapters: str = config["adapters"]

# metadata is imported from ./config/samplesheet.csv file
# Read sample data from samplesheet and skip comments
metadata = pd.read_csv(config["metadata"], comment='#', sep=',', header=0, dtype=str)
SAMPLES = metadata['sample'].unique().tolist()

def get_fastqs(sample):
    row = metadata[(metadata['sample'] == sample)]
    return [
        "data/reads/" + row.iloc[0]['fastq_1'],
        "data/reads/" + row.iloc[0]['fastq_2']
    ]


# Set what rules to run locally
localrules: all 

rule all:
    # IMPORTANT: output file for all rules has to match the name specified in the output file
    # and not include suffixes that the command use might add to it.
    input: "results/8-multiqc/multiqc_report.html"  

# 1- Trim reads
include: "rules/trim_reads.smk"

# 2- Map trimmed reads to reference with bwa-mem2
include: "rules/map_reads.smk"

# 3- Flag duplicated reads with GATK MarkDuplicates
include: "rules/mark_duplicates.smk"

# 4- 
include: "rules/apply_bqsr.smk"

# 5- 
include: "rules/haplotype_caller.smk"

# 6-
include: "rules/combine.smk"

# 7- 
include: "rules/filter_variants.smk"

# 8- Run multiqc report
include: "rules/multiqc.smk"