# GATK_mod_org workflow v1.0
# Hernan Lorenzi
# hernan.lorenzi@nih.gov
# Workflow requires to configure the config.yml file accordingly to include all metadatata required.
# Exome Germline short variant discovery (SNPs + Indels) pipeline
# Based on GATK Best Practices:
# https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels

import os
import glob
import pandas as pd
from scripts.helper_functions import get_chr_ids,process_db_action,fetch_genomeDB_id

# Read sample data from samplesheet and skip comments
metadata = pd.read_csv('./config/samplesheet.csv', comment='#', sep=',', header=0, dtype=str)

# Read config file
configfile: "./config/config.yaml"
bwadb_path: str = config["bwadb"]
annotation: str = config["reference"]["gtf"]
adapters: str = config["adapters"]
reads_dir: str = config["reads"]
genome: str = config["reference"]["fasta"]
genome_dict: str = config["reference"]["dict"]
dbsnp: str = config["reference"]["known_sites"]
intervals: str = config["reference"]["intervals"]
whole_chr_intervals: str = config["reference"]["whole_chr_intervals"]
db_action: str = config["db_action"]

# Set what rules to run locally
localrules: all #,
            #build_abundant_db

rule all:
    # IMPORTANT: output file for all rules has to match the name specified in the output file
    # and not include suffixes that the command use might add to it.
    input: "results/14-multiqc/multiqc_report.html"  

# 1- Trim reads with fastp
include: "rules/trim_reads.smk"

# 2- Map trimmed reads to reference with bwa-mem2
include: "rules/map_reads.smk"

# 3- Flag duplicated reads with GATK MarkDuplicates
include: "rules/mark_duplicates.smk"

    # 4- Base recallibration with GATK BaseRecalibrator
include: "rules/base_recalibrator.smk"

    # 5- Call variants for single samples with GATK HaplotypeCaller (GVCF format)
include: "rules/gvcf_caller.smk"

# 6- Merge by-chr-VCFs and split variants into SNPs and INDELs
include: "rules/merge_select_variants.smk"

# 7- Filter variants with GATK VariantFiltration
include: "rules/filter_variants.smk"

# 8- Run multiqc report
include: "rules/multiqc.smk"