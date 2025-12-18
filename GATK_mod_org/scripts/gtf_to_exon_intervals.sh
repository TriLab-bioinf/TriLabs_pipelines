#!/bin/bash

set -o errexit

module load bedtools/2.31.1

GTF_FILE=$1
FASTA_GENOME=$2

if [[ ! -e "$GTF_FILE" ]] || [[ ! -e "$FASTA_GENOME" ]]; then
    echo; echo "Usage: $0 <gtf_file> <fasta_genome>"; echo
    exit 1
fi

# Extract chromosome names from the FASTA genome
grep '^>' $FASTA_GENOME |sed 's/>//'|awk '{print $1}' > tmp.chroms.txt

# Convert GTF to exon intervals in BED format
awk '$3 == "exon" {print $1"\t"$4-1"\t"$5}' ${GTF_FILE} | \
    grep -F -w -f tmp.chroms.txt | \
    mergeBed -i - 

# Cleanup
rm tmp.chroms.txt
    