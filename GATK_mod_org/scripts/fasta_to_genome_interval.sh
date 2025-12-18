#!/bin/bash

set -o errexit

module load samtools

FASTA_FILE=$1

if [[ ! -e "$FASTA_FILE" ]]; then
    echo; echo "Usage: $0 <fasta_file>"; echo
    exit 1
fi

# generate temporary dictionary file with samtools and convert to genome intervals
samtools dict $FASTA_FILE | \
    grep '^@SQ' | \
    awk -F '[\t:]' '{print $3"\t0\t"$5}'
