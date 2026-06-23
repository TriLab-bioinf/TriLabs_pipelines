#!/usr/bin/env bash
# Download fastq files from SRA using fasterq-dump
# Usage: download_fastq_from_SRA.sh <SRA_accession> <output_directory
set -o errexit

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <SRA_accession> <output_directory>"
    exit 1
fi

SRA_ACCESSION=$1
OUTPUT_DIR=$2

module load sratoolkit/3.3.0

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"
# Download fastq files using fasterq-dump
echo "Downloading fastq files for ${SRA_ACCESSION} to ${OUTPUT_DIR}..."

fasterq-dump ${SRA_ACCESSION} -O ${OUTPUT_DIR} --split-files --threads 6 --progress 

# Compress output fastq files
for file in ${OUTPUT_DIR}/${SRA_ACCESSION}*.fastq; do
    echo "Compressing fastq file ${file}..."
    gzip "$file"
done

# rename files to match expected naming convention
echo "Renaming files for ${SRA_ACCESSION}..."
for file in ${OUTPUT_DIR}/${SRA_ACCESSION}*.fastq.gz; do
    if [[ $file == *"${SRA_ACCESSION}_1.fastq.gz" ]]; then
        mv "$file" "${file/${SRA_ACCESSION}_1.fastq.gz/${OUTPUT_DIR}.R1.fastq.gz}"
    elif [[ $file == *"${SRA_ACCESSION}_2.fastq.gz" ]]; then
        mv "$file" "${file/${SRA_ACCESSION}_2.fastq.gz/${OUTPUT_DIR}.R2.fastq.gz}"
    fi
done

echo "Download completed for ${SRA_ACCESSION}. Fastq files are saved in ${OUTPUT_DIR}."
