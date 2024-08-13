#!/usr/bin/env bash

#this script performs quality check on SRA files using FastQC Tool
cd ~/brcamethyl/dataset/bisulfite-seq 

# fastq files directory
FASTQ_DIR="~/brcamethyl/dataset/bisulfite-seq"

for file in $FASTQ_DIR/*.fastq; do
        mkdir qc_reports && \
        fastqc ${file} -o qc_reports
done

multiqc .
