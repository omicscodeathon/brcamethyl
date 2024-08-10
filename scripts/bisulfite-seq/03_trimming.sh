#!/usr/bin/bash

FASTQ_DIR="~/data/bisulphite-seq"

SAMPLES="$(cat ../accessions/bisulphite.txt)"

for SAMPLE in $SAMPLES; do

trim_galore --paired \
--cores 2 \
--phred33 \
--fastqc \
${FASTQ_DIR}/${SAMPLE}_1.fastq.gz \
${FASTQ_DIR}/${SAMPLE}_2.fastq.gz

done
