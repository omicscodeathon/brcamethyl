#!/usr/bin/env bash

SAMPLE_DIR="../pipeline/RNA-analysis/"
GENOME_DIR="/media/biodataA/brcamethyl/data/refDir/indexDir"
ALIGN_DIR="/media/biodataA/brcamethyl/data/mapping"
SAMPLE_ID="SRR26949155 SRR26949156 SRR26949157 SRR26949158 SRR26949159 SRR26949160"

for SAMPLE in $SAMPLE_ID; do
        STAR \
        --runThreadN 6 \
        --readFilesIn ${SAMPLE_DIR}/${SAMPLE}_1.fastq.gz ${SAMPLE_DIR}/${SAMPLE}_2.fastq.gz \
        --genomeDir ${GENOME_DIR} \
        --readFilesCommand zcat \
        --outSAMattrRGline ID:${SAMPLE} SM:${SAMPLE} PL:ILLUMINA \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${ALIGN_DIR}/${SAMPLE}
done
