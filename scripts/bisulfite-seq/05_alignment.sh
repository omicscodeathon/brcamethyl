#!/usr/bin/env bash

TRIM_DIR="/home/jupyter/methylb/bisulfite-seq/trimbseq"
GENOME_DIR="/home/jupyter/methylb/data/refDir"
OUT_DIR="$/home/jupyter/methylb/bisulfite-seq/aligned"
SAMPLE="SRR26949323"

for sample in "${SAMPLES}"; do

bismark \
        --parallel 2 \
        --genome "${GENOME_DIR}" \
        -1 "${TRIM_DIR}"/"${sample}"_1_val_1.fq.gz \
        -2 "${TRIM_DIR}"/"${sample}"_2_val_2.fq.gz
        -o "${OUT_DIR}"
