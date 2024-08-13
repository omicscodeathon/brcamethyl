#!/usr/bin/env bash

SAMPLE_DIR="../ME-pipeline/Bisulfite-analysis/bismark_files"
SAMPLES="SRR26949323 SRR26949324 SRR26949325 SRR26949326 SRR26949327 SRR26949328"

for SAMPLE in $SAMPLES; do

deduplicate_bismark -p \
    --output_dir ${SAMPLE_DIR} \
    --bam \
    ${SAMPLE_DIR}/${SAMPLE}_1_val_1_bismark_bt2_pe.bam

done