#!/usr/bin/env bash

SAMPLE_DIR="../ME-pipeline/Bisulfite-analysis/bismark_files"
SAMPLES="SRR26949323 SRR26949324 SRR26949325 SRR26949326 SRR26949327 SRR26949328"

for SAMPLE in $SAMPLES; do

bismark_methylation_extractor -p \
    --gzip \
    --comprehensive \
    --cytosine_report \
    --merge_non_CpG \
    --bedGraph \
    --multicore 2 \
    --counts \
    --buffer_size 10G ${SAMPLE_DIR}/${SAMPLE}_1_val_1_bismark_bt2_pe.deduplicated.bam

done


