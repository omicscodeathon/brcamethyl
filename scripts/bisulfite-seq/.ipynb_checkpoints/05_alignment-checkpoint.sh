#!/usr/bin/env bash

TRIM_DIR="../ME-pipeline/Bisulfite-analysis/trimmed_files"
GENOME_DIR="../ME-pipeline/Bisulfite-analysis/ref_genome"
OUT_DIR="../ME-pipeline/Bisulfite-analysis/bismark_files"
SAMPLES="SRR26949323 SRR26949324 SRR26949325 SRR26949326 SRR26949327 SRR26949328"

for sample in ${SAMPLES}; do

bismark \
    --parallel 2 \
    --genome "${GENOME_DIR}" \
    -1 "${TRIM_DIR}"/"${sample}"_1_val_1.fq.gz \
    -2 "${TRIM_DIR}"/"${sample}"_2_val_2.fq.gz \
    -o "${OUT_DIR}"

done
