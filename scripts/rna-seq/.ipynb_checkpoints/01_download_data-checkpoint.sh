#!/usr/bin/env bash

SAMPLE_ID="SRR26949155 SRR26949156 SRR26949157 SRR26949158 SRR26949159 SRR26949160"

for sample in $SAMPLES;
do

fastq-dump --gzip --split-files ${sample}

done