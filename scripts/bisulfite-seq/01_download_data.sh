#!/usr/bin/env bash

#this script is used to download (convert from SRA files) fastq files in a gzipped format
#dependencies: sra-toolkit

SAMPLES=$(cat ../accessions/bisulphite-seq.txt)

for sample in $SAMPLES;
do

prefetch ${sample} --max-size u && fasterq-dump --gzip --split-files ${sample} && rm -Rv ${sample}

done
