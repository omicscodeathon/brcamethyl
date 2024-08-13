#!/usr/bin/env bash

SAMPLES="SRR26949323 SRR26949324 SRR26949325 SRR26949326 SRR26949327 SRR26949328"

for sample in $SAMPLES;
do

fasterq-dump --split-files ${sample}  

done
