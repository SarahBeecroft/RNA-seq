#!/bin/bash
##wrapper around Broad's run_STAR.py script downloaded from their git repo and used to process the Gtex samples

data_dir=/data/RNAseq/new_data
ID_list=$(ls $data_dir | grep 'fastq.gz' | cut -d'_' -f1 | uniq)

for sample in $(echo $ID_list)
do

fastq_input=$(ls | grep $sample | paste -s -d, -)

python3 run_STAR.py \
        $data_dir/$fastq_input \
        $sample \
        --threads 3 \
        --output_dir $data_dir/star_out
done
