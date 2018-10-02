#!/bin/bash
base_dir=/data/RNAseq/Broad_gtex
data_dir=$base_dir/data

set -x
 echo "$1 file"
  python run_RSEM.py \
        $base_dir/rsem_reference \
        $data_dir/star_out/$1.Aligned.toTranscriptome.out.bam \
        $data_dir/$1 \
        --threads 4